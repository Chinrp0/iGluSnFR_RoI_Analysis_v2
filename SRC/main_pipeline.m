function results = main_pipeline(folder, options)
    % MAIN_PIPELINE - Orchestrates fluorescent imaging analysis
    % Implements Option 1: Process per file, then combine
    %
    % Inputs:
    %   folder  - Path to folder containing CSV files
    %   options - Configuration struct (optional)
    %             .verbose (logical, default true)
    %             .useParallel (logical, default true)
    %             .continueOnError (logical, default true)
    %
    % Outputs:
    %   results - Struct containing analysis results
    %             .fileResults - Cell array of per-file results
    %             .metadata - Array of file metadata
    %             .errors - Information about failed files
    %             .summary - Overall processing statistics
    
    % Handle optional arguments
    if nargin < 2
        options = struct();
    end
    
    % Set defaults
    if ~isfield(options, 'verbose'), options.verbose = true; end
    if ~isfield(options, 'useParallel'), options.useParallel = true; end
    if ~isfield(options, 'continueOnError'), options.continueOnError = true; end
    
    if options.verbose
        fprintf('=== Fluorescent Imaging Analysis Pipeline ===\n');
        fprintf('Processing folder: %s\n', folder);
    end
    
    % Initialize results structure
    results = struct();
    results.fileResults = {};
    results.metadata = [];
    results.errors = struct('failedFiles', {}, 'errorMessages', {});
    results.summary = struct();
    
    try
        %% Step 1: Load CSV Data (Parallel)
        if options.verbose
            fprintf('\nStep 1: Loading CSV files...\n');
        end
        
        tic;
        [dataCell, metadataArray] = loadDataFiles(folder, options);
        loadTime = toc;
        
        if isempty(dataCell)
            error('No valid CSV files found or all files failed to load');
        end
        
        if options.verbose
            fprintf('  Loaded %d files in %.3f s\n', length(dataCell), loadTime);
        end
        
        %% Step 2: Process Each File Individually (Option 1 Strategy)
        if options.verbose
            fprintf('\nStep 2: Processing files individually...\n');
        end
        
        tic;
        [fileResults, processingErrors] = processFilesIndividually(dataCell, metadataArray, options);
        processTime = toc;
        
        if options.verbose
            fprintf('  Processed %d files in %.3f s\n', length(fileResults), processTime);
        end
        
        %% Step 3: Compile Results
        results.fileResults = fileResults;
        results.metadata = metadataArray;
        results.errors = processingErrors;
        
        % Create summary statistics
        results.summary = createSummaryStats(fileResults, metadataArray, loadTime, processTime);
        
        if options.verbose
            printSummary(results.summary, results.errors);
        end
        
    catch ME
        if options.verbose
            fprintf('ERROR in main pipeline: %s\n', ME.message);
        end
        results.errors.pipelineError = ME.message;
        rethrow(ME);
    end
end

function [dataCell, metadataArray] = loadDataFiles(folder, options)
    % Load CSV files using optimized loader
    
    try
        % Use the proven csv_loader_v2
        loader = csv_loader_v2();
        
        % Configure loading options
        loadOptions = struct();
        loadOptions.useParallel = options.useParallel;
        loadOptions.dataType = 'single';  % Memory efficient
        loadOptions.expectedFrames = 1200;
        
        [dataCell, metadataArray] = loader.loadBatchFiles(folder, loadOptions);
        
    catch ME
        warning('Failed to load data files: %s', ME.message);
        dataCell = {};
        metadataArray = [];
    end
end

function [fileResults, errors] = processFilesIndividually(dataCell, metadataArray, options)
    % Process each file individually (Option 1 strategy)
    % This is where baseline calculation and dF/F will be implemented
    
    numFiles = length(dataCell);
    fileResults = cell(numFiles, 1);
    
    % Track errors
    errors = struct();
    errors.failedFiles = {};
    errors.errorMessages = {};
    
    if options.verbose
        fprintf('  Processing %d files individually...\n', numFiles);
    end
    
    for i = 1:numFiles
        try
            % Get file data and metadata
            data = dataCell{i};          % [1200 x N_ROIs] matrix
            metadata = metadataArray(i); % File metadata
            
            if options.verbose && mod(i, 2) == 0  % Progress every 2 files
                fprintf('    Processing file %d/%d (%s)\n', i, numFiles, metadata.filename);
            end
            
            % Process this file
            fileResult = processSingleFile(data, metadata, options);
            fileResults{i} = fileResult;
            
        catch ME
            % Handle individual file errors
            if options.continueOnError
                if options.verbose
                    fprintf('    WARNING: Failed to process file %d (%s): %s\n', ...
                        i, metadataArray(i).filename, ME.message);
                end
                
                errors.failedFiles{end+1} = metadataArray(i).filename;
                errors.errorMessages{end+1} = ME.message;
                fileResults{i} = [];
            else
                % Re-throw if not continuing on error
                rethrow(ME);
            end
        end
    end
    
    % Remove empty results from failed files
    validResults = ~cellfun(@isempty, fileResults);
    fileResults = fileResults(validResults);
end

function fileResult = processSingleFile(data, metadata, options)
    % Process a single file's data matrix
    % Input: data is [1200 x N_ROIs] matrix for this file
    %
    % TODO: This is where we'll implement:
    %   1. Baseline calculation
    %   2. dF/F computation  
    %   3. ROI-specific processing
    %
    % For now, return basic file information
    
    [numFrames, numROIs] = size(data);
    
    % Placeholder processing - will be replaced with baseline calculation
    fileResult = struct();
    fileResult.filename = metadata.filename;
    fileResult.numFrames = numFrames;
    fileResult.numROIs = numROIs;
    fileResult.dataSize_MB = numel(data) * 4 / (1024^2);  % Single precision
    
    % TODO: Add baseline calculation here
    % fileResult.baseline = calculate_baseline(data);
    % fileResult.dfof = calculate_dfof(data, baseline);
    
    if options.verbose && numROIs > 2000
        fprintf('      Large file: %d ROIs, %.1f MB\n', numROIs, fileResult.dataSize_MB);
    end
end

function summary = createSummaryStats(fileResults, metadataArray, loadTime, processTime)
    % Create overall processing summary
    
    if isempty(fileResults)
        summary = struct('totalFiles', 0, 'validFiles', 0);
        return;
    end
    
    % Extract statistics from valid results
    numROIs = arrayfun(@(r) r.numROIs, [fileResults{:}]);
    dataSizes = arrayfun(@(r) r.dataSize_MB, [fileResults{:}]);
    
    summary = struct();
    summary.totalFiles = length(metadataArray);
    summary.validFiles = length(fileResults);
    summary.failedFiles = summary.totalFiles - summary.validFiles;
    
    summary.timing = struct();
    summary.timing.loadTime_s = loadTime;
    summary.timing.processTime_s = processTime;
    summary.timing.totalTime_s = loadTime + processTime;
    
    summary.data = struct();
    summary.data.totalROIs = sum(numROIs);
    summary.data.totalSize_MB = sum(dataSizes);
    summary.data.avgROIsPerFile = round(mean(numROIs));
    summary.data.roiRange = [min(numROIs), max(numROIs)];
end

function printSummary(summary, errors)
    % Print formatted processing summary
    
    fprintf('\n=== Pipeline Processing Summary ===\n');
    fprintf('Files: %d/%d processed successfully\n', summary.validFiles, summary.totalFiles);
    
    if summary.failedFiles > 0
        fprintf('FAILED: %d files\n', summary.failedFiles);
        for i = 1:length(errors.failedFiles)
            fprintf('  - %s: %s\n', errors.failedFiles{i}, errors.errorMessages{i});
        end
    end
    
    fprintf('Data: %d total ROIs (%.1f MB)\n', summary.data.totalROIs, summary.data.totalSize_MB);
    fprintf('ROIs per file: %d avg (range: %d-%d)\n', ...
        summary.data.avgROIsPerFile, summary.data.roiRange(1), summary.data.roiRange(2));
    
    fprintf('Timing: %.3f s load + %.3f s process = %.3f s total\n', ...
        summary.timing.loadTime_s, summary.timing.processTime_s, summary.timing.totalTime_s);
    
    fprintf('=====================================\n');
end