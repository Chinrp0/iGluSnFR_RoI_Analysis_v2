function results = main_pipeline(folder, options)
    % MAIN_PIPELINE - Orchestrates fluorescent imaging analysis with baseline calculation
    % Implements Option 1: Process per file, then combine
    % NOW INCLUDES: Iterative rolling median baseline calculation and dF/F normalization
    %
    % Inputs:
    %   folder  - Path to folder containing CSV files OR single CSV file path
    %   options - Configuration struct (optional)
    %             .verbose (logical, default true)
    %             .useParallel (logical, default true)
    %             .continueOnError (logical, default true)
    %             .createPlots (logical, default true)
    %             .savePlots (logical, default false)
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
    if ~isfield(options, 'createPlots'), options.createPlots = true; end
    if ~isfield(options, 'savePlots'), options.savePlots = false; end
    
    % Determine if single file or folder processing
    if isfile(folder) && endsWith(folder, '.csv')
        processingMode = 'single_file';
        if options.verbose
            fprintf('=== Single File Analysis Mode ===\n');
            fprintf('Processing file: %s\n', folder);
        end
    else
        processingMode = 'batch';
        if options.verbose
            fprintf('=== Batch Analysis Mode ===\n');
            fprintf('Processing folder: %s\n', folder);
        end
    end
    
    % Initialize results structure
    results = struct();
    results.fileResults = {};
    results.metadata = [];
    results.errors = struct('failedFiles', {}, 'errorMessages', {});
    results.summary = struct();
    results.processingMode = processingMode;
    
    try
        %% Step 1: Load CSV Data
        if options.verbose
            fprintf('\nStep 1: Loading CSV files...\n');
        end
        
        tic;
        [dataCell, metadataArray] = loadDataFiles(folder, options, processingMode);
        loadTime = toc;
        
        if isempty(dataCell)
            error('No valid CSV files found or all files failed to load');
        end
        
        if options.verbose
            fprintf('  Loaded %d files in %.3f s\n', length(dataCell), loadTime);
        end
        
        %% Step 2: Process Each File with Baseline Calculation
        if options.verbose
            fprintf('\nStep 2: Processing files with baseline calculation...\n');
        end
        
        tic;
        [fileResults, processingErrors] = processFilesWithBaseline(dataCell, metadataArray, options);
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

function [dataCell, metadataArray] = loadDataFiles(folder, options, processingMode)
    % Load CSV files using optimized loader
    
    try
        % Use the proven csv_loader_v2
        loader = csv_loader_v2();
        
        % Configure loading options
        loadOptions = struct();
        loadOptions.useParallel = options.useParallel && strcmp(processingMode, 'batch');
        loadOptions.dataType = 'single';  % Memory efficient
        loadOptions.expectedFrames = 1200;
        
        if strcmp(processingMode, 'single_file')
            % Load single file
            [singleData, singleMetadata] = loader.loadSingleFile(folder, loadOptions);
            dataCell = {singleData};
            metadataArray = singleMetadata;
        else
            % Load batch files
            [dataCell, metadataArray] = loader.loadBatchFiles(folder, loadOptions);
        end
        
    catch ME
        warning('Failed to load data files: %s', ME.message);
        dataCell = {};
        metadataArray = [];
    end
end

function [fileResults, errors] = processFilesWithBaseline(dataCell, metadataArray, options)
    % Process each file with baseline calculation and dF/F normalization
    
    numFiles = length(dataCell);
    fileResults = cell(numFiles, 1);
    
    % Track errors
    errors = struct();
    errors.failedFiles = {};
    errors.errorMessages = {};
    
    if options.verbose
        fprintf('  Processing %d files with baseline calculation...\n', numFiles);
    end
    
    for i = 1:numFiles
        try
            % Get file data and metadata
            data = dataCell{i};          % [1200 x N_ROIs] matrix
            metadata = metadataArray(i); % File metadata
            
            if options.verbose && numFiles > 1 && mod(i, 2) == 0  % Progress every 2 files in batch mode
                fprintf('    Processing file %d/%d (%s)\n', i, numFiles, metadata.filename);
            elseif options.verbose && numFiles == 1  % Always show for single file
                fprintf('    Processing file: %s\n', metadata.filename);
            end
            
            % Process this file with baseline calculation
            fileResult = processSingleFileWithBaseline(data, metadata, options);
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

function fileResult = processSingleFileWithBaseline(data, metadata, options)
    % Process a single file's data matrix with baseline calculation
    % Input: data is [1200 x N_ROIs] matrix for this file
    %
    % NEW: Implements iterative rolling median baseline calculation
    %      and dF/F normalization with quality control
    
    [numFrames, numROIs] = size(data);
    
    % Load baseline calculation configuration
    config = tracenorm_config();
    config.verbose = options.verbose && numROIs > 2000;  % Verbose for large files
    
    if options.verbose && numROIs > 2000
        fprintf('      Large file: %d ROIs, %.1f MB - applying baseline calculation\n', ...
            numROIs, numel(data) * 4 / (1024^2));
    end
    
    %% === Baseline Calculation ===
    tic;
    [baseline, outlier_mask, baseline_stats] = baseline_detector(data, config);
    baseline_time = toc;
    
    %% === dF/F Calculation ===
    tic;
    [dfof_data, dfof_stats] = dfof_calculator(data, baseline, config);
    dfof_time = toc;
    
    %% === Create Visualization (if requested) ===
    plot_handles = [];
    plot_time = 0;
    
    if options.createPlots
        tic;
        plot_handles = baseline_plotter(data, baseline, dfof_data, outlier_mask, ...
            baseline_stats, metadata, config);
        plot_time = toc;
        
        % Save plots if requested
        if options.savePlots
            save_baseline_plots(plot_handles, metadata, options);
        end
    end
    
    %% === Compile File Results ===
    fileResult = struct();
    
    % Basic file information
    fileResult.filename = metadata.filename;
    fileResult.numFrames = numFrames;
    fileResult.numROIs = numROIs;
    fileResult.dataSize_MB = numel(data) * 4 / (1024^2);  % Single precision
    
    % Baseline calculation results
    fileResult.baseline = baseline;                % [1200 x N_ROIs]
    fileResult.outlier_mask = outlier_mask;        % [1200 x N_ROIs] logical
    fileResult.baseline_stats = baseline_stats;    % Baseline quality metrics
    
    % dF/F normalization results  
    fileResult.dfof_data = dfof_data;              % [1200 x N_ROIs] normalized
    fileResult.dfof_stats = dfof_stats;            % dF/F quality metrics
    
    % Processing timing
    fileResult.timing = struct();
    fileResult.timing.baseline_time = baseline_time;
    fileResult.timing.dfof_time = dfof_time;
    fileResult.timing.plot_time = plot_time;
    fileResult.timing.total_time = baseline_time + dfof_time + plot_time;
    
    % Visualization
    fileResult.plot_handles = plot_handles;
    
    % Quality summary
    fileResult.quality = struct();
    fileResult.quality.fraction_good_baseline = baseline_stats.fraction_good_rois;
    fileResult.quality.fraction_good_dfof = dfof_stats.fraction_good_rois;
    fileResult.quality.mean_snr = dfof_stats.mean_snr;
    fileResult.quality.transport_rois = baseline_stats.num_transport_rois;
    fileResult.quality.active_rois = dfof_stats.event_summary.rois_with_events;
    
    if options.verbose && (baseline_stats.num_transport_rois > 0 || dfof_stats.event_summary.rois_with_events > 0)
        fprintf('      Results: %d transport ROIs, %d active ROIs, SNR=%.2f\n', ...
            baseline_stats.num_transport_rois, dfof_stats.event_summary.rois_with_events, ...
            dfof_stats.mean_snr);
    end
end

function save_baseline_plots(plot_handles, metadata, options)
    % Save baseline validation plots to files
    
    if isempty(plot_handles)
        return;
    end
    
    % Create output directory
    [~, filename_base, ~] = fileparts(metadata.filename);
    output_dir = fullfile('baseline_plots', filename_base);
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    % Save each plot
    plot_names = fieldnames(plot_handles);
    for i = 1:length(plot_names)
        plot_name = plot_names{i};
        fig_handle = plot_handles.(plot_name);
        
        if ishandle(fig_handle)
            output_file = fullfile(output_dir, sprintf('%s_%s.png', filename_base, plot_name));
            saveas(fig_handle, output_file, 'png');
            
            if options.verbose
                fprintf('        Saved: %s\n', output_file);
            end
        end
    end
end

function summary = createSummaryStats(fileResults, metadataArray, loadTime, processTime)
    % Create overall processing summary including baseline statistics
    
    if isempty(fileResults)
        summary = struct('totalFiles', 0, 'validFiles', 0);
        return;
    end
    
    % Extract statistics from valid results
    numROIs = arrayfun(@(r) r.numROIs, [fileResults{:}]);
    dataSizes = arrayfun(@(r) r.dataSize_MB, [fileResults{:}]);
    
    % Baseline quality metrics
    baseline_quality = arrayfun(@(r) r.quality.fraction_good_baseline, [fileResults{:}]);
    dfof_quality = arrayfun(@(r) r.quality.fraction_good_dfof, [fileResults{:}]);
    mean_snrs = arrayfun(@(r) r.quality.mean_snr, [fileResults{:}]);
    transport_counts = arrayfun(@(r) r.quality.transport_rois, [fileResults{:}]);
    active_counts = arrayfun(@(r) r.quality.active_rois, [fileResults{:}]);
    
    % Processing times
    baseline_times = arrayfun(@(r) r.timing.baseline_time, [fileResults{:}]);
    dfof_times = arrayfun(@(r) r.timing.dfof_time, [fileResults{:}]);
    plot_times = arrayfun(@(r) r.timing.plot_time, [fileResults{:}]);
    
    summary = struct();
    summary.totalFiles = length(metadataArray);
    summary.validFiles = length(fileResults);
    summary.failedFiles = summary.totalFiles - summary.validFiles;
    
    summary.timing = struct();
    summary.timing.loadTime_s = loadTime;
    summary.timing.processTime_s = processTime;
    summary.timing.totalTime_s = loadTime + processTime;
    summary.timing.avgBaselineTime_s = mean(baseline_times);
    summary.timing.avgDfofTime_s = mean(dfof_times);
    summary.timing.avgPlotTime_s = mean(plot_times);
    
    summary.data = struct();
    summary.data.totalROIs = sum(numROIs);
    summary.data.totalSize_MB = sum(dataSizes);
    summary.data.avgROIsPerFile = round(mean(numROIs));
    summary.data.roiRange = [min(numROIs), max(numROIs)];
    
    summary.quality = struct();
    summary.quality.avgBaselineQuality = mean(baseline_quality);
    summary.quality.avgDfofQuality = mean(dfof_quality);
    summary.quality.avgSNR = mean(mean_snrs);
    summary.quality.totalTransportROIs = sum(transport_counts);
    summary.quality.totalActiveROIs = sum(active_counts);
    summary.quality.fractionActiveROIs = sum(active_counts) / sum(numROIs);
end

function printSummary(summary, errors)
    % Print formatted processing summary including baseline results
    
    fprintf('\n=== Pipeline Processing Summary ===\n');
    fprintf('Files: %d/%d processed successfully\n', summary.validFiles, summary.totalFiles);
    
    if isfield(summary, 'failedFiles') && summary.failedFiles > 0
        fprintf('FAILED: %d files\n', summary.failedFiles);
        if isfield(errors, 'failedFiles') && ~isempty(errors.failedFiles)
            for i = 1:length(errors.failedFiles)
                fprintf('  - %s: %s\n', errors.failedFiles{i}, errors.errorMessages{i});
            end
        end
    end
    
    fprintf('Data: %d total ROIs (%.1f MB)\n', summary.data.totalROIs, summary.data.totalSize_MB);
    fprintf('ROIs per file: %d avg (range: %d-%d)\n', ...
        summary.data.avgROIsPerFile, summary.data.roiRange(1), summary.data.roiRange(2));
    
    fprintf('Timing: %.3f s load + %.3f s process = %.3f s total\n', ...
        summary.timing.loadTime_s, summary.timing.processTime_s, summary.timing.totalTime_s);
    
    % Only show detailed timing if we have valid results
    if summary.validFiles > 0 && isfield(summary.timing, 'avgBaselineTime_s')
        fprintf('  Baseline: %.3f s avg, dF/F: %.3f s avg, Plots: %.3f s avg\n', ...
            summary.timing.avgBaselineTime_s, summary.timing.avgDfofTime_s, summary.timing.avgPlotTime_s);
    end
    
    if isfield(summary, 'quality') && summary.validFiles > 0
        fprintf('\n=== Baseline & dF/F Quality ===\n');
        fprintf('Baseline Quality: %.1f%% good ROIs\n', 100 * summary.quality.avgBaselineQuality);
        fprintf('dF/F Quality: %.1f%% good ROIs\n', 100 * summary.quality.avgDfofQuality);
        fprintf('Average SNR: %.2f\n', summary.quality.avgSNR);
        fprintf('Transport ROIs: %d total\n', summary.quality.totalTransportROIs);
        fprintf('Active ROIs: %d total (%.1f%% of all ROIs)\n', ...
            summary.quality.totalActiveROIs, 100 * summary.quality.fractionActiveROIs);
    end
    
    fprintf('=====================================\n');
end