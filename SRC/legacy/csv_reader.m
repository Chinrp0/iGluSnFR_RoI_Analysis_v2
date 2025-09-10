function reader = csv_reader()
    % CSV_READER - Optimized CSV file reading for fluorescent imaging data
    % Supports both single-file and batch parallel processing
    
    reader = struct();
    reader.loadSingleFile = @csv_loadSingleFile;
    reader.loadBatchFiles = @csv_loadBatchFiles;
    reader.getCSVFiles = @csv_getCSVFiles;
    reader.generateROIHeaders = @csv_generateROIHeaders;
end

function [data, metadata] = csv_loadSingleFile(filepath, options)
    % Load single CSV file with optimized I/O
    % 
    % Inputs:
    %   filepath - Path to CSV file
    %   options  - struct with fields (optional):
    %              .skipHeaders (logical, default true)
    %              .expectedFrames (numeric, default 1200)
    %              .dataType (char, default 'single')
    %
    % Outputs:
    %   data     - Numeric matrix [frames x ROIs]
    %   metadata - struct with file info and dimensions
    
    % Handle optional options argument
    if nargin < 2
        options = struct();
    end
    
    % Set defaults
    if ~isfield(options, 'skipHeaders'), options.skipHeaders = true; end
    if ~isfield(options, 'expectedFrames'), options.expectedFrames = 1200; end
    if ~isfield(options, 'dataType'), options.dataType = 'single'; end
    
    if ~exist(filepath, 'file')
        error('File does not exist: %s', filepath);
    end
    
    try
        if options.skipHeaders
            % Direct numeric loading - skip row 1, use column 2:end
            fullData = readmatrix(filepath);
            data = fullData(2:end, 2:end);  % Skip header row and frame column
            
            % Verify dimensions
            if size(data, 1) ~= options.expectedFrames
                warning('Expected %d frames, found %d in %s', ...
                    options.expectedFrames, size(data, 1), filepath);
            end
        else
            % Load with headers for validation (same result, but explicit)
            fullData = readmatrix(filepath);
            data = fullData(2:end, 2:end);  % Skip frame column and header row
        end
        
        % Convert to specified data type for memory efficiency
        data = cast(data, options.dataType);
        
        % Create metadata
        [~, filename, ~] = fileparts(filepath);
        metadata = struct(...
            'filename', filename, ...
            'filepath', filepath, ...
            'numFrames', size(data, 1), ...
            'numROIs', size(data, 2), ...
            'dataType', options.dataType, ...
            'loadTime', datetime('now'));
            
    catch ME
        error('Failed to load CSV file %s: %s', filepath, ME.message);
    end
end

function [allData, allMetadata] = csv_loadBatchFiles(folder, options)
    % Load multiple CSV files in parallel
    %
    % Inputs:
    %   folder  - Directory containing CSV files
    %   options - Same as loadSingleFile, plus:
    %             .useParallel (logical, default true)
    %             .maxWorkers (numeric, default auto)
    %
    % Outputs:
    %   allData     - Cell array of data matrices
    %   allMetadata - Array of metadata structs
    
    % Handle optional options argument
    if nargin < 2
        options = struct();
    end
    
    % Set defaults
    if ~isfield(options, 'skipHeaders'), options.skipHeaders = true; end
    if ~isfield(options, 'expectedFrames'), options.expectedFrames = 1200; end
    if ~isfield(options, 'dataType'), options.dataType = 'single'; end
    if ~isfield(options, 'useParallel'), options.useParallel = true; end
    if ~isfield(options, 'maxWorkers'), options.maxWorkers = []; end
    
    % Get all CSV files
    csvFiles = csv_getCSVFiles(folder);
    numFiles = length(csvFiles);
    
    if numFiles == 0
        error('No CSV files found in folder: %s', folder);
    end
    
    fprintf('Loading %d CSV files...\n', numFiles);
    
    % Initialize outputs
    allData = cell(numFiles, 1);
    allMetadata(numFiles) = struct();
    
    if options.useParallel && numFiles > 1
        % Setup parallel pool if needed
        csv_setupParallelPool(options.maxWorkers);
        
        % Parallel loading - one file per worker
        parfor i = 1:numFiles
            filepath = fullfile(csvFiles(i).folder, csvFiles(i).name);
            try
                [allData{i}, allMetadata(i)] = csv_loadSingleFile(filepath, options);
            catch ME
                warning('Failed to load file %s: %s', filepath, ME.message);
                allData{i} = [];
                allMetadata(i) = struct('filename', csvFiles(i).name, 'error', ME.message);
            end
        end
    else
        % Sequential loading
        for i = 1:numFiles
            filepath = fullfile(csvFiles(i).folder, csvFiles(i).name);
            try
                [allData{i}, allMetadata(i)] = csv_loadSingleFile(filepath, options);
            catch ME
                warning('Failed to load file %s: %s', filepath, ME.message);
                allData{i} = [];
                allMetadata(i) = struct('filename', csvFiles(i).name, 'error', ME.message);
            end
        end
    end
    
    % Remove failed loads (check both empty data and error field)
    validIdx = ~cellfun(@isempty, allData) & cellfun(@isempty, {allMetadata.error});
    allData = allData(validIdx);
    allMetadata = allMetadata(validIdx);
    
    fprintf('Successfully loaded %d/%d files\n', sum(validIdx), numFiles);
end

function csvFiles = csv_getCSVFiles(folder)
    % Get all CSV files in folder
    
    if ~exist(folder, 'dir')
        error('Folder does not exist: %s', folder);
    end
    
    csvFiles = dir(fullfile(folder, '*.csv'));
    
    if isempty(csvFiles)
        warning('No CSV files found in folder: %s', folder);
    end
end

function headers = csv_generateROIHeaders(numROIs, startIndex)
    % Generate standardized ROI headers programmatically
    %
    % Inputs:
    %   numROIs    - Number of ROIs
    %   startIndex - Starting ROI number (default 1)
    %
    % Output:
    %   headers    - Cell array of ROI header strings
    
    if nargin < 2
        startIndex = 1;
    end
    
    headers = cell(numROIs, 1);
    for i = 1:numROIs
        headers{i} = sprintf('ROI_%03d', startIndex + i - 1);
    end
end

function csv_setupParallelPool(maxWorkers)
    % Setup parallel pool with specified or optimal number of workers
    
    % Check if parallel pool exists
    pool = gcp('nocreate');
    
    if isempty(pool)
        if isempty(maxWorkers)
            % Use optimal number based on available cores
            numCores = feature('numcores');
            maxWorkers = min(numCores - 1, 8); % Reserve 1 core, cap at 8
        end
        
        try
            parpool('local', maxWorkers);
            fprintf('Started parallel pool with %d workers\n', maxWorkers);
        catch ME
            warning('Failed to start parallel pool: %s', ME.message);
        end
    else
        fprintf('Using existing parallel pool with %d workers\n', pool.NumWorkers);
    end
end