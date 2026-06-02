function loader = csv_loader_v2()
    % CSV_LOADER_V2 - Fixed version without caching issues
    
    loader = struct();
    loader.loadSingleFile = @loadSingle;
    loader.loadBatchFiles = @loadBatch;
    loader.generateROIHeaders = @generateHeaders;
end

function [data, metadata] = loadSingle(filepath, options)
    if nargin < 2, options = struct(); end
    if ~isfield(options, 'skipHeaders'), options.skipHeaders = true; end
    if ~isfield(options, 'expectedFrames'), options.expectedFrames = 1200; end
    if ~isfield(options, 'dataType'), options.dataType = 'single'; end
    
    if ~exist(filepath, 'file')
        error('File does not exist: %s', filepath);
    end
    
    try
        fullData = readmatrix(filepath);
        data = fullData(2:end, 2:end);
        data = cast(data, options.dataType);
        
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

function [allData, allMetadata] = loadBatch(folder, options)
    if nargin < 2, options = struct(); end
    if ~isfield(options, 'useParallel'), options.useParallel = true; end
    if ~isfield(options, 'maxWorkers'), options.maxWorkers = []; end
    
    % Get data files (.csv / .xlsx / .xls)
    csvFiles = list_data_files(folder);
    numFiles = length(csvFiles);

    if numFiles == 0
        error('No CSV/XLSX files found in folder: %s', folder);
    end
    
    fprintf('Loading %d CSV files...\n', numFiles);
    
    % Use cell arrays throughout
    allData = cell(numFiles, 1);
    allMetadata = cell(numFiles, 1);
    
    if options.useParallel && numFiles > 1
        setupPool(options.maxWorkers);
        parfor i = 1:numFiles
            filepath = fullfile(csvFiles(i).folder, csvFiles(i).name);
            try
                [allData{i}, allMetadata{i}] = loadSingle(filepath, options);
            catch ME
                warning('Failed to load file %s: %s', filepath, ME.message);
                allData{i} = [];
                allMetadata{i} = struct('filename', csvFiles(i).name, 'error', ME.message);
            end
        end
    else
        for i = 1:numFiles
            filepath = fullfile(csvFiles(i).folder, csvFiles(i).name);
            try
                [allData{i}, allMetadata{i}] = loadSingle(filepath, options);
            catch ME
                warning('Failed to load file %s: %s', filepath, ME.message);
                allData{i} = [];
                allMetadata{i} = struct('filename', csvFiles(i).name, 'error', ME.message);
            end
        end
    end
    
    % Filter and convert
    validIdx = ~cellfun(@isempty, allData);
    allData = allData(validIdx);
    validMeta = allMetadata(validIdx);
    
    if ~isempty(validMeta)
        allMetadata = [validMeta{:}];
    else
        allMetadata = struct.empty();
    end
    
    fprintf('Successfully loaded %d/%d files\n', sum(validIdx), numFiles);
end

function headers = generateHeaders(numROIs, startIndex)
    if nargin < 2, startIndex = 1; end
    headers = cell(numROIs, 1);
    for i = 1:numROIs
        headers{i} = sprintf('ROI_%03d', startIndex + i - 1);
    end
end

function setupPool(maxWorkers)
    pool = gcp('nocreate');
    if isempty(pool)
        if isempty(maxWorkers)
            numCores = feature('numcores');
            maxWorkers = min(numCores - 1, 8);
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