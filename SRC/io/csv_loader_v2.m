function loader = csv_loader_v2()
    % CSV_LOADER_V2 - Loader for the ROI mean-intensity tables produced by the
    % Fiji macro (5_GluSnFR_icyROI_to_mean_batch). That macro writes, via
    % "Multi Measure" + Table.save, a table shaped:
    %
    %     <blank>, Mean(roi-001), Mean(roi-002), ...      <- header row
    %     1,       <val>,         <val>,         ...       <- frame 1
    %     2,       ...                                     <- frame 2
    %
    % i.e. ONE header row plus a leading FRAME-INDEX column (1..N). The loader is
    % layout-aware rather than slicing fixed offsets: it skips leading header
    % rows and strips the leading index column ONLY when column 1 is a contiguous
    % 1..N integer index. That makes it correct for the macro CSVs (spont / 1AP /
    % 2AP) AND for the older 1AP .xlsx exports, which instead have TWO header
    % rows (a .czi name row + an ROI-name row) and NO index column - so nothing
    % is dropped from those.

    loader = struct();
    loader.loadSingleFile   = @loadSingle;
    loader.loadBatchFiles   = @loadBatch;
    loader.generateROIHeaders = @generateHeaders;
end

function [data, metadata] = loadSingle(filepath, options)
    if nargin < 2, options = struct(); end
    if ~isfield(options, 'dataType'), options.dataType = 'single'; end

    if ~exist(filepath, 'file')
        error('csv_loader_v2:missingFile', 'File does not exist: %s', filepath);
    end

    try
        % readmatrix returns only the numeric block (text header rows become NaN
        % or are skipped, depending on format). We then normalise defensively.
        raw = readmatrix(filepath);

        % Drop any entirely-NaN leading/trailing rows (stray header rows) and
        % all-NaN columns (a trailing delimiter). A real fluorescence frame is
        % never NaN across every ROI, so this only removes non-data.
        raw(all(isnan(raw), 2), :) = [];
        raw(:, all(isnan(raw), 1)) = [];

        if isempty(raw)
            error('csv_loader_v2:noData', 'No numeric data found in %s', filepath);
        end

        % Strip the ImageJ frame-index column if present (column 1 == 1..N).
        if size(raw, 2) >= 2 && looks_like_index_column(raw(:, 1))
            raw = raw(:, 2:end);
        end

        data = cast(raw, options.dataType);

        [~, filename] = fileparts(filepath);
        metadata = struct( ...
            'filename',  filename, ...
            'filepath',  filepath, ...
            'numFrames', size(data, 1), ...
            'numROIs',   size(data, 2), ...
            'dataType',  options.dataType, ...
            'loadTime',  datetime('now'));
    catch ME
        error('csv_loader_v2:loadFailed', 'Failed to load %s: %s', filepath, ME.message);
    end
end

function tf = looks_like_index_column(col)
    % True when col is exactly the sequence 1,2,...,N - the ImageJ frame index.
    n  = numel(col);
    tf = n >= 2 && all(isfinite(col)) && isequal(col(:), (1:n)');
end

function [allData, allMetadata] = loadBatch(folder, options)
    % Load every data file in a folder (serial; the pipeline is fast enough that
    % a parallel pool is not worth its startup cost).
    if nargin < 2, options = struct(); end
    if ~isfield(options, 'dataType'), options.dataType = 'single'; end

    files    = list_data_files(folder);   % .csv / .xlsx / .xls
    numFiles = numel(files);
    if numFiles == 0
        error('csv_loader_v2:noFiles', 'No CSV/XLSX files found in folder: %s', folder);
    end

    fprintf('Loading %d data files...\n', numFiles);
    allData     = cell(numFiles, 1);
    allMetadata = cell(numFiles, 1);

    for i = 1:numFiles
        filepath = fullfile(files(i).folder, files(i).name);
        try
            [allData{i}, allMetadata{i}] = loadSingle(filepath, options);
        catch ME
            warning('csv_loader_v2:fileFailed', 'Failed to load %s: %s', filepath, ME.message);
            allData{i}     = [];
            allMetadata{i} = struct('filename', files(i).name, 'error', ME.message);
        end
    end

    validIdx  = ~cellfun(@isempty, allData);
    allData   = allData(validIdx);
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
