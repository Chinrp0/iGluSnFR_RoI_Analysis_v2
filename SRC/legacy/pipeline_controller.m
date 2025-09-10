function pipeline = pipeline_controller()
    % PIPELINE_CONTROLLER - Main controller for fluorescent imaging analysis
    % Manages single-file and batch processing workflows
    
    % Create the pipeline struct with function handles to nested functions
    pipeline = struct();
    pipeline.processSingleFile = @processSingleFile;
    pipeline.processBatchFiles = @processBatchFiles;
    pipeline.getLoadStats = @getLoadStats;
    pipeline.setOptions = @setOptions;
    pipeline.options = getDefaultOptions();
    
    function options = getDefaultOptions()
        % Default processing options
        options = struct(...
            'skipHeaders', true, ...
            'expectedFrames', 1200, ...
            'dataType', 'single', ...
            'useParallel', true, ...
            'maxWorkers', [], ...
            'verbose', true);
    end

    function result = processSingleFile(filepath, varargin)
        % Process single CSV file and return results with timing
        %
        % Inputs:
        %   filepath - Path to CSV file
        %   options  - Processing options (optional)
        %
        % Output:
        %   result   - struct with data, metadata, and performance info
        
        if nargin >= 2 && isstruct(varargin{1})
            options = varargin{1};
        else
            options = getDefaultOptions();
        end
        
        if options.verbose
            fprintf('Loading single file: %s\n', filepath);
        end
        
        % Initialize CSV reader
        reader = csv_reader();
        
        % Time the loading operation
        tic;
        [data, metadata] = reader.loadSingleFile(filepath, options);
        loadTime = toc;
        
        % Generate ROI headers
        roiHeaders = reader.generateROIHeaders(metadata.numROIs);
        
        % Performance statistics
        dataSize_MB = numel(data) * 4 / (1024^2); % Assuming single precision
        loadRate_MB_s = dataSize_MB / loadTime;
        
        % Package results
        result = struct(...
            'data', data, ...
            'metadata', metadata, ...
            'roiHeaders', {roiHeaders}, ...
            'performance', struct(...
                'loadTime_s', loadTime, ...
                'dataSize_MB', dataSize_MB, ...
                'loadRate_MB_s', loadRate_MB_s, ...
                'numROIs', metadata.numROIs, ...
                'numFrames', metadata.numFrames));
        
        if options.verbose
            fprintf('  Loaded: %d frames × %d ROIs (%.1f MB) in %.3f s (%.1f MB/s)\n', ...
                metadata.numFrames, metadata.numROIs, dataSize_MB, loadTime, loadRate_MB_s);
        end
    end

    function results = processBatchFiles(folder, varargin)
        % Process multiple CSV files and return batch results
        %
        % Inputs:
        %   folder  - Directory containing CSV files
        %   options - Processing options (optional)
        %
        % Output:
        %   results - struct with batch data, metadata, and performance info
        
        if nargin >= 2 && isstruct(varargin{1})
            options = varargin{1};
        else
            options = getDefaultOptions();
        end
        
        if options.verbose
            fprintf('Starting batch processing in folder: %s\n', folder);
        end
        
        % Initialize CSV reader
        reader = csv_reader();
        
        % Time the batch loading operation
        tic;
        [allData, allMetadata] = reader.loadBatchFiles(folder, options);
        totalLoadTime = toc;
        
        % Generate ROI headers for each file
        numFiles = length(allData);
        allROIHeaders = cell(numFiles, 1);
        
        for i = 1:numFiles
            if ~isempty(allData{i})
                allROIHeaders{i} = reader.generateROIHeaders(allMetadata(i).numROIs);
            end
        end
        
        % Calculate batch performance statistics
        batchStats = calculateBatchStats(allData, allMetadata, totalLoadTime, options.verbose);
        
        % Package results
        results = struct(...
            'data', {allData}, ...
            'metadata', allMetadata, ...
            'roiHeaders', {allROIHeaders}, ...
            'batchStats', batchStats, ...
            'numFiles', numFiles);
        
        if options.verbose
            printBatchSummary(batchStats, numFiles);
        end
    end

    function stats = calculateBatchStats(allData, allMetadata, totalLoadTime, verbose)
        % Calculate comprehensive batch loading statistics
        
        validFiles = ~cellfun(@isempty, allData);
        numValid = sum(validFiles);
        
        if numValid == 0
            stats = struct('totalFiles', length(allData), 'validFiles', 0);
            return;
        end
        
        % Extract dimensions from valid files
        numFrames = arrayfun(@(x) x.numFrames, allMetadata(validFiles));
        numROIs = arrayfun(@(x) x.numROIs, allMetadata(validFiles));
        
        % Calculate data sizes (assuming single precision)
        fileSizes_MB = arrayfun(@(data) numel(data{1}) * 4 / (1024^2), allData(validFiles));
        totalSize_MB = sum(fileSizes_MB);
        
        % Performance metrics
        avgLoadRate = totalSize_MB / totalLoadTime;
        
        stats = struct(...
            'totalFiles', length(allData), ...
            'validFiles', numValid, ...
            'failedFiles', length(allData) - numValid, ...
            'totalLoadTime_s', totalLoadTime, ...
            'totalDataSize_MB', totalSize_MB, ...
            'avgLoadRate_MB_s', avgLoadRate, ...
            'frameStats', struct(...
                'min', min(numFrames), ...
                'max', max(numFrames), ...
                'mean', mean(numFrames)), ...
            'roiStats', struct(...
                'min', min(numROIs), ...
                'max', max(numROIs), ...
                'mean', round(mean(numROIs)), ...
                'total', sum(numROIs)));
    end

    function printBatchSummary(stats, numFiles)
        % Print formatted batch processing summary
        
        fprintf('\n=== Batch Processing Summary ===\n');
        fprintf('Files: %d/%d loaded successfully\n', stats.validFiles, stats.totalFiles);
        
        if stats.failedFiles > 0
            fprintf('WARNING: %d files failed to load\n', stats.failedFiles);
        end
        
        fprintf('Data: %.1f MB total\n', stats.totalDataSize_MB);
        fprintf('Time: %.3f s (%.1f MB/s average)\n', ...
            stats.totalLoadTime_s, stats.avgLoadRate_MB_s);
        
        fprintf('ROIs: %d total (min: %d, max: %d, avg: %d per file)\n', ...
            stats.roiStats.total, stats.roiStats.min, stats.roiStats.max, stats.roiStats.mean);
        
        fprintf('Frames: min: %d, max: %d, avg: %.1f\n', ...
            stats.frameStats.min, stats.frameStats.max, stats.frameStats.mean);
        fprintf('===============================\n\n');
    end

    function currentOptions = setOptions(varargin)
        % Set or modify processing options
        %
        % Usage:
        %   options = setOptions()                    % Get current defaults
        %   options = setOptions('verbose', false)    % Set specific option
        %   options = setOptions(existingOptions)     % Merge with existing
        
        persistent savedOptions;
        if isempty(savedOptions)
            savedOptions = getDefaultOptions();
        end
        
        if nargin == 0
            currentOptions = savedOptions;
            return;
        end
        
        if nargin == 1 && isstruct(varargin{1})
            % Merge with existing options
            newOptions = varargin{1};
            fields = fieldnames(newOptions);
            for i = 1:length(fields)
                savedOptions.(fields{i}) = newOptions.(fields{i});
            end
        else
            % Parse name-value pairs
            for i = 1:2:length(varargin)
                if i+1 <= length(varargin)
                    savedOptions.(varargin{i}) = varargin{i+1};
                end
            end
        end
        
        currentOptions = savedOptions;
    end

    function loadStats = getLoadStats(results)
        % Extract key loading statistics from results
        %
        % Works with both single file and batch results
        
        if isfield(results, 'batchStats')
            % Batch results
            stats = results.batchStats;
            loadStats = struct(...
                'type', 'batch', ...
                'numFiles', results.numFiles, ...
                'totalSize_MB', stats.totalDataSize_MB, ...
                'loadTime_s', stats.totalLoadTime_s, ...
                'loadRate_MB_s', stats.avgLoadRate_MB_s, ...
                'totalROIs', stats.roiStats.total);
        else
            % Single file results
            perf = results.performance;
            loadStats = struct(...
                'type', 'single', ...
                'numFiles', 1, ...
                'totalSize_MB', perf.dataSize_MB, ...
                'loadTime_s', perf.loadTime_s, ...
                'loadRate_MB_s', perf.loadRate_MB_s, ...
                'totalROIs', perf.numROIs);
        end
    end
end