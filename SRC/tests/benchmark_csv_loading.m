function pipeline = benchmark_csv_loading()
    % PIPELINE_V2 - Clean version using csv_loader_v2
    
    pipeline = struct();
    pipeline.processSingleFile = @processSingleFile;
    pipeline.processBatchFiles = @processBatchFiles;
    pipeline.getLoadStats = @getLoadStats;
    pipeline.setOptions = @setOptions;
    pipeline.options = getDefaultOptions();
    
    function options = getDefaultOptions()
        options = struct(...
            'skipHeaders', true, ...
            'expectedFrames', 1200, ...
            'dataType', 'single', ...
            'useParallel', true, ...
            'maxWorkers', [], ...
            'verbose', true);
    end

    function result = processSingleFile(filepath, varargin)
        if nargin >= 2 && isstruct(varargin{1})
            options = varargin{1};
        else
            options = getDefaultOptions();
        end
        
        if options.verbose
            fprintf('Loading single file: %s\n', filepath);
        end
        
        % Use NEW loader explicitly
        loader = csv_loader_v2();
        
        tic;
        [data, metadata] = loader.loadSingleFile(filepath, options);
        loadTime = toc;
        
        roiHeaders = loader.generateROIHeaders(metadata.numROIs);
        
        dataSize_MB = numel(data) * 4 / (1024^2);
        loadRate_MB_s = dataSize_MB / loadTime;
        
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
        if nargin >= 2 && isstruct(varargin{1})
            options = varargin{1};
        else
            options = getDefaultOptions();
        end
        
        if options.verbose
            fprintf('Starting batch processing in folder: %s\n', folder);
        end
        
        % Use NEW loader explicitly
        loader = csv_loader_v2();
        
        tic;
        [allData, allMetadata] = loader.loadBatchFiles(folder, options);
        totalLoadTime = toc;
        
        numFiles = length(allData);
        allROIHeaders = cell(numFiles, 1);
        
        for i = 1:numFiles
            if ~isempty(allData{i})
                allROIHeaders{i} = loader.generateROIHeaders(allMetadata(i).numROIs);
            end
        end
        
        batchStats = calculateBatchStats(allData, allMetadata, totalLoadTime, options.verbose);
        
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
        validFiles = ~cellfun(@isempty, allData);
        numValid = sum(validFiles);
        
        if numValid == 0
            stats = struct('totalFiles', length(allData), 'validFiles', 0);
            return;
        end
        
        numFrames = arrayfun(@(x) x.numFrames, allMetadata(validFiles));
        numROIs = arrayfun(@(x) x.numROIs, allMetadata(validFiles));
        
        fileSizes_MB = arrayfun(@(data) numel(data{1}) * 4 / (1024^2), allData(validFiles));
        totalSize_MB = sum(fileSizes_MB);
        
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
        persistent savedOptions;
        if isempty(savedOptions)
            savedOptions = getDefaultOptions();
        end
        
        if nargin == 0
            currentOptions = savedOptions;
            return;
        end
        
        if nargin == 1 && isstruct(varargin{1})
            newOptions = varargin{1};
            fields = fieldnames(newOptions);
            for i = 1:length(fields)
                savedOptions.(fields{i}) = newOptions.(fields{i});
            end
        else
            for i = 1:2:length(varargin)
                if i+1 <= length(varargin)
                    savedOptions.(varargin{i}) = varargin{i+1};
                end
            end
        end
        
        currentOptions = savedOptions;
    end

    function loadStats = getLoadStats(results)
        if isfield(results, 'batchStats')
            stats = results.batchStats;
            loadStats = struct(...
                'type', 'batch', ...
                'numFiles', results.numFiles, ...
                'totalSize_MB', stats.totalDataSize_MB, ...
                'loadTime_s', stats.totalLoadTime_s, ...
                'loadRate_MB_s', stats.avgLoadRate_MB_s, ...
                'totalROIs', stats.roiStats.total);
        else
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