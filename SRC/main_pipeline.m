function results = main_pipeline(folder, options)
    % MAIN_PIPELINE - Orchestrates fluorescent imaging analysis with baseline calculation
    % FIXED VERSION: Better error handling for empty results
    
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
    
    % Determine processing mode
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
        
        % Create summary statistics with IMPROVED error handling
        try
            results.summary = createSummaryStats(fileResults, metadataArray, loadTime, processTime);
        catch ME
            warning('Failed to create summary statistics: %s', E.message);
            % Create minimal summary even if creation fails
            results.summary = createMinimalSummary(fileResults, metadataArray, loadTime, processTime);
        end
        
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
        loader = csv_loader_v2();
        
        loadOptions = struct();
        loadOptions.useParallel = options.useParallel && strcmp(processingMode, 'batch');
        loadOptions.dataType = 'single';
        loadOptions.expectedFrames = 1200;
        
        if strcmp(processingMode, 'single_file')
            [singleData, singleMetadata] = loader.loadSingleFile(folder, loadOptions);
            dataCell = {singleData};
            metadataArray = singleMetadata;
        else
            [dataCell, metadataArray] = loader.loadBatchFiles(folder, loadOptions);
        end
        
    catch ME
        warning('Failed to load data files: %s', E.message);
        dataCell = {};
        metadataArray = [];
    end
end

function [fileResults, errors] = processFilesWithBaseline(dataCell, metadataArray, options)
    % Process each file with baseline calculation and dF/F normalization
    
    numFiles = length(dataCell);
    fileResults = cell(numFiles, 1);
    
    errors = struct();
    errors.failedFiles = {};
    errors.errorMessages = {};
    
    if options.verbose
        fprintf('  Processing %d files with baseline calculation...\n', numFiles);
    end
    
    for i = 1:numFiles
        try
            data = dataCell{i};
            metadata = metadataArray(i);
            
            if options.verbose && numFiles > 1 && mod(i, 2) == 0
                fprintf('    Processing file %d/%d (%s)\n', i, numFiles, metadata.filename);
            elseif options.verbose && numFiles == 1
                fprintf('    Processing file: %s\n', metadata.filename);
            end
            
            % Process this file with IMPROVED error handling
            fileResult = processSingleFileWithBaseline(data, metadata, options);
            fileResults{i} = fileResult;
            
        catch ME
            if options.continueOnError
                if options.verbose
                    fprintf('    WARNING: Failed to process file %d (%s): %s\n', ...
                        i, metadataArray(i).filename, ME.message);
                    % Show more detailed error info
                    if ~isempty(ME.stack)
                        fprintf('      Error in: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
                    end
                end
                
                errors.failedFiles{end+1} = metadataArray(i).filename;
                errors.errorMessages{end+1} = ME.message;
                fileResults{i} = [];
            else
                rethrow(ME);
            end
        end
    end
    
    % Remove empty results from failed files
    validResults = ~cellfun(@isempty, fileResults);
    fileResults = fileResults(validResults);
end

function fileResult = processSingleFileWithBaseline(data, metadata, options)
    % UPDATED VERSION: Uses new modular event detection
    % REPLACE this function in your existing main_pipeline.m
    
    [numFrames, numROIs] = size(data);
    
    % Load configuration (now includes Schmitt trigger parameters)
    config = tracenorm_config();
    config.verbose = options.verbose && numROIs > 2000;
    
    if options.verbose && numROIs > 2000
        fprintf('      Large file: %d ROIs, %.1f MB - applying baseline calculation\n', ...
            numROIs, numel(data) * 4 / (1024^2));
    end
    
    %% === Baseline Calculation (unchanged) ===
    try
        tic;
        [baseline, outlier_mask, baseline_stats] = baseline_detector(data, config);
        baseline_time = toc;
    catch ME
        error('Baseline calculation failed: %s', ME.message);
    end
    
    %% === dF/F Calculation (updated - no event detection) ===
    try
        tic;
        [dfof_data, dfof_stats] = dfof_calculator(data, baseline, config);
        dfof_time = toc;
    catch ME
        error('dF/F calculation failed: %s', ME.message);
    end
    
    %% === Event Detection (NEW - Schmitt trigger) ===
    try
        tic;
        [event_mask, event_stats] = schmitt_event_detector(dfof_data, config);
        event_time = toc;
    catch ME
        error('Event detection failed: %s', ME.message);
    end
    
    %% === Quality Assessment (NEW - comprehensive) ===
    try
        tic;
        quality_metrics = quality_assessor(data, baseline_stats, dfof_stats, event_stats, config);
        quality_time = toc;
    catch ME
        error('Quality assessment failed: %s', ME.message);
    end
    
    %% === Visualization (updated to handle events) ===
    plot_handles = [];
    plot_time = 0;
    
    if options.createPlots
        try
            tic;
            % Use existing baseline_plotter (it should still work)
            % or create enhanced version that shows events
            plot_handles = baseline_plotter(data, baseline, dfof_data, outlier_mask, ...
                baseline_stats, dfof_stats, metadata, config, event_stats);
            plot_time = toc;
            
            if options.savePlots
                save_baseline_plots(plot_handles, metadata, options);
            end
        catch ME
            warning('Plotting failed: %s', ME.message);
            if options.verbose
                fprintf('      Skipping plots due to error: %s\n', ME.message);
            end
            plot_handles = [];
            plot_time = 0;
        end
    end
    
    %% === Compile Results (updated structure) ===
    fileResult = struct();
    
    % Basic info (unchanged)
    fileResult.filename = metadata.filename;
    fileResult.numFrames = numFrames;
    fileResult.numROIs = numROIs;
    fileResult.dataSize_MB = numel(data) * 4 / (1024^2);
    
    % Results (updated with new modules)
    fileResult.baseline = baseline;
    fileResult.outlier_mask = outlier_mask;
    fileResult.baseline_stats = baseline_stats;
    fileResult.dfof_data = dfof_data;
    fileResult.dfof_stats = dfof_stats;
    
    % NEW: Event detection results
    fileResult.event_mask = event_mask;
    fileResult.event_stats = event_stats;
    
    % NEW: Comprehensive quality metrics  
    fileResult.quality_metrics = quality_metrics;
    
    % Timing (updated)
    fileResult.timing = struct();
    fileResult.timing.baseline_time = baseline_time;
    fileResult.timing.dfof_time = dfof_time;
    fileResult.timing.event_time = event_time;        % NEW
    fileResult.timing.quality_time = quality_time;    % NEW
    fileResult.timing.plot_time = plot_time;
    fileResult.timing.total_time = baseline_time + dfof_time + event_time + quality_time + plot_time;
    
    % Visualization
    fileResult.plot_handles = plot_handles;
    
    % Summary metrics (for backward compatibility with existing code)
    fileResult.quality = struct();
    fileResult.quality.fraction_good_baseline = getFieldSafe(baseline_stats, 'fraction_good_rois', 0);
    fileResult.quality.fraction_good_dfof = getFieldSafe(quality_metrics, 'fraction_good_rois', 0);
    fileResult.quality.mean_snr = getFieldSafe(dfof_stats, 'mean_snr', 0);
    fileResult.quality.transport_rois = getFieldSafe(baseline_stats, 'num_transport_rois', 0);
    
    % NEW: Enhanced quality metrics
    fileResult.quality.active_rois = getFieldSafe(quality_metrics.summary, 'active_rois', 0);
    fileResult.quality.total_events = getFieldSafe(event_stats, 'total_events', 0);
    fileResult.quality.mean_event_amplitude = getFieldSafe(event_stats, 'mean_event_amplitude', 0);
    
    if options.verbose && (fileResult.quality.transport_rois > 0 || fileResult.quality.active_rois > 0)
        fprintf('      Results: %d transport ROIs, %d active ROIs, %d total events, SNR=%.2f\n', ...
            fileResult.quality.transport_rois, fileResult.quality.active_rois, ...
            fileResult.quality.total_events, fileResult.quality.mean_snr);
    end
end

% Keep this helper function (unchanged)
function value = getFieldSafe(structure, fieldName, defaultValue)
    % Safely get field value with default fallback
    if isfield(structure, fieldName)
        value = structure.(fieldName);
    else
        value = defaultValue;
    end
end

function summary = createMinimalSummary(fileResults, metadataArray, loadTime, processTime)
    % Create minimal summary when full summary creation fails
    
    summary = struct();
    summary.totalFiles = length(metadataArray);
    summary.validFiles = length(fileResults);
    summary.failedFiles = summary.totalFiles - summary.validFiles;
    
    summary.timing = struct();
    summary.timing.loadTime_s = loadTime;
    summary.timing.processTime_s = processTime;
    summary.timing.totalTime_s = loadTime + processTime;
    
    if summary.validFiles > 0
        try
            numROIs = arrayfun(@(r) r.numROIs, [fileResults{:}]);
            dataSizes = arrayfun(@(r) r.dataSize_MB, [fileResults{:}]);
            
            summary.data = struct();
            summary.data.totalROIs = sum(numROIs);
            summary.data.totalSize_MB = sum(dataSizes);
            summary.data.avgROIsPerFile = round(mean(numROIs));
            summary.data.roiRange = [min(numROIs), max(numROIs)];
        catch
            summary.data = struct('totalROIs', 0, 'totalSize_MB', 0, 'avgROIsPerFile', 0, 'roiRange', [0, 0]);
        end
    else
        summary.data = struct('totalROIs', 0, 'totalSize_MB', 0, 'avgROIsPerFile', 0, 'roiRange', [0, 0]);
    end
end

function summary = createSummaryStats(fileResults, metadataArray, loadTime, processTime)
    % Create comprehensive summary with better error handling
    
    if isempty(fileResults)
        summary = createMinimalSummary(fileResults, metadataArray, loadTime, processTime);
        return;
    end
    
    try
        % Extract basic metrics
        numROIs = arrayfun(@(r) r.numROIs, [fileResults{:}]);
        dataSizes = arrayfun(@(r) r.dataSize_MB, [fileResults{:}]);
        
        % Quality metrics with safe extraction
        baseline_quality = arrayfun(@(r) getFieldSafe(r.quality, 'fraction_good_baseline', 0), [fileResults{:}]);
        dfof_quality = arrayfun(@(r) getFieldSafe(r.quality, 'fraction_good_dfof', 0), [fileResults{:}]);
        mean_snrs = arrayfun(@(r) getFieldSafe(r.quality, 'mean_snr', 0), [fileResults{:}]);
        transport_counts = arrayfun(@(r) getFieldSafe(r.quality, 'transport_rois', 0), [fileResults{:}]);
        active_counts = arrayfun(@(r) getFieldSafe(r.quality, 'active_rois', 0), [fileResults{:}]);
        
        % Timing metrics
        baseline_times = arrayfun(@(r) getFieldSafe(r.timing, 'baseline_time', 0), [fileResults{:}]);
        dfof_times = arrayfun(@(r) getFieldSafe(r.timing, 'dfof_time', 0), [fileResults{:}]);
        plot_times = arrayfun(@(r) getFieldSafe(r.timing, 'plot_time', 0), [fileResults{:}]);
        
        % Build summary
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
        
    catch ME
        warning('Error creating detailed summary: %s', E.message);
        summary = createMinimalSummary(fileResults, metadataArray, loadTime, processTime);
    end
end

function save_baseline_plots(plot_handles, metadata, options)
    % Save baseline validation plots to files
    
    if isempty(plot_handles)
        return;
    end
    
    try
        [~, filename_base, ~] = fileparts(metadata.filename);
        output_dir = fullfile('baseline_plots', filename_base);
        if ~exist(output_dir, 'dir')
            mkdir(output_dir);
        end
        
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
    catch ME
        warning('Failed to save plots: %s', E.message);
    end
end

function printSummary(summary, errors)
    % Print formatted summary with better error handling
    
    fprintf('\n=== Pipeline Processing Summary ===\n');
    
    % Handle case where summary might be incomplete
    if isfield(summary, 'totalFiles') && isfield(summary, 'validFiles')
        fprintf('Files: %d/%d processed successfully\n', summary.validFiles, summary.totalFiles);
        
        if isfield(summary, 'failedFiles') && summary.failedFiles > 0
            fprintf('FAILED: %d files\n', summary.failedFiles);
            if isfield(errors, 'failedFiles') && ~isempty(errors.failedFiles)
                for i = 1:length(errors.failedFiles)
                    fprintf('  - %s: %s\n', errors.failedFiles{i}, errors.errorMessages{i});
                end
            end
        end
    end
    
    % Data summary
    if isfield(summary, 'data')
        fprintf('Data: %d total ROIs (%.1f MB)\n', summary.data.totalROIs, summary.data.totalSize_MB);
        fprintf('ROIs per file: %d avg (range: %d-%d)\n', ...
            summary.data.avgROIsPerFile, summary.data.roiRange(1), summary.data.roiRange(2));
    end
    
    % Timing summary
    if isfield(summary, 'timing')
        fprintf('Timing: %.3f s load + %.3f s process = %.3f s total\n', ...
            summary.timing.loadTime_s, summary.timing.processTime_s, summary.timing.totalTime_s);
        
        if isfield(summary.timing, 'avgBaselineTime_s')
            fprintf('  Baseline: %.3f s avg, dF/F: %.3f s avg, Plots: %.3f s avg\n', ...
                summary.timing.avgBaselineTime_s, summary.timing.avgDfofTime_s, summary.timing.avgPlotTime_s);
        end
    end
    
    % Quality summary
    if isfield(summary, 'quality')
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