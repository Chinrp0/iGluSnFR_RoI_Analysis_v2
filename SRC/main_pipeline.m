function results = main_pipeline(folder, options)
    % MAIN_PIPELINE - Orchestrates fluorescent imaging analysis with corrected Schmitt trigger
    % UPDATED VERSION: Integrates corrected_schmitt_detector with proper noise estimation
    
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
    if ~isfield(options, 'useCorrectedSchmitt'), options.useCorrectedSchmitt = true; end  % NEW: Use corrected version by default
    
    % Determine processing mode
    if isfile(folder) && endsWith(folder, {'.csv', '.xlsx', '.xls'}, 'IgnoreCase', true)
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
    
    % Display detector being used
    if options.verbose
        if options.useCorrectedSchmitt
            fprintf('Using: CORRECTED Schmitt trigger (includes outliers in noise)\n');
        else
            fprintf('Using: LEGACY Schmitt trigger (excludes outliers from noise)\n');
        end
    end
    
    % Initialize results structure
    results = struct();
    results.fileResults = {};
    results.metadata = [];
    results.errors = struct('failedFiles', {}, 'errorMessages', {});
    results.summary = struct();
    results.processingMode = processingMode;
    if options.useCorrectedSchmitt
        results.detectorUsed = 'corrected_schmitt';
    else
        results.detectorUsed = 'legacy_schmitt';
    end
    
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
            warning('Failed to create summary statistics: %s', ME.message);
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
        loadOptions.dataType = 'single';
        
        if strcmp(processingMode, 'single_file')
            [singleData, singleMetadata] = loader.loadSingleFile(folder, loadOptions);
            dataCell = {singleData};
            metadataArray = singleMetadata;
        else
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
            
            % Process this file with CORRECTED detector
            fileResult = processSingleFileWithCorrectedDetector(data, metadata, options);
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

function fileResult = processSingleFileWithCorrectedDetector(data, metadata, options)
    % UPDATED VERSION: Uses corrected Schmitt detector with proper noise estimation
    % Maintains backward compatibility while using improved algorithm
    
    [numFrames, numROIs] = size(data);
    
    % Load configuration
    config = tracenorm_config();
    config.verbose = options.verbose && numROIs > 2000;
    
    if options.verbose && numROIs > 2000
        fprintf('      Large file: %d ROIs, %.1f MB - applying corrected Schmitt detector\n', ...
            numROIs, numel(data) * 4 / (1024^2));
    end
    
    %% === Baseline Calculation ===
    try
        tic;
        [baseline, outlier_mask, baseline_stats] = baseline_detector(data, config);
        baseline_time = toc;
        
        if config.verbose
            fprintf('      Baseline: %.3f s, %.1f%% outliers detected\n', ...
                baseline_time, 100 * baseline_stats.mean_outlier_fraction);
        end
    catch ME
        error('Baseline calculation failed: %s', ME.message);
    end
    
    %% === FIX 3: Ensure baseline_stats has outlier_mask BEFORE event detection ===
    if ~isfield(baseline_stats, 'outlier_mask')
        if config.verbose
            fprintf('      DEBUG: Adding missing outlier_mask to baseline_stats\n');
        end
        baseline_stats.outlier_mask = outlier_mask;  % Add the outlier_mask that was returned
    end
    
    %% === FIX 2: Debug output (optional - you can remove this once working) ===
    if config.verbose && numROIs > 1000
        fprintf('      baseline_stats fields: %s\n', strjoin(fieldnames(baseline_stats), ', '));
    end
    
    %% === dF/F Calculation ===
    try
        tic;
        [dfof_data, dfof_stats] = dfof_calculator(data, baseline, config);
        dfof_time = toc;
        
        if config.verbose
            fprintf('      dF/F: %.3f s, mean SNR = %.2f\n', ...
                dfof_time, dfof_stats.mean_snr);
        end
    catch ME
        error('dF/F calculation failed: %s', ME.message);
    end
    
    %% === Event Detection - CORRECTED or LEGACY ===
    try
        tic;
        
        if options.useCorrectedSchmitt
            % Use NEW corrected detector (includes outliers in noise calculation)
            [event_mask, event_stats] = corrected_schmitt_detector(dfof_data, config, baseline_stats);
            detector_used = 'corrected_schmitt';
        else
            % Use LEGACY detector (excludes outliers from noise calculation)
            [event_mask, event_stats] = pure_schmitt_trigger_detector(dfof_data, config, baseline_stats);
            detector_used = 'legacy_schmitt';
        end
        
        event_time = toc;
        
        % CRITICAL: Ensure event_mask is properly included
        if ~isfield(event_stats, 'event_mask')
            event_stats.event_mask = event_mask;
        end
        
        if config.verbose
            fprintf('      Events (%s): %.3f s, %d events across %d ROIs\n', ...
                detector_used, event_time, event_stats.total_events, event_stats.rois_with_events);
        end
    catch ME
        error('Event detection failed: %s', ME.message);
    end
    
    %% === Quality Assessment ===
    try
        tic;
        quality_metrics = quality_assessor(data, baseline_stats, dfof_stats, event_stats, config);
        quality_time = toc;
        
        if config.verbose
            fprintf('      Quality: %.3f s, %d active ROIs (%.1f%%)\n', ...
                quality_time, quality_metrics.summary.active_rois, ...
                100 * quality_metrics.summary.fraction_active);
        end
    catch ME
        warning('Quality assessment failed: %s', ME.message);
        % Create minimal quality metrics for backward compatibility
        quality_metrics = create_minimal_quality_metrics(baseline_stats, dfof_stats, event_stats);
        quality_time = 0;
    end
    
    %% === Visualization ===
    plot_handles = [];
    plot_time = 0;
    
    if options.createPlots
        try
            tic;
            % FIXED: Pass all required arguments with consistent data structure
            plot_handles = baseline_plotter_clean(data, baseline, dfof_data, outlier_mask, ...
                baseline_stats, dfof_stats, metadata, config, event_stats);
            plot_time = toc;
            
            if options.savePlots
                save_baseline_plots(plot_handles, metadata, options);
            end
            
            if config.verbose
                fprintf('      Plots: %.3f s, %d figures created\n', ...
                    plot_time, length(fieldnames(plot_handles)));
            end
        catch ME
            warning('Plotting failed: %s', ME.message);
            if config.verbose
                fprintf('      Skipping plots due to error: %s\n', ME.message);
                % Show more detailed error info
                if ~isempty(ME.stack)
                    fprintf('        Error in: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
                end
            end
            plot_handles = [];
            plot_time = 0;
        end
    end
    
    %% === Compile Results with CONSISTENT structure ===
    fileResult = struct();
    
    % Basic info
    fileResult.filename = metadata.filename;
    fileResult.numFrames = numFrames;
    fileResult.numROIs = numROIs;
    fileResult.dataSize_MB = numel(data) * 4 / (1024^2);
    fileResult.detectorUsed = detector_used;  % NEW: Track which detector was used
    
    % Core results
    fileResult.baseline = baseline;
    fileResult.outlier_mask = outlier_mask;
    fileResult.baseline_stats = baseline_stats;
    fileResult.dfof_data = dfof_data;
    fileResult.dfof_stats = dfof_stats;
    
    % Event detection results
    fileResult.event_mask = event_mask;
    fileResult.event_stats = event_stats;
    
    % Quality metrics
    fileResult.quality_metrics = quality_metrics;
    
    % Timing information
    fileResult.timing = struct();
    fileResult.timing.baseline_time = baseline_time;
    fileResult.timing.dfof_time = dfof_time;
    fileResult.timing.event_time = event_time;
    fileResult.timing.quality_time = quality_time;
    fileResult.timing.plot_time = plot_time;
    fileResult.timing.total_time = baseline_time + dfof_time + event_time + quality_time + plot_time;
    
    % Visualization
    fileResult.plot_handles = plot_handles;
    
    % Summary metrics (for backward compatibility)
    fileResult.quality = struct();
    fileResult.quality.fraction_good_baseline = getFieldSafe(baseline_stats, 'fraction_good_rois', 0);
    fileResult.quality.fraction_good_dfof = getFieldSafe(quality_metrics, 'fraction_good_rois', 0);
    fileResult.quality.mean_snr = getFieldSafe(dfof_stats, 'mean_snr', 0);
    fileResult.quality.transport_rois = getFieldSafe(baseline_stats, 'num_transport_rois', 0);
    
    % Enhanced quality metrics with consistent naming
    fileResult.quality.active_rois = getFieldSafe(quality_metrics.summary, 'active_rois', 0);
    fileResult.quality.total_events = getFieldSafe(event_stats, 'total_events', 0);
    fileResult.quality.mean_event_amplitude = getFieldSafe(event_stats, 'mean_event_amplitude', 0);
    
    if options.verbose && (fileResult.quality.transport_rois > 0 || fileResult.quality.active_rois > 0)
        fprintf('      Final Results: %d transport, %d active ROIs, %d events, SNR=%.2f (%s)\n', ...
            fileResult.quality.transport_rois, fileResult.quality.active_rois, ...
            fileResult.quality.total_events, fileResult.quality.mean_snr, detector_used);
    end
end

function quality_metrics = create_minimal_quality_metrics(baseline_stats, dfof_stats, event_stats)
    % Create minimal quality metrics when full assessment fails
    
    quality_metrics = struct();
    
    % Basic summary
    quality_metrics.summary = struct();
    quality_metrics.summary.active_rois = getFieldSafe(event_stats, 'rois_with_events', 0);
    quality_metrics.summary.total_rois = size(dfof_stats.dfof_mean, 2);
    quality_metrics.summary.fraction_active = quality_metrics.summary.active_rois / quality_metrics.summary.total_rois;
    
    % Backward compatibility
    quality_metrics.fraction_good_rois = getFieldSafe(baseline_stats, 'fraction_good_rois', 0);
    quality_metrics.num_low_quality = 0;
end

function value = getFieldSafe(structure, fieldName, defaultValue)
    % Safely get field value with default fallback - handles nested fields
    
    if contains(fieldName, '.')
        % Handle nested field access (e.g., 'summary.active_rois')
        field_parts = split(fieldName, '.');
        current_struct = structure;
        
        for i = 1:length(field_parts)
            if isfield(current_struct, field_parts{i})
                current_struct = current_struct.(field_parts{i});
            else
                value = defaultValue;
                return;
            end
        end
        value = current_struct;
    else
        % Handle simple field access
        if isfield(structure, fieldName)
            value = structure.(fieldName);
        else
            value = defaultValue;
        end
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
        warning('Failed to save plots: %s', ME.message);
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
        
        % NEW: Detector tracking
        detectors_used = cell(1, length(fileResults));
        for i = 1:length(fileResults)
            if ~isempty(fileResults{i}) && isfield(fileResults{i}, 'detectorUsed')
                detectors_used{i} = fileResults{i}.detectorUsed;
            else
                detectors_used{i} = 'unknown';
            end
        end

        corrected_count = sum(strcmp(detectors_used, 'corrected_schmitt'));
        legacy_count = sum(strcmp(detectors_used, 'legacy_schmitt'));

        
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
        
        % NEW: Detector usage tracking
        summary.detection = struct();
        summary.detection.correctedSchmittFiles = corrected_count;
        summary.detection.legacySchmittFiles = legacy_count;
        if corrected_count > legacy_count
            summary.detection.primaryDetector = 'corrected_schmitt';
        else
            summary.detection.primaryDetector = 'legacy_schmitt';
        end
        
    catch ME
        warning('Error creating detailed summary: %s', ME.message);
        summary = createMinimalSummary(fileResults, metadataArray, loadTime, processTime);
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
    
    % NEW: Detector usage summary
    if isfield(summary, 'detection')
        fprintf('\n=== Event Detection Method ===\n');
        fprintf('Corrected Schmitt: %d files (includes outliers in noise)\n', summary.detection.correctedSchmittFiles);
        fprintf('Legacy Schmitt: %d files (excludes outliers from noise)\n', summary.detection.legacySchmittFiles);
        fprintf('Primary method: %s\n', summary.detection.primaryDetector);
    end
    
    % Timing summary
    if isfield(summary, 'timing')
        fprintf('\n=== Timing ===\n');
        fprintf('Total: %.3f s load + %.3f s process = %.3f s total\n', ...
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