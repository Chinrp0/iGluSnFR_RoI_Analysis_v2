function results = fixed_batch_condition_analysis(folder_path, options)
    % FIXED_BATCH_CONDITION_ANALYSIS - Corrected version with proper file tracking
    % 
    % Fixes the issue where R213W files were processed but lost during aggregation
    
    if nargin < 2
        options = struct();
    end
    
    % Set defaults
    if ~isfield(options, 'verbose'), options.verbose = true; end
    if ~isfield(options, 'createPlots'), options.createPlots = true; end
    if ~isfield(options, 'useParallel'), options.useParallel = true; end
    if ~isfield(options, 'continueOnError'), options.continueOnError = true; end
    if ~isfield(options, 'savePlots'), options.savePlots = false; end
    
    fprintf('=== FIXED Batch Condition Analysis: WT vs R213W ===\n');
    fprintf('Folder: %s\n', folder_path);
    
    try
        %% Step 1: Find and classify CSV files
        if options.verbose
            fprintf('\nStep 1: Finding and classifying CSV files...\n');
        end
        
        [wt_files, mut_files, all_files] = classify_files_by_condition_fixed(folder_path);
        
        fprintf('  Found %d WT files, %d R213W files (%d total)\n', ...
            length(wt_files), length(mut_files), length(all_files));
        
        if isempty(wt_files) || isempty(mut_files)
            error('Need both WT and R213W files for comparison. Found: %d WT, %d R213W', ...
                length(wt_files), length(mut_files));
        end
        
        %% Step 2: Process files with proper tracking
        if options.verbose
            fprintf('\nStep 2: Processing files with condition tracking...\n');
        end
        
        tic;
        [wt_results, mut_results] = process_files_with_tracking(wt_files, mut_files, options);
        process_time = toc;
        
        wt_success = sum(~cellfun(@isempty, wt_results));
        mut_success = sum(~cellfun(@isempty, mut_results));
        
        fprintf('  Processed WT: %d/%d files successfully\n', wt_success, length(wt_files));
        fprintf('  Processed R213W: %d/%d files successfully\n', mut_success, length(mut_files));
        fprintf('  Total processing time: %.3f s\n', process_time);
        
        %% Step 3: Aggregate data by condition
        if options.verbose
            fprintf('\nStep 3: Aggregating data by condition...\n');
        end
        
        % Remove failed files
        wt_valid = wt_results(~cellfun(@isempty, wt_results));
        mut_valid = mut_results(~cellfun(@isempty, mut_results));
        
        fprintf('  Valid results - WT: %d files, R213W: %d files\n', ...
            length(wt_valid), length(mut_valid));
        
        % Get corresponding file names for valid results
        wt_valid_files = wt_files(~cellfun(@isempty, wt_results));
        mut_valid_files = mut_files(~cellfun(@isempty, mut_results));
        
        % Aggregate condition data
        wt_data = aggregate_condition_data_fixed(wt_valid, wt_valid_files, 'WT');
        mut_data = aggregate_condition_data_fixed(mut_valid, mut_valid_files, 'R213W');
        
        fprintf('  WT: %d files, %d total ROIs, %d total events\n', ...
            wt_data.num_files, wt_data.total_rois, wt_data.total_events);
        fprintf('  R213W: %d files, %d total ROIs, %d total events\n', ...
            mut_data.num_files, mut_data.total_rois, mut_data.total_events);
        
        %% Step 4: Statistical comparison
        if options.verbose
            fprintf('\nStep 4: Statistical comparison...\n');
        end
        
        comparison = compare_conditions_statistical_fixed(wt_data, mut_data);
        
        %% Step 5: Create comparison visualizations
        plot_handles = [];
        if options.createPlots
            if options.verbose
                fprintf('\nStep 5: Creating comparison plots...\n');
            end
            
            plot_handles = create_condition_comparison_plots_fixed(wt_data, mut_data, comparison);
            
            if options.savePlots
                save_condition_plots(plot_handles, folder_path);
            end
        end
        
        %% Compile results
        results = struct();
        results.wt_data = wt_data;
        results.mut_data = mut_data;
        results.comparison = comparison;
        results.file_results = struct('wt_results', {wt_results}, 'mut_results', {mut_results});
        results.plot_handles = plot_handles;
        results.processing_info = struct();
        results.processing_info.folder_path = folder_path;
        results.processing_info.wt_files = wt_files;
        results.processing_info.mut_files = mut_files;
        results.processing_info.process_time = process_time;
        results.processing_info.timestamp = datetime('now');
        
        % Print summary
        if options.verbose
            print_condition_summary_fixed(results);
        end
        
    catch ME
        fprintf('ERROR in fixed batch condition analysis: %s\n', ME.message);
        if ~isempty(ME.stack)
            fprintf('  In: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
        end
        rethrow(ME);
    end
end

function [wt_files, mut_files, all_files] = classify_files_by_condition_fixed(folder_path)
    % Fixed version with better debugging
    
    if isfile(folder_path) && endsWith(folder_path, '.csv')
        [~, filename, ~] = fileparts(folder_path);
        if contains(filename, 'WT', 'IgnoreCase', true)
            wt_files = {folder_path};
            mut_files = {};
        elseif contains(filename, 'R213W', 'IgnoreCase', true)
            wt_files = {};
            mut_files = {folder_path};
        else
            error('Cannot determine condition from filename: %s', filename);
        end
        all_files = {folder_path};
        return;
    end
    
    csv_files = dir(fullfile(folder_path, '*.csv'));
    
    if isempty(csv_files)
        error('No CSV files found in folder: %s', folder_path);
    end
    
    all_files = {};
    wt_files = {};
    mut_files = {};
    
    for i = 1:length(csv_files)
        file_path = fullfile(csv_files(i).folder, csv_files(i).name);
        filename = csv_files(i).name;
        all_files{end+1} = file_path;
        
        if contains(filename, 'Doc2b-WT', 'IgnoreCase', true) || ...
           (contains(filename, 'WT', 'IgnoreCase', true) && ~contains(filename, 'R213W', 'IgnoreCase', true))
            wt_files{end+1} = file_path;
            fprintf('  Classified as WT: %s\n', filename);
            
        elseif contains(filename, 'Doc2b-R213W', 'IgnoreCase', true) || ...
               contains(filename, 'R213W', 'IgnoreCase', true)
            mut_files{end+1} = file_path;
            fprintf('  Classified as R213W: %s\n', filename);
            
        else
            warning('Cannot classify file (skipping): %s', filename);
        end
    end
end

function [wt_results, mut_results] = process_files_with_tracking(wt_files, mut_files, options)
    % Process files separately by condition to avoid indexing issues
    
    % Process WT files
    fprintf('  Processing WT files...\n');
    wt_results = cell(length(wt_files), 1);
    
    for i = 1:length(wt_files)
        file_path = wt_files{i};
        [~, filename, ~] = fileparts(file_path);
        
        try
            if options.verbose
                fprintf('    WT %d/%d: %s\n', i, length(wt_files), filename);
            end
            
            pipeline_options = options;
            pipeline_options.createPlots = false;
            pipeline_options.verbose = false;
            
            result = main_pipeline(file_path, pipeline_options);
            
            if isfield(result, 'fileResults') && ~isempty(result.fileResults) && ...
               ~isempty(result.fileResults{1})
                wt_results{i} = result.fileResults{1};
                if options.verbose
                    fprintf('      SUCCESS: %d ROIs, %d events\n', ...
                        wt_results{i}.numROIs, wt_results{i}.event_stats.total_events);
                end
            else
                wt_results{i} = [];
                fprintf('      FAILED: No valid file result\n');
            end
            
        catch ME
            if options.continueOnError
                fprintf('    ERROR processing WT file %s: %s\n', filename, ME.message);
                wt_results{i} = [];
            else
                rethrow(ME);
            end
        end
    end
    
    % Process R213W files
    fprintf('  Processing R213W files...\n');
    mut_results = cell(length(mut_files), 1);
    
    for i = 1:length(mut_files)
        file_path = mut_files{i};
        [~, filename, ~] = fileparts(file_path);
        
        try
            if options.verbose
                fprintf('    R213W %d/%d: %s\n', i, length(mut_files), filename);
            end
            
            pipeline_options = options;
            pipeline_options.createPlots = false;
            pipeline_options.verbose = false;
            
            result = main_pipeline(file_path, pipeline_options);
            
            if isfield(result, 'fileResults') && ~isempty(result.fileResults) && ...
               ~isempty(result.fileResults{1})
                mut_results{i} = result.fileResults{1};
                if options.verbose
                    fprintf('      SUCCESS: %d ROIs, %d events\n', ...
                        mut_results{i}.numROIs, mut_results{i}.event_stats.total_events);
                end
            else
                mut_results{i} = [];
                fprintf('      FAILED: No valid file result\n');
            end
            
        catch ME
            if options.continueOnError
                fprintf('    ERROR processing R213W file %s: %s\n', filename, ME.message);
                mut_results{i} = [];
            else
                rethrow(ME);
            end
        end
    end
end

function condition_data = aggregate_condition_data_fixed(file_results, file_paths, condition_name)
    % Fixed aggregation with better error checking
    
    condition_data = struct();
    condition_data.condition = condition_name;
    condition_data.file_paths = file_paths;
    condition_data.num_files = length(file_results);
    
    fprintf('    Aggregating %s data from %d files...\n', condition_name, length(file_results));
    
    % Initialize aggregated arrays
    all_roi_frequencies = [];
    all_roi_amplitudes = [];
    all_individual_events = [];
    all_roi_event_counts = [];
    
    file_summaries = [];
    total_rois = 0;
    total_events = 0;
    
    for i = 1:length(file_results)
        file_result = file_results{i};
        
        if isempty(file_result)
            fprintf('      File %d: EMPTY - skipping\n', i);
            continue;
        end
        
        if ~isfield(file_result, 'event_stats')
            fprintf('      File %d: Missing event_stats - skipping\n', i);
            continue;
        end
        
        event_stats = file_result.event_stats;
        [~, filename, ~] = fileparts(file_paths{i});
        
        % Aggregate per-ROI data
        if isfield(event_stats, 'frequency_hz')
            all_roi_frequencies = [all_roi_frequencies, event_stats.frequency_hz];
            fprintf('      File %d (%s): Added %d ROI frequencies\n', i, filename, length(event_stats.frequency_hz));
        end
        
        if isfield(event_stats, 'mean_peak_amplitude_per_roi')
            all_roi_amplitudes = [all_roi_amplitudes, event_stats.mean_peak_amplitude_per_roi];
        end
        
        if isfield(event_stats, 'all_event_peaks')
            all_individual_events = [all_individual_events, event_stats.all_event_peaks];
            fprintf('      File %d (%s): Added %d individual events\n', i, filename, length(event_stats.all_event_peaks));
        end
        
        if isfield(event_stats, 'events_per_roi')
            all_roi_event_counts = [all_roi_event_counts, event_stats.events_per_roi];
        end
        
        % File summary
        file_summary = struct();
        file_summary.filename = filename;
        file_summary.num_rois = file_result.numROIs;
        file_summary.num_events = event_stats.total_events;
        file_summary.active_rois = sum(event_stats.events_per_roi > 0);
        file_summary.mean_frequency = mean(event_stats.frequency_hz(event_stats.frequency_hz > 0));
        if isempty(file_summary.mean_frequency) || isnan(file_summary.mean_frequency)
            file_summary.mean_frequency = 0;
        end
        
        file_summaries = [file_summaries; file_summary];
        total_rois = total_rois + file_summary.num_rois;
        total_events = total_events + file_summary.num_events;
    end
    
    % Store aggregated data
    condition_data.all_roi_frequencies = all_roi_frequencies;
    condition_data.all_roi_amplitudes = all_roi_amplitudes;
    condition_data.all_individual_events = all_individual_events;
    condition_data.all_roi_event_counts = all_roi_event_counts;
    condition_data.file_summaries = file_summaries;
    condition_data.total_rois = total_rois;
    condition_data.total_events = total_events;
    
    fprintf('    %s aggregation complete: %d ROIs, %d events\n', ...
        condition_name, total_rois, total_events);
    
    % Calculate condition-level statistics
    condition_data.stats = calculate_condition_statistics_fixed(condition_data);
end

function stats = calculate_condition_statistics_fixed(condition_data)
    % Fixed statistics calculation
    
    stats = struct();
    
    % ROI-level frequency statistics (includes zeros)
    freqs = condition_data.all_roi_frequencies;
    if ~isempty(freqs)
        stats.frequency.all_rois.mean = mean(freqs);
        stats.frequency.all_rois.median = median(freqs);
        stats.frequency.all_rois.std = std(freqs);
        stats.frequency.all_rois.count = length(freqs);
    else
        stats.frequency.all_rois = struct('mean', 0, 'median', 0, 'std', 0, 'count', 0);
    end
    
    % Active ROI frequency statistics (excludes zeros)
    active_freqs = freqs(freqs > 0);
    if ~isempty(active_freqs)
        stats.frequency.active_rois.mean = mean(active_freqs);
        stats.frequency.active_rois.median = median(active_freqs);
        stats.frequency.active_rois.std = std(active_freqs);
        stats.frequency.active_rois.count = length(active_freqs);
        stats.frequency.active_rois.fraction = length(active_freqs) / length(freqs);
    else
        stats.frequency.active_rois = struct('mean', 0, 'median', 0, 'std', 0, 'count', 0, 'fraction', 0);
    end
    
    % Event amplitude statistics
    amps = condition_data.all_individual_events;
    if ~isempty(amps)
        stats.amplitude.mean = mean(amps);
        stats.amplitude.median = median(amps);
        stats.amplitude.std = std(amps);
        stats.amplitude.count = length(amps);
        stats.amplitude.range = [min(amps), max(amps)];
    else
        stats.amplitude = struct('mean', 0, 'median', 0, 'std', 0, 'count', 0, 'range', [0, 0]);
    end
    
    % Event count statistics per ROI
    counts = condition_data.all_roi_event_counts;
    if ~isempty(counts)
        stats.events_per_roi.mean = mean(counts);
        stats.events_per_roi.median = median(counts);
        stats.events_per_roi.std = std(counts);
        stats.events_per_roi.max = max(counts);
    else
        stats.events_per_roi = struct('mean', 0, 'median', 0, 'std', 0, 'max', 0);
    end
end

function comparison = compare_conditions_statistical_fixed(wt_data, mut_data)
    % Fixed statistical comparison
    
    comparison = struct();
    comparison.test_type = 'Mann-Whitney U (ranksum)';
    comparison.alpha = 0.05;
    
    fprintf('    Running statistical comparisons...\n');
    
    % Event Frequency Comparison
    wt_freqs = wt_data.all_roi_frequencies(wt_data.all_roi_frequencies > 0);
    mut_freqs = mut_data.all_roi_frequencies(mut_data.all_roi_frequencies > 0);
    
    if ~isempty(wt_freqs) && ~isempty(mut_freqs)
        [comparison.frequency_pvalue, comparison.frequency_h] = ranksum(wt_freqs, mut_freqs);
        fprintf('    Frequency: WT=%.4f Hz, R213W=%.4f Hz, p=%.4f\n', ...
            mean(wt_freqs), mean(mut_freqs), comparison.frequency_pvalue);
    else
        comparison.frequency_pvalue = NaN;
        comparison.frequency_h = 0;
        fprintf('    Frequency: Insufficient data\n');
    end
    
    % Event Amplitude Comparison
    wt_amps = wt_data.all_individual_events;
    mut_amps = mut_data.all_individual_events;
    
    if ~isempty(wt_amps) && ~isempty(mut_amps)
        [comparison.amplitude_pvalue, comparison.amplitude_h] = ranksum(wt_amps, mut_amps);
        fprintf('    Amplitude: WT=%.4f dF/F, R213W=%.4f dF/F, p=%.4f\n', ...
            mean(wt_amps), mean(mut_amps), comparison.amplitude_pvalue);
    else
        comparison.amplitude_pvalue = NaN;
        comparison.amplitude_h = 0;
        fprintf('    Amplitude: Insufficient data\n');
    end
    
    % Summary with significance counting
    comparison.summary = struct();
    comparison.summary.wt_n_files = wt_data.num_files;
    comparison.summary.mut_n_files = mut_data.num_files;
    comparison.summary.wt_n_rois = wt_data.total_rois;
    comparison.summary.mut_n_rois = mut_data.total_rois;
    comparison.summary.wt_n_events = wt_data.total_events;
    comparison.summary.mut_n_events = mut_data.total_events;
    
    % Count significant tests
    p_values = [comparison.frequency_pvalue, comparison.amplitude_pvalue];
    comparison.summary.significant_tests = sum(p_values < 0.05, 'omitnan');
    comparison.summary.total_tests = sum(~isnan(p_values));
end

function print_condition_summary_fixed(results)
    % Fixed summary printing
    
    fprintf('\n=== FIXED CONDITION COMPARISON SUMMARY ===\n');
    
    wt = results.wt_data;
    mut = results.mut_data;
    comp = results.comparison;
    
    fprintf('Dataset Overview:\n');
    fprintf('  WT: %d files, %d ROIs, %d events\n', wt.num_files, wt.total_rois, wt.total_events);
    fprintf('  R213W: %d files, %d ROIs, %d events\n', mut.num_files, mut.total_rois, mut.total_events);
    
    fprintf('\nEvent Frequency (Active ROIs only):\n');
    fprintf('  WT: %.4f ± %.4f Hz (n=%d ROIs)\n', ...
        wt.stats.frequency.active_rois.mean, wt.stats.frequency.active_rois.std, wt.stats.frequency.active_rois.count);
    fprintf('  R213W: %.4f ± %.4f Hz (n=%d ROIs)\n', ...
        mut.stats.frequency.active_rois.mean, mut.stats.frequency.active_rois.std, mut.stats.frequency.active_rois.count);
    
    if isfield(comp, 'frequency_pvalue') && ~isnan(comp.frequency_pvalue)
        fprintf('  Statistical difference: p = %.4f', comp.frequency_pvalue);
        if comp.frequency_pvalue < 0.05
            fprintf(' *SIGNIFICANT*');
        end
        fprintf('\n');
    end
    
    fprintf('\nEvent Amplitude (All Events):\n');
    fprintf('  WT: %.4f ± %.4f dF/F (n=%d events)\n', ...
        wt.stats.amplitude.mean, wt.stats.amplitude.std, wt.stats.amplitude.count);
    fprintf('  R213W: %.4f ± %.4f dF/F (n=%d events)\n', ...
        mut.stats.amplitude.mean, mut.stats.amplitude.std, mut.stats.amplitude.count);
    
    if isfield(comp, 'amplitude_pvalue') && ~isnan(comp.amplitude_pvalue)
        fprintf('  Statistical difference: p = %.4f', comp.amplitude_pvalue);
        if comp.amplitude_pvalue < 0.05
            fprintf(' *SIGNIFICANT*');
        end
        fprintf('\n');
    end
    
    fprintf('\nActive ROI Fraction:\n');
    fprintf('  WT: %.1f%% (%d/%d ROIs)\n', ...
        100 * wt.stats.frequency.active_rois.fraction, wt.stats.frequency.active_rois.count, wt.total_rois);
    fprintf('  R213W: %.1f%% (%d/%d ROIs)\n', ...
        100 * mut.stats.frequency.active_rois.fraction, mut.stats.frequency.active_rois.count, mut.total_rois);
    
    fprintf('\n===============================================\n');
end