function results = batch_condition_analysis(folder_path, options)
    % BATCH_CONDITION_ANALYSIS - Compare WT vs R213W conditions across multiple files
    % 
    % Usage:
    %   results = batch_condition_analysis(folder_path);
    %   results = batch_condition_analysis(folder_path, options);
    %
    % Inputs:
    %   folder_path - Path to folder containing CSV files
    %   options     - Optional settings structure
    %
    % Outputs:
    %   results.wt_data      - Aggregated WT condition data
    %   results.mut_data     - Aggregated R213W mutant data  
    %   results.comparison   - Statistical comparison results
    %   results.file_results - Individual file processing results
    %   results.plot_handles - Condition comparison plots
    
    if nargin < 2
        options = struct();
    end
    
    % Set defaults
    if ~isfield(options, 'verbose'), options.verbose = true; end
    if ~isfield(options, 'createPlots'), options.createPlots = true; end
    if ~isfield(options, 'useParallel'), options.useParallel = true; end
    if ~isfield(options, 'continueOnError'), options.continueOnError = true; end
    if ~isfield(options, 'savePlots'), options.savePlots = false; end
    
    fprintf('=== Batch Condition Analysis: WT vs R213W ===\n');
    fprintf('Folder: %s\n', folder_path);
    
    try
        %% Step 1: Find and classify CSV files
        if options.verbose
            fprintf('\nStep 1: Finding and classifying CSV files...\n');
        end
        
        [wt_files, mut_files, all_files] = classify_files_by_condition(folder_path);
        
        fprintf('  Found %d WT files, %d R213W files (%d total)\n', ...
            length(wt_files), length(mut_files), length(all_files));
        
        if isempty(wt_files) || isempty(mut_files)
            error('Need both WT and R213W files for comparison. Found: %d WT, %d R213W', ...
                length(wt_files), length(mut_files));
        end
        
        %% Step 2: Process all files using existing pipeline
        if options.verbose
            fprintf('\nStep 2: Processing files with existing pipeline...\n');
        end
        
        tic;
        file_results = process_condition_files([wt_files; mut_files], options);
        process_time = toc;
        
        successful_files = sum(~cellfun(@isempty, file_results));
        fprintf('  Processed %d/%d files successfully in %.2f s\n', ...
            successful_files, length(all_files), process_time);
        
        %% Step 3: Aggregate data by condition
        if options.verbose
            fprintf('\nStep 3: Aggregating data by condition...\n');
        end
        
        % Separate results by condition
        wt_results = file_results(1:length(wt_files));
        mut_results = file_results(length(wt_files)+1:end);
        
        % Remove failed files
        wt_results = wt_results(~cellfun(@isempty, wt_results));
        mut_results = mut_results(~cellfun(@isempty, mut_results));
        
        % Aggregate condition data
        wt_data = aggregate_condition_data(wt_results, wt_files, 'WT');
        mut_data = aggregate_condition_data(mut_results, mut_files, 'R213W');
        
        fprintf('  WT: %d files, %d total ROIs, %d total events\n', ...
            wt_data.num_files, wt_data.total_rois, wt_data.total_events);
        fprintf('  R213W: %d files, %d total ROIs, %d total events\n', ...
            mut_data.num_files, mut_data.total_rois, mut_data.total_events);
        
        %% Step 4: Statistical comparison
        if options.verbose
            fprintf('\nStep 4: Statistical comparison...\n');
        end
        
        comparison = compare_conditions_statistical(wt_data, mut_data);
        
        %% Step 5: Create comparison visualizations
        plot_handles = [];
        if options.createPlots
            if options.verbose
                fprintf('\nStep 5: Creating comparison plots...\n');
            end
            
            plot_handles = create_condition_comparison_plots(wt_data, mut_data, comparison);
            
            if options.savePlots
                save_condition_plots(plot_handles, folder_path);
            end
        end
        
        %% Compile results
        results = struct();
        results.wt_data = wt_data;
        results.mut_data = mut_data;
        results.comparison = comparison;
        results.file_results = file_results;
        results.plot_handles = plot_handles;
        results.processing_info = struct();
        results.processing_info.folder_path = folder_path;
        results.processing_info.wt_files = wt_files;
        results.processing_info.mut_files = mut_files;
        results.processing_info.process_time = process_time;
        results.processing_info.timestamp = datetime('now');
        
        % Print summary
        if options.verbose
            print_condition_summary(results);
        end
        
    catch ME
        fprintf('ERROR in batch condition analysis: %s\n', ME.message);
        if ~isempty(ME.stack)
            fprintf('  In: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
        end
        rethrow(ME);
    end
end

function [wt_files, mut_files, all_files] = classify_files_by_condition(folder_path)
    % Classify CSV files as WT or R213W based on filename patterns
    
    % Find all CSV files
    if isfile(folder_path) && endsWith(folder_path, '.csv')
        % Single file mode - determine condition
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
    
    % Folder mode - find all CSV files
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
        
        % Classify based on filename patterns
        if contains(filename, 'Doc2b-WT', 'IgnoreCase', true) || ...
           contains(filename, 'WT', 'IgnoreCase', true)
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

function file_results = process_condition_files(file_list, options)
    % Process files using existing main_pipeline with IMPROVED error handling
    
    file_results = cell(length(file_list), 1);
    
    for i = 1:length(file_list)
        file_path = file_list{i};
        [~, filename, ~] = fileparts(file_path);
        
        try
            if options.verbose
                fprintf('  Processing: %s\n', filename);
            end
            
            % Use existing pipeline with modified options
            pipeline_options = options;
            pipeline_options.createPlots = false;  % Skip individual plots
            pipeline_options.verbose = false;      % Reduce verbosity
            
            result = main_pipeline(file_path, pipeline_options);
            
            % IMPROVED: Extract file result with better error checking
            if isfield(result, 'fileResults') && ~isempty(result.fileResults) && ...
               ~isempty(result.fileResults{1})
                
                file_results{i} = result.fileResults{1};  % Single file result
                
                % Verify essential fields are present
                if ~isfield(file_results{i}, 'event_stats') || ...
                   ~isfield(file_results{i}.event_stats, 'frequency_hz')
                    warning('File %s: Missing essential event_stats fields', filename);
                    file_results{i} = [];
                else
                    if options.verbose
                        fprintf('    SUCCESS: %d ROIs, %d events\n', ...
                            file_results{i}.numROIs, file_results{i}.event_stats.total_events);
                    end
                end
            else
                file_results{i} = [];
                if options.verbose
                    fprintf('    WARNING: No valid file result returned for: %s\n', filename);
                    
                    % Debug: Show what was actually returned
                    if isfield(result, 'fileResults')
                        if isempty(result.fileResults)
                            fprintf('      fileResults is empty\n');
                        elseif isempty(result.fileResults{1})
                            fprintf('      fileResults{1} is empty\n');
                        end
                    else
                        fprintf('      No fileResults field in result\n');
                        fprintf('      Available fields: %s\n', strjoin(fieldnames(result), ', '));
                    end
                end
            end
            
        catch ME
            if options.continueOnError
                fprintf('    ERROR: Failed to process %s: %s\n', filename, ME.message);
                if ~isempty(ME.stack) && options.verbose
                    fprintf('      In: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
                end
                file_results{i} = [];
            else
                rethrow(ME);
            end
        end
    end
end

function condition_data = aggregate_condition_data(file_results, file_paths, condition_name)
    % Aggregate event data across files for a single condition
    
    condition_data = struct();
    condition_data.condition = condition_name;
    condition_data.file_paths = file_paths;
    condition_data.num_files = length(file_results);
    
    % Initialize aggregated arrays
    all_roi_frequencies = [];           % All ROI frequencies across files
    all_roi_amplitudes = [];            % All ROI mean amplitudes across files  
    all_individual_events = [];         % All individual event peaks across files
    all_roi_event_counts = [];          % All ROI event counts across files
    
    file_summaries = [];
    total_rois = 0;
    total_events = 0;
    
    for i = 1:length(file_results)
        file_result = file_results{i};
        
        if isempty(file_result)
            continue;  % Skip failed files
        end
        
        % Extract event statistics
        if isfield(file_result, 'event_stats')
            event_stats = file_result.event_stats;
            
            % Aggregate per-ROI data
            if isfield(event_stats, 'frequency_hz')
                all_roi_frequencies = [all_roi_frequencies, event_stats.frequency_hz];
            end
            
            if isfield(event_stats, 'mean_peak_amplitude_per_roi')
                all_roi_amplitudes = [all_roi_amplitudes, event_stats.mean_peak_amplitude_per_roi];
            end
            
            if isfield(event_stats, 'all_event_peaks')
                all_individual_events = [all_individual_events, event_stats.all_event_peaks];
            end
            
            if isfield(event_stats, 'events_per_roi')
                all_roi_event_counts = [all_roi_event_counts, event_stats.events_per_roi];
            end
            
            % File summary
            [~, filename, ~] = fileparts(file_paths{i});
            file_summary = struct();
            file_summary.filename = filename;
            file_summary.num_rois = length(event_stats.frequency_hz);
            file_summary.num_events = length(event_stats.all_event_peaks);
            file_summary.active_rois = sum(event_stats.events_per_roi > 0);
            file_summary.mean_frequency = mean(event_stats.frequency_hz(event_stats.frequency_hz > 0));
            if isempty(file_summary.mean_frequency) || isnan(file_summary.mean_frequency)
                file_summary.mean_frequency = 0;
            end
            
            file_summaries = [file_summaries; file_summary];
            total_rois = total_rois + file_summary.num_rois;
            total_events = total_events + file_summary.num_events;
        end
    end
    
    % Store aggregated data
    condition_data.all_roi_frequencies = all_roi_frequencies;
    condition_data.all_roi_amplitudes = all_roi_amplitudes;
    condition_data.all_individual_events = all_individual_events;
    condition_data.all_roi_event_counts = all_roi_event_counts;
    condition_data.file_summaries = file_summaries;
    condition_data.total_rois = total_rois;
    condition_data.total_events = total_events;
    
    % Calculate condition-level statistics
    condition_data.stats = calculate_condition_statistics(condition_data);
end

function stats = calculate_condition_statistics(condition_data)
    % Calculate descriptive statistics for a condition
    
    stats = struct();
    
    % ROI-level frequency statistics (includes zeros)
    freqs = condition_data.all_roi_frequencies;
    stats.frequency.all_rois.mean = mean(freqs);
    stats.frequency.all_rois.median = median(freqs);
    stats.frequency.all_rois.std = std(freqs);
    stats.frequency.all_rois.count = length(freqs);
    
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
    stats.events_per_roi.mean = mean(counts);
    stats.events_per_roi.median = median(counts);
    stats.events_per_roi.std = std(counts);
    stats.events_per_roi.max = max(counts);
end

function print_condition_summary(results)
    % Print formatted summary of condition comparison
    
    fprintf('\n=== CONDITION COMPARISON SUMMARY ===\n');
    
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
    
    if isfield(comp, 'frequency_pvalue')
        fprintf('  Statistical difference: p = %.4f\n', comp.frequency_pvalue);
    end
    
    fprintf('\nEvent Amplitude (All Events):\n');
    fprintf('  WT: %.4f ± %.4f dF/F (n=%d events)\n', ...
        wt.stats.amplitude.mean, wt.stats.amplitude.std, wt.stats.amplitude.count);
    fprintf('  R213W: %.4f ± %.4f dF/F (n=%d events)\n', ...
        mut.stats.amplitude.mean, mut.stats.amplitude.std, mut.stats.amplitude.count);
    
    if isfield(comp, 'amplitude_pvalue')
        fprintf('  Statistical difference: p = %.4f\n', comp.amplitude_pvalue);
    end
    
    fprintf('\nActive ROI Fraction:\n');
    fprintf('  WT: %.1f%% (%d/%d ROIs)\n', ...
        100 * wt.stats.frequency.active_rois.fraction, wt.stats.frequency.active_rois.count, wt.total_rois);
    fprintf('  R213W: %.1f%% (%d/%d ROIs)\n', ...
        100 * mut.stats.frequency.active_rois.fraction, mut.stats.frequency.active_rois.count, mut.total_rois);
    
    fprintf('\n======================================\n');
end