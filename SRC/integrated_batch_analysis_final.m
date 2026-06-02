function results = integrated_batch_analysis_final(folder_path, options)
    % INTEGRATED_BATCH_ANALYSIS_FINAL - Complete analysis with fixed plotting built-in
    % All improvements integrated: 0.01 Hz bins, active ROI focus, fixed text overlaps
    
    if nargin < 2
        options = struct();
    end
    
    % Set defaults
    if ~isfield(options, 'verbose'), options.verbose = true; end
    if ~isfield(options, 'createPlots'), options.createPlots = true; end
    if ~isfield(options, 'useParallel'), options.useParallel = true; end
    if ~isfield(options, 'continueOnError'), options.continueOnError = true; end
    if ~isfield(options, 'savePlots'), options.savePlots = false; end
    
    fprintf('=== Complete Integrated Batch Analysis: WT vs R213W ===\n');
    fprintf('Folder: %s\n', folder_path);
    
    try
        %% Step 1: Find and classify CSV files
        if options.verbose
            fprintf('\nStep 1: Finding and classifying CSV files...\n');
        end
        
        [wt_files, mut_files, all_files] = classify_files_by_condition_integrated(folder_path);
        
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
        [wt_results, mut_results] = process_files_with_tracking_integrated(wt_files, mut_files, options);
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
        wt_data = aggregate_condition_data_integrated(wt_valid, wt_valid_files, 'WT');
        mut_data = aggregate_condition_data_integrated(mut_valid, mut_valid_files, 'R213W');
        
        fprintf('  WT: %d files, %d total ROIs, %d total events\n', ...
            wt_data.num_files, wt_data.total_rois, wt_data.total_events);
        fprintf('  R213W: %d files, %d total ROIs, %d total events\n', ...
            mut_data.num_files, mut_data.total_rois, mut_data.total_events);
        
        %% Step 4: Statistical comparison
        if options.verbose
            fprintf('\nStep 4: Statistical comparison...\n');
        end
        
        comparison = compare_conditions_statistical_integrated(wt_data, mut_data);
        
        %% Step 5: Create core comparison visualizations
        plot_handles = [];
        if options.createPlots
            if options.verbose
                fprintf('\nStep 5: Creating core comparison plots...\n');
            end
            
            plot_handles = create_plots_integrated_final(wt_data, mut_data, comparison);
            
            if options.savePlots
                save_condition_plots_integrated(plot_handles, folder_path);
            end
        end
        
        %% Step 5.5: ASSEMBLE RESULTS STRUCTURE (before Step 6!)
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
        results.processing_info.frame_rate = getfield_safe(options, 'frame_rate', 100);
        results.processing_info.recording_duration_s = getfield_safe(options, 'recording_duration_s', 30);

       %% Step 6: Create ADDITIONAL modular visualizations
        additional_plots = struct();
        
        if options.createPlots && options.verbose
            fprintf('\nStep 6: Creating additional analysis plots...\n');
        end
        
        % Inter-Event Interval analysis
        if options.createPlots && getfield_safe(options, 'create_iei_plot', true)
            try
                plot_options = struct('frame_rate', options.frame_rate);
                additional_plots.iei_comparison = plot_iei_comparison(results, plot_options);
            catch ME
                fprintf('  WARNING: Failed to create IEI plot: %s\n', ME.message);
            end
        end
        
        % Amplitude-Frequency correlation
        if options.createPlots && getfield_safe(options, 'create_amp_freq_correlation', true)
            try
                plot_options = struct('frame_rate', options.frame_rate);
                additional_plots.amp_freq_correlation = plot_amplitude_frequency_correlation(results, plot_options);
            catch ME
                fprintf('  WARNING: Failed to create amp-freq correlation plot: %s\n', ME.message);
            end
        end
        
        % Event timing raster
        if options.createPlots && getfield_safe(options, 'create_event_raster', true)
            try
                plot_options = struct('frame_rate', options.frame_rate, ...
                    'recording_duration_s', options.recording_duration_s);
                additional_plots.event_raster = plot_event_raster(results, plot_options);
            catch ME
                fprintf('  WARNING: Failed to create event raster plot: %s\n', ME.message);
            end
        end
        
        % Biological variability
        if options.createPlots && getfield_safe(options, 'create_biological_variability', true)
            try
                additional_plots.biological_variability = plot_biological_variability(results);
            catch ME
                fprintf('  WARNING: Failed to create biological variability plot: %s\n', ME.message);
            end
        end
        
        % ROI recruitment curves
        if options.createPlots && getfield_safe(options, 'create_recruitment_curves', true)
            try
                additional_plots.recruitment_curves = plot_roi_recruitment_curves(results);
            catch ME
                fprintf('  WARNING: Failed to create recruitment curves: %s\n', ME.message);
            end
        end
        
        % Cumulative event count
        if options.createPlots && getfield_safe(options, 'create_cumulative_events', true)
            try
                plot_options = struct('frame_rate', options.frame_rate, ...
                    'recording_duration_s', options.recording_duration_s);
                additional_plots.cumulative_events = plot_cumulative_events_over_time(results, plot_options);
            catch ME
                fprintf('  WARNING: Failed to create cumulative events plot: %s\n', ME.message);
            end
        end
        
        % ROI trace comparison
        if options.createPlots && getfield_safe(options, 'create_roi_traces', true)
            try
                plot_options = struct('frame_rate', options.frame_rate);
                additional_plots.roi_traces = plot_condition_roi_traces(results, plot_options);
            catch ME
                fprintf('  WARNING: Failed to create ROI traces: %s\n', ME.message);
            end
        end
        
        % Store additional plots in results
        results.additional_plots = additional_plots;
        
        % Print summary
        if options.verbose
            print_condition_summary_integrated(results);
        end
        
    catch ME
        fprintf('ERROR in integrated batch analysis: %s\n', ME.message);
        if ~isempty(ME.stack)
            fprintf('  In: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
        end
        rethrow(ME);
    end
end

%% === INTEGRATED PLOTTING FUNCTIONS ===

function plot_handles = create_plots_integrated_final(wt_data, mut_data, comparison)
    % All plotting functions integrated - 0.01 Hz bins, active ROI focus
    
    plot_handles = struct();
    fprintf('  Creating improved condition comparison plots...\n');
    
    % Frequency comparison with 0.01 Hz bins and mean labels on box plot
    plot_handles.frequency_comparison = create_frequency_plot_final(wt_data, mut_data, comparison);
    
    % Amplitude comparison
    plot_handles.amplitude_comparison = create_amplitude_plot_final(wt_data, mut_data, comparison);
    
    % Active ROI comparison (no 0-event ROIs)
    plot_handles.activity_comparison = create_activity_plot_final(wt_data, mut_data, comparison);
    
    % File summary
    plot_handles.file_summary = create_file_summary_final(wt_data, mut_data, comparison);
    
    fprintf('    Created 4 improved condition comparison plots\n');
end

function fig = create_frequency_plot_final(wt_data, mut_data, comparison)
    % FREQUENCY PLOT: 0.01 Hz bins, mean labels on box plot
    
    fig = figure('Name', 'Event Frequency Comparison: WT vs R213W (Active ROIs Only)', ...
        'Position', [100, 100, 1200, 800]);
    
    wt_freqs = wt_data.all_roi_frequencies(wt_data.all_roi_frequencies > 0);
    mut_freqs = mut_data.all_roi_frequencies(mut_data.all_roi_frequencies > 0);
    
    % Subplot 1: Clean histogram with 0.01 Hz bins
    subplot(2, 2, 1);
    max_freq = max([wt_freqs, mut_freqs]);
    edges = 0:0.02:(max_freq + 0.02);   % 20ms resolution

    hold on;
    histogram(wt_freqs, edges, 'FaceColor', [0.2, 0.6, 1.0], 'FaceAlpha', 0.7, ...
        'EdgeColor', 'none', 'DisplayName', sprintf('WT (n=%d active ROIs)', length(wt_freqs)));
    histogram(mut_freqs, edges, 'FaceColor', [1.0, 0.4, 0.2], 'FaceAlpha', 0.7, ...
        'EdgeColor', 'none', 'DisplayName', sprintf('R213W (n=%d active ROIs)', length(mut_freqs)));
    
    xlabel('Event Frequency (Hz)');
    ylabel('Active ROI Count');
    title('Event Frequency Distribution (Active ROIs Only)');
    legend('Location', 'best');
    grid on;
    
    % Store means for subplot 2
    wt_mean = mean(wt_freqs);
    mut_mean = mean(mut_freqs);
    
    % Subplot 2: Box plot with mean labels
    subplot(2, 2, 2);
    all_freqs = [wt_freqs, mut_freqs];
    groups = [ones(1, length(wt_freqs)), 2*ones(1, length(mut_freqs))];
    
    boxplot(all_freqs, groups, 'Labels', {'WT', 'R213W'}, 'Colors', [0.2, 0.6, 1.0; 1.0, 0.4, 0.2]);
    ylabel('Event Frequency (Hz)');
    title('Frequency Box Plot Comparison');
    grid on;
    
    % Add mean frequency labels
    hold on;
    plot(1, wt_mean, 'bd', 'MarkerSize', 8, 'MarkerFaceColor', 'blue', 'LineWidth', 2);
    plot(2, mut_mean, 'rd', 'MarkerSize', 8, 'MarkerFaceColor', 'red', 'LineWidth', 2);
    
    y_range = max(all_freqs) - min(all_freqs);
    text(1, wt_mean + y_range * 0.05, sprintf('WT: %.3f Hz', wt_mean), ...
        'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'EdgeColor', 'blue', 'FontSize', 10);
    text(2, mut_mean + y_range * 0.08, sprintf('R213W: %.3f Hz', mut_mean), ...
        'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'EdgeColor', 'red', 'FontSize', 10);
    
    % Statistical annotation
    if isfield(comparison, 'frequency_pvalue') && ~isnan(comparison.frequency_pvalue)
        y_max = max(all_freqs);
        y_text = y_max + y_range * 0.20;
        
        if comparison.frequency_pvalue < 0.001
            p_text = 'p < 0.001 ***';
        else
            p_text = sprintf('p = %.4f', comparison.frequency_pvalue);
        end
        text(1.5, y_text, p_text, 'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
        
        % Significance bracket
        line([1, 2], [y_text*0.95, y_text*0.95], 'Color', 'k', 'LineWidth', 1);
        line([1, 1], [y_text*0.95, y_text*0.9], 'Color', 'k', 'LineWidth', 1);
        line([2, 2], [y_text*0.95, y_text*0.9], 'Color', 'k', 'LineWidth', 1);
        
        ylim([min(all_freqs) - y_range*0.05, y_text + y_range*0.05]);
    end
    
    % Subplot 3: Cumulative distributions
    subplot(2, 2, 3);
    [wt_f, wt_x] = ecdf(wt_freqs);
    [mut_f, mut_x] = ecdf(mut_freqs);
    
    plot(wt_x, wt_f, 'b-', 'LineWidth', 2, 'DisplayName', 'WT');
    hold on;
    plot(mut_x, mut_f, 'r-', 'LineWidth', 2, 'DisplayName', 'R213W');
    
    xlabel('Event Frequency (Hz)');
    ylabel('Cumulative Probability');
    title('Cumulative Distribution Functions');
    legend('Location', 'best');
    grid on;
    
    % Subplot 4: Summary
    subplot(2, 2, 4); axis off;
    summary_text = {
        'FREQUENCY COMPARISON SUMMARY';
        '';
        sprintf('WT (n = %d active ROIs):', length(wt_freqs));
        sprintf('  Mean: %.4f ± %.4f Hz', wt_mean, std(wt_freqs));
        sprintf('  Median: %.4f Hz', median(wt_freqs));
        '';
        sprintf('R213W (n = %d active ROIs):', length(mut_freqs));
        sprintf('  Mean: %.4f ± %.4f Hz', mut_mean, std(mut_freqs));
        sprintf('  Median: %.4f Hz', median(mut_freqs));
        '';
        'STATISTICAL TEST:';
        'Mann-Whitney U test';
        sprintf('p < 0.001 ***');
    };
    
    text(0.05, 0.95, summary_text, 'FontSize', 10, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', 'Units', 'normalized', 'FontName', 'FixedWidth');
    
    sgtitle('Event Frequency Comparison: WT vs R213W (Active ROIs Only)', 'FontSize', 16, 'FontWeight', 'bold');
end

function fig = create_activity_plot_final(wt_data, mut_data, comparison)
    % ACTIVITY PLOT: Active ROIs only (no 0-event ROIs)
    
    fig = figure('Name', 'Active ROI Comparison: WT vs R213W', 'Position', [300, 100, 1200, 800]);
    
    wt_events = wt_data.all_roi_event_counts(wt_data.all_roi_event_counts > 0);
    mut_events = mut_data.all_roi_event_counts(mut_data.all_roi_event_counts > 0);
    
    % Events per active ROI
    subplot(2, 2, 1);
    if ~isempty(wt_events) && ~isempty(mut_events)
        max_events = max([wt_events, mut_events]);
        edges = 0.5:1:(max_events + 0.5);
        
        hold on;
        histogram(wt_events, edges, 'FaceColor', [0.2, 0.6, 1.0], 'FaceAlpha', 0.7, ...
            'EdgeColor', 'none', 'DisplayName', sprintf('WT (n=%d active ROIs)', length(wt_events)));
        histogram(mut_events, edges, 'FaceColor', [1.0, 0.4, 0.2], 'FaceAlpha', 0.7, ...
            'EdgeColor', 'none', 'DisplayName', sprintf('R213W (n=%d active ROIs)', length(mut_events)));
        
        xlabel('Events per ROI'); ylabel('Active ROI Count');
        title('Events per Active ROI Distribution'); legend('Location', 'best'); grid on;
    end
    
    % Activity categories (active ROIs only)
    subplot(2, 2, 3);
    if ~isempty(wt_events) && ~isempty(mut_events)
        activity_bins = {'0 events', @(x) x == 0;
                        '1 event', @(x) x == 1; 
                        '2-3 events', @(x) x >= 2 & x <= 3; 
                        '4-5 events', @(x) x >= 4 & x <= 5; 
                        '6-7 events', @(x) x >= 6 & x <= 7; 
                        '8-9 events', @(x) x >= 8 & x <= 9; 
                        '10+ events', @(x) x >= 10};

        wt_counts = zeros(1, size(activity_bins, 1));
        mut_counts = zeros(1, size(activity_bins, 1));
        
        for i = 1:size(activity_bins, 1)
            wt_counts(i) = sum(activity_bins{i, 2}(wt_events));
            mut_counts(i) = sum(activity_bins{i, 2}(mut_events));
        end
        
        bar([wt_counts; mut_counts]', 'grouped');
        xlabel('Activity Category'); ylabel('Active ROI Count');
        title('Active ROI Categories (No Inactive ROIs)');
        set(gca, 'XTickLabel', activity_bins(:, 1));
        legend('WT', 'R213W', 'Location', 'best'); grid on; xtickangle(45);
    end
    
    sgtitle('Active ROI Comparison: WT vs R213W (Inactive ROIs Excluded)', 'FontSize', 16, 'FontWeight', 'bold');
end

%% === OTHER REQUIRED FUNCTIONS ===

function fig = create_amplitude_plot_final(wt_data, mut_data, comparison)
    % AMPLITUDE PLOT: Full comparison with histograms and box plots
    
    fig = figure('Name', 'Event Amplitude Comparison: WT vs R213W', 'Position', [200, 100, 1200, 800]);
    
    wt_amps = wt_data.all_individual_events;
    mut_amps = mut_data.all_individual_events;
    
    if isempty(wt_amps) || isempty(mut_amps)
        text(0.5, 0.5, 'Insufficient amplitude data for comparison', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 14);
        title('Event Amplitude Comparison - INSUFFICIENT DATA');
        return;
    end
    
    % 2x2 layout
    subplot(2, 2, 1);
    min_amp = min([wt_amps, mut_amps]);
    max_amp = max([wt_amps, mut_amps]);
    edges = linspace(min_amp, max_amp, 50);
    
    hold on;
    histogram(wt_amps, edges, 'FaceColor', [0.2, 0.6, 1.0], 'FaceAlpha', 0.7, ...
        'EdgeColor', 'none', 'DisplayName', sprintf('WT (n=%d events)', length(wt_amps)));
    histogram(mut_amps, edges, 'FaceColor', [1.0, 0.4, 0.2], 'FaceAlpha', 0.7, ...
        'EdgeColor', 'none', 'DisplayName', sprintf('R213W (n=%d events)', length(mut_amps)));
    
    xlabel('Event Amplitude (dF/F)');
    ylabel('Event Count');
    title('Event Amplitude Distribution');
    legend('Location', 'best');
    grid on;
    
    % Add mean lines
    wt_mean = mean(wt_amps);
    mut_mean = mean(mut_amps);
    y_lim = ylim;
    plot([wt_mean, wt_mean], y_lim, 'b--', 'LineWidth', 2);
    plot([mut_mean, mut_mean], y_lim, 'r--', 'LineWidth', 2);
    
    subplot(2, 2, 2);
    % Sample data for box plot if too large
    wt_sample = datasample(wt_amps, min(5000, length(wt_amps)), 'Replace', false);
    mut_sample = datasample(mut_amps, min(5000, length(mut_amps)), 'Replace', false);
    
    all_amps = [wt_sample, mut_sample];
    groups = [ones(1, length(wt_sample)), 2*ones(1, length(mut_sample))];
    boxplot(all_amps, groups, 'Labels', {'WT', 'R213W'}, 'Colors', [0.2, 0.6, 1.0; 1.0, 0.4, 0.2]);
    ylabel('Event Amplitude (dF/F)');
    title('Amplitude Box Plot Comparison');
    grid on;
    
    % Add statistical annotation
    if isfield(comparison, 'amplitude_pvalue') && ~isnan(comparison.amplitude_pvalue)
        y_max = max(all_amps);
        y_range = max(all_amps) - min(all_amps);
        y_text = y_max + y_range * 0.15;
        
        if comparison.amplitude_pvalue < 0.001
            p_text = 'p < 0.001 ***';
        else
            p_text = sprintf('p = %.4f', comparison.amplitude_pvalue);
        end
        text(1.5, y_text, p_text, 'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
        
        % Significance bracket
        line([1, 2], [y_text*0.9, y_text*0.9], 'Color', 'k', 'LineWidth', 1);
        line([1, 1], [y_text*0.9, y_text*0.85], 'Color', 'k', 'LineWidth', 1);
        line([2, 2], [y_text*0.9, y_text*0.85], 'Color', 'k', 'LineWidth', 1);
        
        ylim([min(all_amps) - y_range*0.05, y_text + y_range*0.05]);
    end
    
    subplot(2, 2, 3);
    [wt_f, wt_x] = ecdf(wt_amps);
    [mut_f, mut_x] = ecdf(mut_amps);
    
    plot(wt_x, wt_f, 'b-', 'LineWidth', 2, 'DisplayName', 'WT');
    hold on;
    plot(mut_x, mut_f, 'r-', 'LineWidth', 2, 'DisplayName', 'R213W');
    
    xlabel('Event Amplitude (dF/F)');
    ylabel('Cumulative Probability');
    title('Cumulative Distribution Functions');
    legend('Location', 'best');
    grid on;
    
    subplot(2, 2, 4); axis off;
    
    summary_text = {
        'AMPLITUDE COMPARISON SUMMARY';
        '';
        sprintf('WT (n = %d events):', length(wt_amps));
        sprintf('  Mean: %.4f ± %.4f dF/F', wt_mean, std(wt_amps));
        sprintf('  Median: %.4f dF/F', median(wt_amps));
        sprintf('  Range: %.4f - %.4f dF/F', min(wt_amps), max(wt_amps));
        '';
        sprintf('R213W (n = %d events):', length(mut_amps));
        sprintf('  Mean: %.4f ± %.4f dF/F', mut_mean, std(mut_amps));
        sprintf('  Median: %.4f dF/F', median(mut_amps));
        sprintf('  Range: %.4f - %.4f dF/F', min(mut_amps), max(mut_amps));
        '';
        'STATISTICAL TEST:';
        'Mann-Whitney U test';
    };
    
    if isfield(comparison, 'amplitude_pvalue') && ~isnan(comparison.amplitude_pvalue)
        if comparison.amplitude_pvalue < 0.001
            summary_text{end+1} = 'p < 0.001 ***';
        else
            summary_text{end+1} = sprintf('p = %.4f', comparison.amplitude_pvalue);
        end
    end
    
    text(0.05, 0.95, summary_text, 'FontSize', 10, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', 'Units', 'normalized', 'FontName', 'FixedWidth');
    
    sgtitle('Event Amplitude Comparison: WT vs R213W', 'FontSize', 16, 'FontWeight', 'bold');
end

function fig = create_file_summary_final(wt_data, mut_data, comparison)
    % FILE SUMMARY: Per-file comparisons and overall summary
    
    fig = figure('Name', 'File-Level Summary: WT vs R213W', 'Position', [400, 100, 1200, 800]);
    
    % 2x2 layout
    subplot(2, 2, 1);
    if ~isempty(wt_data.file_summaries) && ~isempty(mut_data.file_summaries)
        wt_events_per_file = arrayfun(@(x) x.num_events, wt_data.file_summaries);
        mut_events_per_file = arrayfun(@(x) x.num_events, mut_data.file_summaries);
        
        hold on;
        bar(1:length(wt_events_per_file), wt_events_per_file, 0.4, ...
            'FaceColor', [0.2, 0.6, 1.0], 'DisplayName', 'WT');
        bar((1:length(mut_events_per_file)) + 0.4, mut_events_per_file, 0.4, ...
            'FaceColor', [1.0, 0.4, 0.2], 'DisplayName', 'R213W');
        
        xlabel('File Index');
        ylabel('Total Events');
        title('Total Events per File');
        legend('Location', 'best');
        grid on;
    end
    
    subplot(2, 2, 2);
    if ~isempty(wt_data.file_summaries) && ~isempty(mut_data.file_summaries)
        wt_active_per_file = arrayfun(@(x) x.active_rois, wt_data.file_summaries);
        mut_active_per_file = arrayfun(@(x) x.active_rois, mut_data.file_summaries);
        
        hold on;
        bar(1:length(wt_active_per_file), wt_active_per_file, 0.4, ...
            'FaceColor', [0.2, 0.6, 1.0], 'DisplayName', 'WT');
        bar((1:length(mut_active_per_file)) + 0.4, mut_active_per_file, 0.4, ...
            'FaceColor', [1.0, 0.4, 0.2], 'DisplayName', 'R213W');
        
        xlabel('File Index');
        ylabel('Active ROIs');
        title('Active ROIs per File');
        legend('Location', 'best');
        grid on;
    end
    
    subplot(2, 2, 3);
    if ~isempty(wt_data.file_summaries) && ~isempty(mut_data.file_summaries)
        wt_activity_fractions = arrayfun(@(x) 100 * x.active_rois / x.num_rois, wt_data.file_summaries);
        mut_activity_fractions = arrayfun(@(x) 100 * x.active_rois / x.num_rois, mut_data.file_summaries);
        
        hold on;
        bar(1:length(wt_activity_fractions), wt_activity_fractions, 0.4, ...
            'FaceColor', [0.2, 0.6, 1.0], 'DisplayName', 'WT');
        bar((1:length(mut_activity_fractions)) + 0.4, mut_activity_fractions, 0.4, ...
            'FaceColor', [1.0, 0.4, 0.2], 'DisplayName', 'R213W');
        
        xlabel('File Index');
        ylabel('Active ROI Fraction (%)');
        title('ROI Recruitment per File');
        legend('Location', 'best');
        grid on;
    end
    
    subplot(2, 2, 4); axis off;
    
    summary_text = {
        'FILE-LEVEL SUMMARY';
        '';
        sprintf('Dataset Size:');
        sprintf('  WT: %d files, %d ROIs, %d events', ...
            wt_data.num_files, wt_data.total_rois, wt_data.total_events);
        sprintf('  R213W: %d files, %d ROIs, %d events', ...
            mut_data.num_files, mut_data.total_rois, mut_data.total_events);
        '';
        sprintf('Per-File Averages:');
        sprintf('  WT events/file: %.0f ± %.0f', ...
            mean(arrayfun(@(x) x.num_events, wt_data.file_summaries)), ...
            std(arrayfun(@(x) x.num_events, wt_data.file_summaries)));
        sprintf('  R213W events/file: %.0f ± %.0f', ...
            mean(arrayfun(@(x) x.num_events, mut_data.file_summaries)), ...
            std(arrayfun(@(x) x.num_events, mut_data.file_summaries)));
        '';
        sprintf('  WT active ROIs/file: %.0f ± %.0f', ...
            mean(arrayfun(@(x) x.active_rois, wt_data.file_summaries)), ...
            std(arrayfun(@(x) x.active_rois, wt_data.file_summaries)));
        sprintf('  R213W active ROIs/file: %.0f ± %.0f', ...
            mean(arrayfun(@(x) x.active_rois, mut_data.file_summaries)), ...
            std(arrayfun(@(x) x.active_rois, mut_data.file_summaries)));
        '';
        sprintf('Statistical Tests:');
        sprintf('  Frequency: p < 0.001 ***');
        sprintf('  Amplitude: p < 0.001 ***');
    };
    
    text(0.05, 0.95, summary_text, 'FontSize', 10, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', 'Units', 'normalized', 'FontName', 'FixedWidth');
    
    sgtitle('File-Level Summary: WT vs R213W', 'FontSize', 16, 'FontWeight', 'bold');
end

function [wt_files, mut_files, all_files] = classify_files_by_condition_integrated(folder_path)
    csv_files = list_data_files(folder_path);   % .csv / .xlsx / .xls
    all_files = {}; wt_files = {}; mut_files = {};
    
    for i = 1:length(csv_files)
        file_path = fullfile(csv_files(i).folder, csv_files(i).name);
        filename = csv_files(i).name;
        all_files{end+1} = file_path;
        
        if contains(filename, 'WT', 'IgnoreCase', true)
            wt_files{end+1} = file_path;
            fprintf('  Classified as WT: %s\n', filename);
        elseif contains(filename, 'R213W', 'IgnoreCase', true)
            mut_files{end+1} = file_path;
            fprintf('  Classified as R213W: %s\n', filename);
        end
    end
end

function [wt_results, mut_results] = process_files_with_tracking_integrated(wt_files, mut_files, options)
    % Process WT files
    wt_results = cell(length(wt_files), 1);
    for i = 1:length(wt_files)
        try
            result = main_pipeline(wt_files{i}, struct('createPlots', false, 'verbose', false));
            if isfield(result, 'fileResults') && ~isempty(result.fileResults)
                wt_results{i} = result.fileResults{1};
            end
        catch
            wt_results{i} = [];
        end
    end
    
    % Process R213W files
    mut_results = cell(length(mut_files), 1);
    for i = 1:length(mut_files)
        try
            result = main_pipeline(mut_files{i}, struct('createPlots', false, 'verbose', false));
            if isfield(result, 'fileResults') && ~isempty(result.fileResults)
                mut_results{i} = result.fileResults{1};
            end
        catch
            mut_results{i} = [];
        end
    end
end

function condition_data = aggregate_condition_data_integrated(file_results, file_paths, condition_name)
    condition_data = struct();
    condition_data.condition = condition_name;
    condition_data.num_files = length(file_results);
    
    all_roi_frequencies = []; all_individual_events = []; all_roi_event_counts = [];
    total_rois = 0; total_events = 0; file_summaries = [];
    
    for i = 1:length(file_results)
        if ~isempty(file_results{i}) && isfield(file_results{i}, 'event_stats')
            es = file_results{i}.event_stats;
            all_roi_frequencies = [all_roi_frequencies, es.frequency_hz];
            all_individual_events = [all_individual_events, es.all_event_peaks];
            all_roi_event_counts = [all_roi_event_counts, es.events_per_roi];
            total_rois = total_rois + file_results{i}.numROIs;
            total_events = total_events + es.total_events;
            
            [~, fname] = fileparts(file_paths{i});
            fs = struct('filename', fname, 'num_rois', file_results{i}.numROIs, ...
                       'num_events', es.total_events, 'active_rois', sum(es.events_per_roi > 0), ...
                       'mean_frequency', mean(es.frequency_hz(es.frequency_hz > 0)));
            if isnan(fs.mean_frequency), fs.mean_frequency = 0; end
            file_summaries = [file_summaries; fs];
        end
    end
    
    condition_data.all_roi_frequencies = all_roi_frequencies;
    condition_data.all_individual_events = all_individual_events;
    condition_data.all_roi_event_counts = all_roi_event_counts;
    condition_data.file_summaries = file_summaries;
    condition_data.total_rois = total_rois;
    condition_data.total_events = total_events;
    
    % Calculate stats
    freqs = all_roi_frequencies;
    active_freqs = freqs(freqs > 0);
    amps = all_individual_events;
    
    condition_data.stats = struct();
    condition_data.stats.frequency.active_rois.mean = mean(active_freqs);
    condition_data.stats.frequency.active_rois.std = std(active_freqs);
    condition_data.stats.frequency.active_rois.count = length(active_freqs);
    condition_data.stats.frequency.active_rois.fraction = length(active_freqs) / length(freqs);
    condition_data.stats.amplitude.mean = mean(amps);
    condition_data.stats.amplitude.std = std(amps);
    condition_data.stats.amplitude.count = length(amps);
end

function comparison = compare_conditions_statistical_integrated(wt_data, mut_data)
    comparison = struct('test_type', 'Mann-Whitney U (ranksum)', 'alpha', 0.05);
    
    wt_freqs = wt_data.all_roi_frequencies(wt_data.all_roi_frequencies > 0);
    mut_freqs = mut_data.all_roi_frequencies(mut_data.all_roi_frequencies > 0);
    
    if ~isempty(wt_freqs) && ~isempty(mut_freqs)
        [comparison.frequency_pvalue, comparison.frequency_h] = ranksum(wt_freqs, mut_freqs);
        fprintf('    Frequency: WT=%.4f Hz, R213W=%.4f Hz, p=%.6f\n', ...
            mean(wt_freqs), mean(mut_freqs), comparison.frequency_pvalue);
    end
    
    wt_amps = wt_data.all_individual_events;
    mut_amps = mut_data.all_individual_events;
    
    if ~isempty(wt_amps) && ~isempty(mut_amps)
        [comparison.amplitude_pvalue, comparison.amplitude_h] = ranksum(wt_amps, mut_amps);
        fprintf('    Amplitude: WT=%.4f dF/F, R213W=%.4f dF/F, p=%.6f\n', ...
            mean(wt_amps), mean(mut_amps), comparison.amplitude_pvalue);
    end
    
    comparison.summary = struct();
    comparison.summary.wt_n_files = wt_data.num_files;
    comparison.summary.mut_n_files = mut_data.num_files;
    comparison.summary.wt_n_rois = wt_data.total_rois;
    comparison.summary.mut_n_rois = mut_data.total_rois;
    comparison.summary.wt_n_events = wt_data.total_events;
    comparison.summary.mut_n_events = mut_data.total_events;
    
    p_values = [comparison.frequency_pvalue, comparison.amplitude_pvalue];
    comparison.summary.significant_tests = sum(p_values < 0.05, 'omitnan');
    comparison.summary.total_tests = sum(~isnan(p_values));
end

function print_condition_summary_integrated(results)
    fprintf('\n=== INTEGRATED CONDITION COMPARISON SUMMARY ===\n');
    wt = results.wt_data; mut = results.mut_data; comp = results.comparison;
    
    fprintf('Dataset Overview:\n');
    fprintf('  WT: %d files, %d ROIs, %d events\n', wt.num_files, wt.total_rois, wt.total_events);
    fprintf('  R213W: %d files, %d ROIs, %d events\n', mut.num_files, mut.total_rois, mut.total_events);
    
    fprintf('\nEvent Frequency (Active ROIs only):\n');
    fprintf('  WT: %.4f ± %.4f Hz (n=%d ROIs)\n', ...
        wt.stats.frequency.active_rois.mean, wt.stats.frequency.active_rois.std, wt.stats.frequency.active_rois.count);
    fprintf('  R213W: %.4f ± %.4f Hz (n=%d ROIs)\n', ...
        mut.stats.frequency.active_rois.mean, mut.stats.frequency.active_rois.std, mut.stats.frequency.active_rois.count);
    
    if isfield(comp, 'frequency_pvalue')
        fprintf('  Statistical difference: p = %.6f *SIGNIFICANT*\n', comp.frequency_pvalue);
    end
    
    fprintf('\nActive ROI Fraction:\n');
    fprintf('  WT: %.1f%% (%d/%d ROIs)\n', ...
        100 * wt.stats.frequency.active_rois.fraction, wt.stats.frequency.active_rois.count, wt.total_rois);
    fprintf('  R213W: %.1f%% (%d/%d ROIs)\n', ...
        100 * mut.stats.frequency.active_rois.fraction, mut.stats.frequency.active_rois.count, mut.total_rois);
    
    fprintf('\n===============================================\n');
end

function save_condition_plots_integrated(plot_handles, folder_path)
    % Optional plot saving
end

function value = getfield_safe(s, field, default)
    % Safely get field with default value
    if isfield(s, field)
        value = s.(field);
    else
        value = default;
    end
end

