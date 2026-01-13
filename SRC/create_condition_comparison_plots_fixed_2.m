function plot_handles = create_condition_comparison_plots_fixed(wt_data, mut_data, comparison)
    % CREATE_CONDITION_COMPARISON_PLOTS_FIXED - Generate WT vs R213W comparison visualizations
    % UPDATED: Active ROI focus, fine binning, fixed overlapping text
    
    plot_handles = struct();
    
    fprintf('  Creating improved condition comparison plots...\n');
    
    %% Plot 1: Event Frequency Comparison (IMPROVED)
    plot_handles.frequency_comparison = create_frequency_comparison_plot_fixed(wt_data, mut_data, comparison);
    
    %% Plot 2: Event Amplitude Comparison  
    plot_handles.amplitude_comparison = create_amplitude_comparison_plot_fixed(wt_data, mut_data, comparison);
    
    %% Plot 3: ROI Activity Comparison (ACTIVE ROIs ONLY)
    plot_handles.activity_comparison = create_activity_comparison_plot_fixed(wt_data, mut_data, comparison);
    
    %% Plot 4: File-Level Summary (FIXED overlapping text)
    plot_handles.file_summary = create_file_summary_plot_fixed(wt_data, mut_data, comparison);
    
    fprintf('    Created %d improved condition comparison plots\n', length(fieldnames(plot_handles)));
end

function fig = create_frequency_comparison_plot_fixed(wt_data, mut_data, comparison)
    % FIXED: Fine binning, active ROIs only, no overlapping text
    
    fig = figure('Name', 'Event Frequency Comparison: WT vs R213W (Active ROIs Only)', ...
        'Position', [100, 100, 1200, 800]);
    
    % Extract active ROI frequencies
    wt_freqs = wt_data.all_roi_frequencies(wt_data.all_roi_frequencies > 0);
    mut_freqs = mut_data.all_roi_frequencies(mut_data.all_roi_frequencies > 0);
    
    if isempty(wt_freqs) || isempty(mut_freqs)
        text(0.5, 0.5, 'Insufficient frequency data for comparison', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 14);
        title('Event Frequency Comparison - INSUFFICIENT DATA');
        return;
    end
    
    % 2x2 layout for frequency analysis
    
    % Subplot 1: Histograms overlay with FINE BINNING
    subplot(2, 2, 1);
    max_freq = max([wt_freqs, mut_freqs]);
    
    % Use fine binning for better resolution
    if max_freq <= 0.2
        edges = 0:0.005:max_freq + 0.005;  % 5ms resolution
    elseif max_freq <= 0.5
        edges = 0:0.01:max_freq + 0.01;   % 10ms resolution
    else
        edges = 0:0.02:max_freq + 0.02;   % 20ms resolution
    end
    
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
    
    % Add mean lines with labels
    wt_mean = mean(wt_freqs);
    mut_mean = mean(mut_freqs);
    y_lim = ylim;
    plot([wt_mean, wt_mean], y_lim, 'b--', 'LineWidth', 2);
    plot([mut_mean, mut_mean], y_lim, 'r--', 'LineWidth', 2);
    
    % Add mean labels
    text(wt_mean, y_lim(2)*0.9, sprintf('WT: %.3f Hz', wt_mean), ...
        'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'EdgeColor', 'blue');
    text(mut_mean, y_lim(2)*0.8, sprintf('R213W: %.3f Hz', mut_mean), ...
        'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'EdgeColor', 'red');
    
    % Subplot 2: Box plot comparison with MEAN LABELS MOVED HERE
    subplot(2, 2, 2);
    
    all_freqs = [wt_freqs, mut_freqs];
    groups = [ones(1, length(wt_freqs)), 2*ones(1, length(mut_freqs))];
    
    boxplot(all_freqs, groups, 'Labels', {'WT', 'R213W'}, 'Colors', [0.2, 0.6, 1.0; 1.0, 0.4, 0.2]);
    ylabel('Event Frequency (Hz)');
    title('Frequency Box Plot Comparison');
    grid on;
    
    % Add mean frequency labels on the box plot
    hold on;
    plot(1, wt_mean, 'bd', 'MarkerSize', 8, 'MarkerFaceColor', 'blue', 'LineWidth', 2);
    plot(2, mut_mean, 'rd', 'MarkerSize', 8, 'MarkerFaceColor', 'red', 'LineWidth', 2);
    
    % Add text labels for means
    text(1, wt_mean + (max(all_freqs) - min(all_freqs)) * 0.05, sprintf('WT: %.3f Hz', wt_mean), ...
        'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'EdgeColor', 'blue', 'FontSize', 10);
    text(2, mut_mean + (max(all_freqs) - min(all_freqs)) * 0.08, sprintf('R213W: %.3f Hz', mut_mean), ...
        'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'EdgeColor', 'red', 'FontSize', 10);
    
    % Add statistical annotation with PROPER positioning
    if isfield(comparison, 'frequency_pvalue') && ~isnan(comparison.frequency_pvalue)
        y_max = max(all_freqs);
        y_range = max(all_freqs) - min(all_freqs);
        y_text = y_max + y_range * 0.20;  % Position above the mean labels
        
        if comparison.frequency_pvalue < 0.001
            p_text = 'p < 0.001 ***';
        else
            p_text = sprintf('p = %.4f', comparison.frequency_pvalue);
        end
        text(1.5, y_text, p_text, 'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
        
        % Add significance bracket
        line([1, 2], [y_text*0.95, y_text*0.95], 'Color', 'k', 'LineWidth', 1);
        line([1, 1], [y_text*0.95, y_text*0.9], 'Color', 'k', 'LineWidth', 1);
        line([2, 2], [y_text*0.95, y_text*0.9], 'Color', 'k', 'LineWidth', 1);
        
        % Adjust y-limits to accommodate the annotation
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
    
    % Subplot 4: Summary statistics with NO overlapping text
    subplot(2, 2, 4); axis off;
    
    summary_text = {
        'FREQUENCY COMPARISON SUMMARY';
        '';
        sprintf('WT (n = %d active ROIs):', length(wt_freqs));
        sprintf('  Mean: %.4f ± %.4f Hz', wt_mean, std(wt_freqs));
        sprintf('  Median: %.4f Hz', median(wt_freqs));
        sprintf('  Range: %.4f - %.4f Hz', min(wt_freqs), max(wt_freqs));
        '';
        sprintf('R213W (n = %d active ROIs):', length(mut_freqs));
        sprintf('  Mean: %.4f ± %.4f Hz', mut_mean, std(mut_freqs));
        sprintf('  Median: %.4f Hz', median(mut_freqs));
        sprintf('  Range: %.4f - %.4f Hz', min(mut_freqs), max(mut_freqs));
        '';
        'STATISTICAL TEST:';
        sprintf('Mann-Whitney U test');
    };
    
    if isfield(comparison, 'frequency_pvalue') && ~isnan(comparison.frequency_pvalue)
        if comparison.frequency_pvalue < 0.001
            summary_text{end+1} = 'p < 0.001 ***';
        elseif comparison.frequency_pvalue < 0.01
            summary_text{end+1} = sprintf('p = %.4f **', comparison.frequency_pvalue);
        elseif comparison.frequency_pvalue < 0.05
            summary_text{end+1} = sprintf('p = %.4f *', comparison.frequency_pvalue);
        else
            summary_text{end+1} = sprintf('p = %.4f (ns)', comparison.frequency_pvalue);
        end
    else
        summary_text{end+1} = 'p = N/A';
    end
    
    text(0.05, 0.95, summary_text, 'FontSize', 10, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', 'Units', 'normalized', 'FontName', 'FixedWidth');
    
    sgtitle('Event Frequency Comparison: WT vs R213W (Active ROIs Only)', 'FontSize', 16, 'FontWeight', 'bold');
end

function fig = create_amplitude_comparison_plot_fixed(wt_data, mut_data, comparison)
    % Keep amplitude comparison as is - it was working correctly
    
    fig = figure('Name', 'Event Amplitude Comparison: WT vs R213W', ...
        'Position', [200, 100, 1200, 800]);
    
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
        'EdgeColor', 'none', 'DisplayName', sprintf('WT (n=%d)', length(wt_amps)));
    histogram(mut_amps, edges, 'FaceColor', [1.0, 0.4, 0.2], 'FaceAlpha', 0.7, ...
        'EdgeColor', 'none', 'DisplayName', sprintf('R213W (n=%d)', length(mut_amps)));
    
    xlabel('Event Amplitude (dF/F)');
    ylabel('Event Count');
    title('Event Amplitude Distribution');
    legend('Location', 'best');
    grid on;
    
    subplot(2, 2, 2);
    % Sample data for box plot
    wt_sample = datasample(wt_amps, min(5000, length(wt_amps)), 'Replace', false);
    mut_sample = datasample(mut_amps, min(5000, length(mut_amps)), 'Replace', false);
    
    all_amps = [wt_sample, mut_sample];
    groups = [ones(1, length(wt_sample)), 2*ones(1, length(mut_sample))];
    boxplot(all_amps, groups, 'Labels', {'WT', 'R213W'});
    ylabel('Event Amplitude (dF/F)');
    title('Amplitude Box Plot');
    grid on;
    
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
        sprintf('  Mean: %.4f ± %.4f dF/F', mean(wt_amps), std(wt_amps));
        sprintf('  Median: %.4f dF/F', median(wt_amps));
        '';
        sprintf('R213W (n = %d events):', length(mut_amps));
        sprintf('  Mean: %.4f ± %.4f dF/F', mean(mut_amps), std(mut_amps));
        sprintf('  Median: %.4f dF/F', median(mut_amps));
        '';
        'STATISTICAL TEST:';
        sprintf('Mann-Whitney U test');
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

function fig = create_activity_comparison_plot_fixed(wt_data, mut_data, comparison)
    % ACTIVE ROIs ONLY - completely remove 0-event ROIs
    
    fig = figure('Name', 'Active ROI Comparison: WT vs R213W (No Inactive ROIs)', ...
        'Position', [300, 100, 1200, 800]);
    
    % FILTER OUT ZERO-EVENT ROIs
    wt_events = wt_data.all_roi_event_counts(wt_data.all_roi_event_counts > 0);
    mut_events = mut_data.all_roi_event_counts(mut_data.all_roi_event_counts > 0);
    
    % 2x2 layout
    subplot(2, 2, 1);
    if ~isempty(wt_events) && ~isempty(mut_events)
        max_events = max([wt_events, mut_events]);
        edges = 0.5:1:(max_events + 0.5);  % Bins centered on integers
        
        hold on;
        histogram(wt_events, edges, 'FaceColor', [0.2, 0.6, 1.0], 'FaceAlpha', 0.7, ...
            'EdgeColor', 'none', 'DisplayName', sprintf('WT (n=%d active ROIs)', length(wt_events)));
        histogram(mut_events, edges, 'FaceColor', [1.0, 0.4, 0.2], 'FaceAlpha', 0.7, ...
            'EdgeColor', 'none', 'DisplayName', sprintf('R213W (n=%d active ROIs)', length(mut_events)));
        
        xlabel('Events per ROI');
        ylabel('Active ROI Count');
        title('Events per Active ROI Distribution');
        legend('Location', 'best');
        grid on;
    end
    
    subplot(2, 2, 2);
    % File-level activity fractions (what percentage of ROIs are active per file)
    if ~isempty(wt_data.file_summaries) && ~isempty(mut_data.file_summaries)
        wt_activity_fractions = arrayfun(@(x) x.active_rois / x.num_rois, wt_data.file_summaries) * 100;
        mut_activity_fractions = arrayfun(@(x) x.active_rois / x.num_rois, mut_data.file_summaries) * 100;
        
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
    
    subplot(2, 2, 3);
    % Activity categories (ACTIVE ROIs ONLY)
    if ~isempty(wt_events) && ~isempty(mut_events)
        activity_bins = {
            '1 event', @(x) x == 1;
            '2-3 events', @(x) x >= 2 & x <= 3;
            '4-5 events', @(x) x >= 4 & x <= 5;
            '6-10 events', @(x) x >= 6 & x <= 10;
            '11+ events', @(x) x >= 11
        };
        
        wt_counts = zeros(1, length(activity_bins));
        mut_counts = zeros(1, length(activity_bins));
        
        for i = 1:length(activity_bins)
            wt_counts(i) = sum(activity_bins{i, 2}(wt_events));
            mut_counts(i) = sum(activity_bins{i, 2}(mut_events));
        end
        
        bar_data = [wt_counts; mut_counts]';
        bar(bar_data, 'grouped');
        
        xlabel('Activity Category');
        ylabel('Active ROI Count');
        title('Active ROI Categories (No Inactive ROIs)');
        set(gca, 'XTickLabel', activity_bins(:, 1));
        legend('WT', 'R213W', 'Location', 'best');
        grid on;
        xtickangle(45);
    end
    
    subplot(2, 2, 4); axis off;
    
    summary_text = {
        'ACTIVE ROI COMPARISON';
        '';
        sprintf('WT Active ROIs:');
        sprintf('  Count: %d ROIs', length(wt_events));
        sprintf('  Total ROIs: %d', wt_data.total_rois);
        sprintf('  Activity Rate: %.1f%%', 100 * length(wt_events) / wt_data.total_rois);
        sprintf('  Events/Active ROI: %.1f ± %.1f', mean(wt_events), std(wt_events));
        '';
        sprintf('R213W Active ROIs:');
        sprintf('  Count: %d ROIs', length(mut_events));
        sprintf('  Total ROIs: %d', mut_data.total_rois);
        sprintf('  Activity Rate: %.1f%%', 100 * length(mut_events) / mut_data.total_rois);
        sprintf('  Events/Active ROI: %.1f ± %.1f', mean(mut_events), std(mut_events));
        '';
        sprintf('Recruitment Difference:');
        sprintf('  %.1f percentage points', ...
            100 * (length(wt_events) / wt_data.total_rois - length(mut_events) / mut_data.total_rois));
    };
    
    text(0.05, 0.95, summary_text, 'FontSize', 10, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', 'Units', 'normalized', 'FontName', 'FixedWidth');
    
    sgtitle('Active ROI Comparison: WT vs R213W (Inactive ROIs Excluded)', 'FontSize', 16, 'FontWeight', 'bold');
end

function fig = create_file_summary_plot_fixed(wt_data, mut_data, comparison)
    % FIXED: No overlapping text in summary
    
    fig = figure('Name', 'File-Level Summary: WT vs R213W', ...
        'Position', [400, 100, 1200, 800]);
    
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
        wt_rois_per_file = arrayfun(@(x) x.num_rois, wt_data.file_summaries);
        mut_rois_per_file = arrayfun(@(x) x.num_rois, mut_data.file_summaries);
        
        hold on;
        bar(1:length(wt_rois_per_file), wt_rois_per_file, 0.4, ...
            'FaceColor', [0.2, 0.6, 1.0], 'DisplayName', 'WT');
        bar((1:length(mut_rois_per_file)) + 0.4, mut_rois_per_file, 0.4, ...
            'FaceColor', [1.0, 0.4, 0.2], 'DisplayName', 'R213W');
        
        xlabel('File Index');
        ylabel('Total ROIs');
        title('ROIs per File');
        legend('Location', 'best');
        grid on;
    end
    
    subplot(2, 2, 3);
    if ~isempty(wt_data.file_summaries) && ~isempty(mut_data.file_summaries)
        wt_mean_freq = arrayfun(@(x) x.mean_frequency, wt_data.file_summaries);
        mut_mean_freq = arrayfun(@(x) x.mean_frequency, mut_data.file_summaries);
        
        wt_mean_freq(isnan(wt_mean_freq)) = 0;
        mut_mean_freq(isnan(mut_mean_freq)) = 0;
        
        hold on;
        bar(1:length(wt_mean_freq), wt_mean_freq, 0.4, ...
            'FaceColor', [0.2, 0.6, 1.0], 'DisplayName', 'WT');
        bar((1:length(mut_mean_freq)) + 0.4, mut_mean_freq, 0.4, ...
            'FaceColor', [1.0, 0.4, 0.2], 'DisplayName', 'R213W');
        
        xlabel('File Index');
        ylabel('Mean Frequency (Hz)');
        title('Mean Event Frequency per File');
        legend('Location', 'best');
        grid on;
    end
    
    % Subplot 4: FIXED summary with proper spacing
    subplot(2, 2, 4); axis off;
    
    summary_text = {
        'OVERALL COMPARISON SUMMARY';
        '';
        sprintf('Dataset Size:');
        sprintf('  WT: %d files, %d ROIs, %d events', ...
            wt_data.num_files, wt_data.total_rois, wt_data.total_events);
        sprintf('  R213W: %d files, %d ROIs, %d events', ...
            mut_data.num_files, mut_data.total_rois, mut_data.total_events);
        '';
        sprintf('Key Metrics:');
        sprintf('  Event frequency (active ROIs):');
        sprintf('    WT: %.4f ± %.4f Hz', ...
            wt_data.stats.frequency.active_rois.mean, wt_data.stats.frequency.active_rois.std);
        sprintf('    R213W: %.4f ± %.4f Hz', ...
            mut_data.stats.frequency.active_rois.mean, mut_data.stats.frequency.active_rois.std);
        '';
        sprintf('  Event amplitude (all events):');
        sprintf('    WT: %.4f ± %.4f dF/F', ...
            wt_data.stats.amplitude.mean, wt_data.stats.amplitude.std);
        sprintf('    R213W: %.4f ± %.4f dF/F', ...
            mut_data.stats.amplitude.mean, mut_data.stats.amplitude.std);
        '';
        sprintf('Statistical Tests:');
    };
    
    % Add significance summary with PROPER spacing
    if isfield(comparison, 'summary') && isfield(comparison.summary, 'significant_tests')
        summary_text{end+1} = sprintf('  %d/%d tests significant (p < 0.05)', ...
            comparison.summary.significant_tests, comparison.summary.total_tests);
    else
        sig_count = 0;
        total_count = 0;
        if isfield(comparison, 'frequency_pvalue') && ~isnan(comparison.frequency_pvalue)
            total_count = total_count + 1;
            if comparison.frequency_pvalue < 0.05
                sig_count = sig_count + 1;
            end
        end
        if isfield(comparison, 'amplitude_pvalue') && ~isnan(comparison.amplitude_pvalue)
            total_count = total_count + 1;
            if comparison.amplitude_pvalue < 0.05
                sig_count = sig_count + 1;
            end
        end
        summary_text{end+1} = sprintf('  %d/%d tests significant (p < 0.05)', sig_count, total_count);
    end
    
    text(0.05, 0.95, summary_text, 'FontSize', 9, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', 'Units', 'normalized', 'FontName', 'FixedWidth');
    
    sgtitle('File-Level Summary: WT vs R213W', 'FontSize', 16, 'FontWeight', 'bold');
end

%% APPLY THE FIXES TO YOUR EXISTING RESULTS
fprintf('\n=== Creating Fixed Plots ===\n');

% Check if results exist
if exist('results', 'var') && isfield(results, 'wt_data')
    fprintf('Applying fixes to existing results...\n');
    
    % Create new improved plots
    close all;  % Close old plots
    new_plot_handles = create_condition_comparison_plots_fixed(results.wt_data, results.mut_data, results.comparison);
    
    % Update the results structure
    results.plot_handles = new_plot_handles;
    
    fprintf('✓ Generated 4 fixed plots with:\n');
    fprintf('  - Fine frequency binning (5ms resolution)\n');
    fprintf('  - Active ROI focus (no 0-event ROIs)\n');
    fprintf('  - Fixed overlapping text issues\n');
    fprintf('  - Proper statistical annotation positioning\n');
    
else
    fprintf('No results found. Please run the batch analysis first:\n');
    fprintf('  results = fixed_batch_condition_analysis(data_folder);\n');
    fprintf('Then run this script again.\n');
end