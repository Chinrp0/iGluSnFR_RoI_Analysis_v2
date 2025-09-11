function plot_handles = create_condition_comparison_plots(wt_data, mut_data, comparison)
    % CREATE_CONDITION_COMPARISON_PLOTS - Generate WT vs R213W comparison visualizations
    %
    % Creates comprehensive condition comparison plots:
    %   1. Event frequency distributions and box plots
    %   2. Event amplitude distributions and box plots  
    %   3. ROI activity comparison across files
    %   4. Summary statistics overview
    
    plot_handles = struct();
    
    fprintf('  Creating condition comparison plots...\n');
    
    %% Plot 1: Event Frequency Comparison
    plot_handles.frequency_comparison = create_frequency_comparison_plot(wt_data, mut_data, comparison);
    
    %% Plot 2: Event Amplitude Comparison  
    plot_handles.amplitude_comparison = create_amplitude_comparison_plot(wt_data, mut_data, comparison);
    
    %% Plot 3: ROI Activity Comparison
    plot_handles.activity_comparison = create_activity_comparison_plot(wt_data, mut_data, comparison);
    
    %% Plot 4: File-Level Summary
    plot_handles.file_summary = create_file_summary_plot(wt_data, mut_data, comparison);
    
    fprintf('    Created %d condition comparison plots\n', length(fieldnames(plot_handles)));
end

function fig = create_frequency_comparison_plot(wt_data, mut_data, comparison)
    % Compare event frequencies between WT and R213W (active ROIs only)
    
    fig = figure('Name', 'Event Frequency Comparison: WT vs R213W', ...
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
    
    % Subplot 1: Histograms overlay
    subplot(2, 2, 1);
    edges = 0:0.01:max([wt_freqs, mut_freqs]) + 0.01;
    
    hold on;
    histogram(wt_freqs, edges, 'FaceColor', [0.2, 0.6, 1.0], 'FaceAlpha', 0.7, ...
        'EdgeColor', 'none', 'DisplayName', sprintf('WT (n=%d)', length(wt_freqs)));
    histogram(mut_freqs, edges, 'FaceColor', [1.0, 0.4, 0.2], 'FaceAlpha', 0.7, ...
        'EdgeColor', 'none', 'DisplayName', sprintf('R213W (n=%d)', length(mut_freqs)));
    
    xlabel('Event Frequency (Hz)');
    ylabel('ROI Count');
    title('Event Frequency Distributions');
    legend('Location', 'best');
    grid on;
    
    % Add mean lines
    wt_mean = mean(wt_freqs);
    mut_mean = mean(mut_freqs);
    y_lim = ylim;
    plot([wt_mean, wt_mean], y_lim, 'b--', 'LineWidth', 2);
    plot([mut_mean, mut_mean], y_lim, 'r--', 'LineWidth', 2);
    
    % Subplot 2: Box plot comparison
    subplot(2, 2, 2);
    
    % Prepare data for box plot
    all_freqs = [wt_freqs, mut_freqs];
    groups = [ones(1, length(wt_freqs)), 2*ones(1, length(mut_freqs))];
    
    boxplot(all_freqs, groups, 'Labels', {'WT', 'R213W'}, 'Colors', [0.2, 0.6, 1.0; 1.0, 0.4, 0.2]);
    ylabel('Event Frequency (Hz)');
    title('Frequency Box Plot Comparison');
    grid on;
    
    % Add statistical annotation
    if isfield(comparison, 'frequency_pvalue') && ~isnan(comparison.frequency_pvalue)
        y_max = max(all_freqs) * 1.1;
        if comparison.frequency_pvalue < 0.001
            p_text = 'p < 0.001';
        else
            p_text = sprintf('p = %.3f', comparison.frequency_pvalue);
        end
        text(1.5, y_max, p_text, 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
        
        % Add significance bracket
        line([1, 2], [y_max*0.95, y_max*0.95], 'Color', 'k', 'LineWidth', 1);
        line([1, 1], [y_max*0.95, y_max*0.92], 'Color', 'k', 'LineWidth', 1);
        line([2, 2], [y_max*0.95, y_max*0.92], 'Color', 'k', 'LineWidth', 1);
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
    
    % Subplot 4: Summary statistics
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
            summary_text{end+1} = sprintf('p = %.3f **', comparison.frequency_pvalue);
        elseif comparison.frequency_pvalue < 0.05
            summary_text{end+1} = sprintf('p = %.3f *', comparison.frequency_pvalue);
        else
            summary_text{end+1} = sprintf('p = %.3f (ns)', comparison.frequency_pvalue);
        end
        
        if isfield(comparison, 'frequency_effect_size')
            summary_text{end+1} = sprintf('Effect size (Cohen''s d): %.3f', comparison.frequency_effect_size);
        end
    else
        summary_text{end+1} = 'p = N/A';
    end
    
    text(0.05, 0.95, summary_text, 'FontSize', 10, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', 'Units', 'normalized', 'FontName', 'FixedWidth');
    
    sgtitle('Event Frequency Comparison: WT vs R213W', 'FontSize', 16, 'FontWeight', 'bold');
end

function fig = create_amplitude_comparison_plot(wt_data, mut_data, comparison)
    % Compare event amplitudes between WT and R213W (all individual events)
    
    fig = figure('Name', 'Event Amplitude Comparison: WT vs R213W', ...
        'Position', [200, 100, 1200, 800]);
    
    % Extract individual event amplitudes
    wt_amps = wt_data.all_individual_events;
    mut_amps = mut_data.all_individual_events;
    
    if isempty(wt_amps) || isempty(mut_amps)
        text(0.5, 0.5, 'Insufficient amplitude data for comparison', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 14);
        title('Event Amplitude Comparison - INSUFFICIENT DATA');
        return;
    end
    
    % 2x2 layout for amplitude analysis
    
    % Subplot 1: Histograms overlay
    subplot(2, 2, 1);
    
    % Use logarithmic binning for better visualization of amplitude distributions
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
    title('Event Amplitude Distributions');
    legend('Location', 'best');
    grid on;
    
    % Add mean lines
    wt_mean = mean(wt_amps);
    mut_mean = mean(mut_amps);
    y_lim = ylim;
    plot([wt_mean, wt_mean], y_lim, 'b--', 'LineWidth', 2);
    plot([mut_mean, mut_mean], y_lim, 'r--', 'LineWidth', 2);
    
    % Subplot 2: Box plot comparison
    subplot(2, 2, 2);
    
    % Sample data if too large for visualization
    max_samples = 10000;  % Limit for box plot display
    if length(wt_amps) > max_samples
        wt_amps_plot = datasample(wt_amps, max_samples, 'Replace', false);
    else
        wt_amps_plot = wt_amps;
    end
    
    if length(mut_amps) > max_samples
        mut_amps_plot = datasample(mut_amps, max_samples, 'Replace', false);
    else
        mut_amps_plot = mut_amps;
    end
    
    % Prepare data for box plot
    all_amps = [wt_amps_plot, mut_amps_plot];
    groups = [ones(1, length(wt_amps_plot)), 2*ones(1, length(mut_amps_plot))];
    
    boxplot(all_amps, groups, 'Labels', {'WT', 'R213W'}, 'Colors', [0.2, 0.6, 1.0; 1.0, 0.4, 0.2]);
    ylabel('Event Amplitude (dF/F)');
    title('Amplitude Box Plot Comparison');
    grid on;
    
    % Add statistical annotation
    if isfield(comparison, 'amplitude_pvalue') && ~isnan(comparison.amplitude_pvalue)
        y_max = max(all_amps) * 1.1;
        if comparison.amplitude_pvalue < 0.001
            p_text = 'p < 0.001';
        else
            p_text = sprintf('p = %.3f', comparison.amplitude_pvalue);
        end
        text(1.5, y_max, p_text, 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
        
        % Add significance bracket
        line([1, 2], [y_max*0.95, y_max*0.95], 'Color', 'k', 'LineWidth', 1);
        line([1, 1], [y_max*0.95, y_max*0.92], 'Color', 'k', 'LineWidth', 1);
        line([2, 2], [y_max*0.95, y_max*0.92], 'Color', 'k', 'LineWidth', 1);
    end
    
    % Subplot 3: Log-scale histograms for wide range data
    subplot(2, 2, 3);
    
    % Create log-scale bins if amplitude range is large
    if (max_amp / min_amp) > 100
        log_edges = logspace(log10(min_amp), log10(max_amp), 30);
        hold on;
        histogram(wt_amps, log_edges, 'FaceColor', [0.2, 0.6, 1.0], 'FaceAlpha', 0.7, ...
            'EdgeColor', 'none', 'DisplayName', 'WT');
        histogram(mut_amps, log_edges, 'FaceColor', [1.0, 0.4, 0.2], 'FaceAlpha', 0.7, ...
            'EdgeColor', 'none', 'DisplayName', 'R213W');
        set(gca, 'XScale', 'log');
        xlabel('Event Amplitude (dF/F) - Log Scale');
        ylabel('Event Count');
        title('Amplitude Distributions (Log Scale)');
        legend('Location', 'best');
        grid on;
    else
        % Regular cumulative distributions
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
    end
    
    % Subplot 4: Summary statistics
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
        sprintf('Mann-Whitney U test');
    };
    
    if isfield(comparison, 'amplitude_pvalue') && ~isnan(comparison.amplitude_pvalue)
        if comparison.amplitude_pvalue < 0.001
            summary_text{end+1} = 'p < 0.001 ***';
        elseif comparison.amplitude_pvalue < 0.01
            summary_text{end+1} = sprintf('p = %.3f **', comparison.amplitude_pvalue);
        elseif comparison.amplitude_pvalue < 0.05
            summary_text{end+1} = sprintf('p = %.3f *', comparison.amplitude_pvalue);
        else
            summary_text{end+1} = sprintf('p = %.3f (ns)', comparison.amplitude_pvalue);
        end
        
        if isfield(comparison, 'amplitude_effect_size')
            summary_text{end+1} = sprintf('Effect size (Cohen''s d): %.3f', comparison.amplitude_effect_size);
        end
    else
        summary_text{end+1} = 'p = N/A';
    end
    
    text(0.05, 0.95, summary_text, 'FontSize', 10, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', 'Units', 'normalized', 'FontName', 'FixedWidth');
    
    sgtitle('Event Amplitude Comparison: WT vs R213W', 'FontSize', 16, 'FontWeight', 'bold');
end

function fig = create_activity_comparison_plot(wt_data, mut_data, comparison)
    % Compare ROI activity patterns between WT and R213W
    
    fig = figure('Name', 'ROI Activity Comparison: WT vs R213W', ...
        'Position', [300, 100, 1200, 800]);
    
    % 2x2 layout for activity analysis
    
    % Subplot 1: Events per ROI distribution
    subplot(2, 2, 1);
    
    wt_events = wt_data.all_roi_event_counts;
    mut_events = mut_data.all_roi_event_counts;
    
    if ~isempty(wt_events) && ~isempty(mut_events)
        max_events = max([wt_events, mut_events]);
        edges = -0.5:1:(max_events + 0.5);
        
        hold on;
        histogram(wt_events, edges, 'FaceColor', [0.2, 0.6, 1.0], 'FaceAlpha', 0.7, ...
            'EdgeColor', 'none', 'DisplayName', sprintf('WT (n=%d)', length(wt_events)));
        histogram(mut_events, edges, 'FaceColor', [1.0, 0.4, 0.2], 'FaceAlpha', 0.7, ...
            'EdgeColor', 'none', 'DisplayName', sprintf('R213W (n=%d)', length(mut_events)));
        
        xlabel('Events per ROI');
        ylabel('ROI Count');
        title('Events per ROI Distribution');
        legend('Location', 'best');
        grid on;
    end
    
    % Subplot 2: File-level activity comparison
    subplot(2, 2, 2);
    
    if ~isempty(wt_data.file_summaries) && ~isempty(mut_data.file_summaries)
        % Extract file-level metrics
        wt_activity_fractions = arrayfun(@(x) x.active_rois / x.num_rois, wt_data.file_summaries) * 100;
        mut_activity_fractions = arrayfun(@(x) x.active_rois / x.num_rois, mut_data.file_summaries) * 100;
        
        % Bar plot of activity fractions by file
        file_indices = 1:max(length(wt_activity_fractions), length(mut_activity_fractions));
        
        hold on;
        if ~isempty(wt_activity_fractions)
            bar(1:length(wt_activity_fractions), wt_activity_fractions, 0.4, ...
                'FaceColor', [0.2, 0.6, 1.0], 'DisplayName', 'WT');
        end
        if ~isempty(mut_activity_fractions)
            bar((1:length(mut_activity_fractions)) + 0.4, mut_activity_fractions, 0.4, ...
                'FaceColor', [1.0, 0.4, 0.2], 'DisplayName', 'R213W');
        end
        
        xlabel('File Index');
        ylabel('Active ROI Fraction (%)');
        title('File-Level Activity Comparison');
        legend('Location', 'best');
        grid on;
    end
    
    % Subplot 3: Activity category comparison
    subplot(2, 2, 3);
    
    % Create activity categories for both conditions
    if ~isempty(wt_events) && ~isempty(mut_events)
        % Define activity bins
        activity_bins = {
            '0 events', @(x) x == 0;
            '1 event', @(x) x == 1;
            '2-3 events', @(x) x >= 2 & x <= 3;
            '4-5 events', @(x) x >= 4 & x <= 5;
            '6+ events', @(x) x >= 6
        };
        
        wt_counts = zeros(1, length(activity_bins));
        mut_counts = zeros(1, length(activity_bins));
        
        for i = 1:length(activity_bins)
            wt_counts(i) = sum(activity_bins{i, 2}(wt_events));
            mut_counts(i) = sum(activity_bins{i, 2}(mut_events));
        end
        
        % Create grouped bar chart
        bar_data = [wt_counts; mut_counts]';
        bar(bar_data, 'grouped');
        
        xlabel('Activity Category');
        ylabel('ROI Count');
        title('ROI Activity Categories');
        set(gca, 'XTickLabel', activity_bins(:, 1));
        legend('WT', 'R213W', 'Location', 'best');
        grid on;
        
        % Rotate x-axis labels for better readability
        xtickangle(45);
    end
    
    % Subplot 4: Summary statistics
    subplot(2, 2, 4); axis off;
    
    summary_text = {
        'ROI ACTIVITY COMPARISON';
        '';
        sprintf('WT Summary:');
        sprintf('  Files: %d', wt_data.num_files);
        sprintf('  Total ROIs: %d', wt_data.total_rois);
        sprintf('  Active ROIs: %d (%.1f%%)', ...
            wt_data.stats.frequency.active_rois.count, ...
            100 * wt_data.stats.frequency.active_rois.fraction);
        sprintf('  Events/ROI: %.2f ± %.2f', ...
            wt_data.stats.events_per_roi.mean, wt_data.stats.events_per_roi.std);
        '';
        sprintf('R213W Summary:');
        sprintf('  Files: %d', mut_data.num_files);
        sprintf('  Total ROIs: %d', mut_data.total_rois);
        sprintf('  Active ROIs: %d (%.1f%%)', ...
            mut_data.stats.frequency.active_rois.count, ...
            100 * mut_data.stats.frequency.active_rois.fraction);
        sprintf('  Events/ROI: %.2f ± %.2f', ...
            mut_data.stats.events_per_roi.mean, mut_data.stats.events_per_roi.std);
    };
    
    % Add statistical comparison if available
    if isfield(comparison, 'events_per_roi_pvalue') && ~isnan(comparison.events_per_roi_pvalue)
        summary_text{end+1} = '';
        summary_text{end+1} = 'STATISTICAL COMPARISON:';
        summary_text{end+1} = 'Events per ROI (Mann-Whitney U):';
        
        if comparison.events_per_roi_pvalue < 0.001
            summary_text{end+1} = 'p < 0.001 ***';
        elseif comparison.events_per_roi_pvalue < 0.01
            summary_text{end+1} = sprintf('p = %.3f **', comparison.events_per_roi_pvalue);
        elseif comparison.events_per_roi_pvalue < 0.05
            summary_text{end+1} = sprintf('p = %.3f *', comparison.events_per_roi_pvalue);
        else
            summary_text{end+1} = sprintf('p = %.3f (ns)', comparison.events_per_roi_pvalue);
        end
    end
    
    text(0.05, 0.95, summary_text, 'FontSize', 10, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', 'Units', 'normalized', 'FontName', 'FixedWidth');
    
    sgtitle('Active ROI Comparison: WT vs R213W (Inactive ROIs Excluded)', 'FontSize', 16, 'FontWeight', 'bold');
end

function fig = create_file_summary_plot(wt_data, mut_data, comparison)
    % Create file-level summary comparison plot
    
    fig = figure('Name', 'File-Level Summary: WT vs R213W', ...
        'Position', [400, 100, 1200, 800]);
    
    % 2x2 layout for file summary
    
    % Subplot 1: Events per file
    subplot(2, 2, 1);
    
    if ~isempty(wt_data.file_summaries) && ~isempty(mut_data.file_summaries)
        wt_events_per_file = arrayfun(@(x) x.num_events, wt_data.file_summaries);
        mut_events_per_file = arrayfun(@(x) x.num_events, mut_data.file_summaries);
        
        file_nums = 1:max(length(wt_events_per_file), length(mut_events_per_file));
        
        hold on;
        if ~isempty(wt_events_per_file)
            bar(1:length(wt_events_per_file), wt_events_per_file, 0.4, ...
                'FaceColor', [0.2, 0.6, 1.0], 'DisplayName', 'WT');
        end
        if ~isempty(mut_events_per_file)
            bar((1:length(mut_events_per_file)) + 0.4, mut_events_per_file, 0.4, ...
                'FaceColor', [1.0, 0.4, 0.2], 'DisplayName', 'R213W');
        end
        
        xlabel('File Index');
        ylabel('Total Events');
        title('Total Events per File');
        legend('Location', 'best');
        grid on;
    end
    
    % Subplot 2: ROIs per file
    subplot(2, 2, 2);
    
    if ~isempty(wt_data.file_summaries) && ~isempty(mut_data.file_summaries)
        wt_rois_per_file = arrayfun(@(x) x.num_rois, wt_data.file_summaries);
        mut_rois_per_file = arrayfun(@(x) x.num_rois, mut_data.file_summaries);
        
        hold on;
        if ~isempty(wt_rois_per_file)
            bar(1:length(wt_rois_per_file), wt_rois_per_file, 0.4, ...
                'FaceColor', [0.2, 0.6, 1.0], 'DisplayName', 'WT');
        end
        if ~isempty(mut_rois_per_file)
            bar((1:length(mut_rois_per_file)) + 0.4, mut_rois_per_file, 0.4, ...
                'FaceColor', [1.0, 0.4, 0.2], 'DisplayName', 'R213W');
        end
        
        xlabel('File Index');
        ylabel('Total ROIs');
        title('ROIs per File');
        legend('Location', 'best');
        grid on;
    end
    
    % Subplot 3: Mean frequency per file
    subplot(2, 2, 3);
    
    if ~isempty(wt_data.file_summaries) && ~isempty(mut_data.file_summaries)
        wt_mean_freq = arrayfun(@(x) x.mean_frequency, wt_data.file_summaries);
        mut_mean_freq = arrayfun(@(x) x.mean_frequency, mut_data.file_summaries);
        
        % Handle NaN values
        wt_mean_freq(isnan(wt_mean_freq)) = 0;
        mut_mean_freq(isnan(mut_mean_freq)) = 0;
        
        hold on;
        if ~isempty(wt_mean_freq)
            bar(1:length(wt_mean_freq), wt_mean_freq, 0.4, ...
                'FaceColor', [0.2, 0.6, 1.0], 'DisplayName', 'WT');
        end
        if ~isempty(mut_mean_freq)
            bar((1:length(mut_mean_freq)) + 0.4, mut_mean_freq, 0.4, ...
                'FaceColor', [1.0, 0.4, 0.2], 'DisplayName', 'R213W');
        end
        
        xlabel('File Index');
        ylabel('Mean Frequency (Hz)');
        title('Mean Event Frequency per File');
        legend('Location', 'best');
        grid on;
    end
    
    % Subplot 4: Overall comparison summary
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
    };
    
    % Add significance summary
    if isfield(comparison, 'summary')
        summary_text{end+1} = '';
        summary_text{end+1} = sprintf('Statistical Tests:');
        summary_text{end+1} = sprintf('  %d/%d tests significant (p < 0.05)', ...
            comparison.summary.significant_tests, comparison.summary.total_tests);
    end
    
    text(0.05, 0.95, summary_text, 'FontSize', 10, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', 'Units', 'normalized', 'FontName', 'FixedWidth');
    
    sgtitle('File-Level Summary: WT vs R213W', 'FontSize', 16, 'FontWeight', 'bold');
end