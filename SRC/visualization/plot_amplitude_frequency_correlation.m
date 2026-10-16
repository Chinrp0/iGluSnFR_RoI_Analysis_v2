function fig = plot_amplitude_frequency_correlation(results, options)
    [lblA, lblB, colA, colB] = spont_group_labels(results);   % group labels/colours
    % PLOT_AMPLITUDE_FREQUENCY_CORRELATION - Scatter plot of amplitude vs frequency
    % Tests if high-frequency ROIs have different amplitude characteristics
    %
    % Args:
    %   results - Output from spont_batch_analysis()
    %   options - Unused (kept for call compatibility); frequencies and
    %             amplitudes come precomputed in the file results.
    
    if nargin < 2
        options = struct();
    end
    
    fprintf('Creating amplitude-frequency correlation analysis...\n');
    
    %% Extract per-ROI metrics
    [wt_freq, wt_amp, wt_stats] = extract_roi_metrics(results.file_results.wt_results, options);
    [mut_freq, mut_amp, mut_stats] = extract_roi_metrics(results.file_results.mut_results, options);
    
    fprintf(['  ' lblA ': %d active ROIs\n'], length(wt_freq));
    fprintf(['  ' lblB ': %d active ROIs\n'], length(mut_freq));
    
    %% Create figure
    fig = figure('Name', ['Amplitude-Frequency Correlation: ' lblA ' vs ' lblB], ...
        'Position', [100, 100, 1400, 1000]);
    
    %% Subplot 1: WT Scatter with marginal histograms
    % Main scatter plot
    h1 = subplot(3, 3, [4, 7]);
    scatter(wt_freq, wt_amp, 20, colA, 'filled', 'MarkerFaceAlpha', 0.5);
    xlabel('Event Frequency (Hz)');
    ylabel('Mean Event Amplitude (dF/F)');
    title([lblA ': Amplitude vs Frequency']);
    grid on;
    
    % Add correlation info
    if length(wt_freq) > 3
        [r_wt, p_wt] = corr(wt_freq', wt_amp', 'Type', 'Spearman');
        text(0.05, 0.95, sprintf('r = %.3f\np = %.4f', r_wt, p_wt), ...
            'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize', 10);
    end
    
    % Marginal histogram - frequency (top)
    h2 = subplot(3, 3, 1);
    histogram(wt_freq, 30, 'FaceColor', colA, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    title([lblA ' Frequency Distribution']);
    ylabel('Count');
    set(gca, 'XTickLabel', []);
    grid on;
    xlim(h1.XLim);
    
    % Marginal histogram - amplitude (right)
    h3 = subplot(3, 3, 8);
    histogram(wt_amp, 30, 'FaceColor', colA, 'EdgeColor', 'none', ...
        'FaceAlpha', 0.7, 'Orientation', 'horizontal');
    xlabel('Count');
    set(gca, 'YTickLabel', []);
    grid on;
    ylim(h1.YLim);
    
    %% Subplot 2: R213W Scatter with marginal histograms
    % Main scatter plot
    h4 = subplot(3, 3, [5, 8]);
    scatter(mut_freq, mut_amp, 20, colB, 'filled', 'MarkerFaceAlpha', 0.5);
    xlabel('Event Frequency (Hz)');
    ylabel('Mean Event Amplitude (dF/F)');
    title([lblB ': Amplitude vs Frequency']);
    grid on;
    
    % Add correlation info
    if length(mut_freq) > 3
        [r_mut, p_mut] = corr(mut_freq', mut_amp', 'Type', 'Spearman');
        text(0.05, 0.95, sprintf('r = %.3f\np = %.4f', r_mut, p_mut), ...
            'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize', 10);
    end
    
    % Marginal histogram - frequency (top)
    h5 = subplot(3, 3, 2);
    histogram(mut_freq, 30, 'FaceColor', colB, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    title([lblB ' Frequency Distribution']);
    ylabel('Count');
    set(gca, 'XTickLabel', []);
    grid on;
    xlim(h4.XLim);
    
    % Marginal histogram - amplitude (right)
    h6 = subplot(3, 3, 9);
    histogram(mut_amp, 30, 'FaceColor', colB, 'EdgeColor', 'none', ...
        'FaceAlpha', 0.7, 'Orientation', 'horizontal');
    xlabel('Count');
    set(gca, 'YTickLabel', []);
    grid on;
    ylim(h4.YLim);
    
    %% Subplot 3: Combined overlay comparison
    h7 = subplot(3, 3, [3, 6, 9]);
    hold on;
    scatter(wt_freq, wt_amp, 30, colA, 'filled', 'MarkerFaceAlpha', 0.3, ...
        'DisplayName', sprintf([lblA ' (n=%d)'], length(wt_freq)));
    scatter(mut_freq, mut_amp, 30, colB, 'filled', 'MarkerFaceAlpha', 0.3, ...
        'DisplayName', sprintf([lblB ' (n=%d)'], length(mut_freq)));
    
    xlabel('Event Frequency (Hz)');
    ylabel('Mean Event Amplitude (dF/F)');
    title('Combined Comparison');
    legend('Location', 'best');
    grid on;
    
    % Add trend lines if sufficient data
    if length(wt_freq) > 10
        fit_wt = polyfit(wt_freq, wt_amp, 1);
        x_fit = linspace(min(wt_freq), max(wt_freq), 100);
        y_fit_wt = polyval(fit_wt, x_fit);
        plot(x_fit, y_fit_wt, '-', 'Color', colA, 'LineWidth', 2, 'HandleVisibility', 'off');
    end
    
    if length(mut_freq) > 10
        fit_mut = polyfit(mut_freq, mut_amp, 1);
        x_fit = linspace(min(mut_freq), max(mut_freq), 100);
        y_fit_mut = polyval(fit_mut, x_fit);
        plot(x_fit, y_fit_mut, '-', 'Color', colB, 'LineWidth', 2, 'HandleVisibility', 'off');
    end
    
    %% Summary text box
    summary_text = {
        'CORRELATION SUMMARY';
        '';
        [lblA ':'];
    };
    
    if length(wt_freq) > 3
        summary_text{end+1} = sprintf('  r = %.3f, p = %.4f', r_wt, p_wt);
    else
        summary_text{end+1} = '  Insufficient data';
    end
    
    summary_text{end+1} = '';
    summary_text{end+1} = [lblB ':'];
    
    if length(mut_freq) > 3
        summary_text{end+1} = sprintf('  r = %.3f, p = %.4f', r_mut, p_mut);
    else
        summary_text{end+1} = '  Insufficient data';
    end
    
    text(0.98, 0.02, summary_text, 'Units', 'normalized', ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
        'BackgroundColor', 'white', 'EdgeColor', 'black', ...
        'FontSize', 10, 'FontName', 'FixedWidth');
    
    sgtitle(['Amplitude-Frequency Correlation: ' lblA ' vs ' lblB], 'FontSize', 16, 'FontWeight', 'bold');
    
    fprintf('  Created amplitude-frequency correlation plot\n');
end

function [frequencies, amplitudes, stats] = extract_roi_metrics(file_results, options)
    % Extract per-ROI frequency and mean amplitude metrics
    
    frequencies = [];
    amplitudes = [];
    
    for file_idx = 1:length(file_results)
        if isempty(file_results{file_idx})
            continue;
        end
        
        file_result = file_results{file_idx};
        
        if ~isfield(file_result, 'event_stats') || isempty(file_result.event_stats)
            continue;
        end
        
        es = file_result.event_stats;
        
        % Get frequency for each ROI (only active ROIs)
        if isfield(es, 'frequency_hz')
            roi_freqs = es.frequency_hz;
            active_mask = roi_freqs > 0;
            
            % Get corresponding mean amplitudes for active ROIs
            if isfield(es, 'mean_peak_amplitude_per_roi')
                roi_amps = es.mean_peak_amplitude_per_roi;
                
                % Only include ROIs that have events
                frequencies = [frequencies, roi_freqs(active_mask)];
                amplitudes = [amplitudes, roi_amps(active_mask)];
            end
        end
    end
    
    stats = struct();
    stats.n_rois = length(frequencies);
    
    if length(frequencies) > 3
        [stats.correlation, stats.p_value] = corr(frequencies', amplitudes', 'Type', 'Spearman');
    else
        stats.correlation = NaN;
        stats.p_value = NaN;
    end
end
