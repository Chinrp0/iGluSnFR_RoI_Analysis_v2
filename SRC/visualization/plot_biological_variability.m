function fig = plot_biological_variability(results, options)
    % PLOT_BIOLOGICAL_VARIABILITY - Shows variability across biological replicates
    % Violin plots, individual data points, and coefficient of variation
    %
    % Args:
    %   results - Output from integrated_batch_analysis_final()
    %   options - Optional struct with plotting preferences
    
    if nargin < 2
        options = struct();
    end
    
    fprintf('Creating biological variability analysis...\n');
    
    %% Extract per-file metrics
    [wt_metrics, wt_filenames] = extract_file_metrics(results.file_results.wt_results, ...
        results.processing_info.wt_files);
    [mut_metrics, mut_filenames] = extract_file_metrics(results.file_results.mut_results, ...
        results.processing_info.mut_files);
    
    % FIXED: Use length of field array, not struct itself
    fprintf('  WT: %d files\n', length(wt_metrics.mean_frequency));
    fprintf('  R213W: %d files\n', length(mut_metrics.mean_frequency));
    
    %% Create figure
    fig = figure('Name', 'Biological Variability Analysis: WT vs R213W', ...
        'Position', [100, 100, 1400, 1000]);
    
    %% Subplot 1: Mean Event Frequency per File
    subplot(2, 3, 1);
    plot_metric_comparison(wt_metrics.mean_frequency, mut_metrics.mean_frequency, ...
        'WT', 'R213W', 'Mean Event Frequency (Hz)', [0.2, 0.6, 1.0], [1.0, 0.4, 0.2]);
    
    %% Subplot 2: Mean Event Amplitude per File
    subplot(2, 3, 2);
    plot_metric_comparison(wt_metrics.mean_amplitude, mut_metrics.mean_amplitude, ...
        'WT', 'R213W', 'Mean Event Amplitude (dF/F)', [0.2, 0.6, 1.0], [1.0, 0.4, 0.2]);
    
    %% Subplot 3: Active ROI Fraction per File
    subplot(2, 3, 3);
    plot_metric_comparison(wt_metrics.active_roi_fraction * 100, mut_metrics.active_roi_fraction * 100, ...
        'WT', 'R213W', 'Active ROI Fraction (%)', [0.2, 0.6, 1.0], [1.0, 0.4, 0.2]);
    
    %% Subplot 4: Total Events per File
    subplot(2, 3, 4);
    plot_metric_comparison(wt_metrics.total_events, mut_metrics.total_events, ...
        'WT', 'R213W', 'Total Events', [0.2, 0.6, 1.0], [1.0, 0.4, 0.2]);
    
    %% Subplot 5: Events per Active ROI
    subplot(2, 3, 5);
    wt_events_per_active = wt_metrics.total_events ./ wt_metrics.active_rois;
    mut_events_per_active = mut_metrics.total_events ./ mut_metrics.active_rois;
    
    plot_metric_comparison(wt_events_per_active, mut_events_per_active, ...
        'WT', 'R213W', 'Events per Active ROI', [0.2, 0.6, 1.0], [1.0, 0.4, 0.2]);
    
    %% Subplot 6: Coefficient of Variation Summary
    subplot(2, 3, 6); axis off;
    
    % Calculate CVs
    cv_wt_freq = std(wt_metrics.mean_frequency) / mean(wt_metrics.mean_frequency);
    cv_mut_freq = std(mut_metrics.mean_frequency) / mean(mut_metrics.mean_frequency);
    cv_wt_amp = std(wt_metrics.mean_amplitude) / mean(wt_metrics.mean_amplitude);
    cv_mut_amp = std(mut_metrics.mean_amplitude) / mean(mut_metrics.mean_amplitude);
    cv_wt_active = std(wt_metrics.active_roi_fraction) / mean(wt_metrics.active_roi_fraction);
    cv_mut_active = std(mut_metrics.active_roi_fraction) / mean(mut_metrics.active_roi_fraction);
    
    summary_text = {
        'BIOLOGICAL VARIABILITY SUMMARY';
        '(Coefficient of Variation)';
        '';
        'Event Frequency:';
        sprintf('  WT: %.3f', cv_wt_freq);
        sprintf('  R213W: %.3f', cv_mut_freq);
        '';
        'Event Amplitude:';
        sprintf('  WT: %.3f', cv_wt_amp);
        sprintf('  R213W: %.3f', cv_mut_amp);
        '';
        'Active ROI Fraction:';
        sprintf('  WT: %.3f', cv_wt_active);
        sprintf('  R213W: %.3f', cv_mut_active);
        '';
        'Sample Sizes:';
        sprintf('  WT: %d files', length(wt_metrics.mean_frequency));
        sprintf('  R213W: %d files', length(mut_metrics.mean_frequency));
        '';
        'NOTE:';
        'Lower CV = more consistent';
        'Higher CV = more variable';
    };
    
    text(0.05, 0.95, summary_text, 'FontSize', 10, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', 'Units', 'normalized', 'FontName', 'FixedWidth');
    
    sgtitle('Biological Variability Analysis: WT vs R213W', 'FontSize', 16, 'FontWeight', 'bold');
    
    fprintf('  Created biological variability analysis\n');
end

function plot_metric_comparison(wt_data, mut_data, wt_label, mut_label, ylabel_text, wt_color, mut_color)
    % Helper function to create consistent comparison plots
    
    % Combine data
    all_data = [wt_data; mut_data];
    n_wt = length(wt_data);
    n_mut = length(mut_data);
    groups = [ones(n_wt, 1); 2 * ones(n_mut, 1)];
    
    % Box plot
    boxplot(all_data, groups, 'Labels', {wt_label, mut_label}, ...
        'Colors', [wt_color; mut_color], 'Symbol', '');
    hold on;
    
    % Overlay individual points with jitter
    jitter_amount = 0.15;
    wt_x = ones(n_wt, 1) + (rand(n_wt, 1) - 0.5) * jitter_amount;
    mut_x = 2 * ones(n_mut, 1) + (rand(n_mut, 1) - 0.5) * jitter_amount;
    
    scatter(wt_x, wt_data, 50, wt_color, 'filled', 'MarkerFaceAlpha', 0.6);
    scatter(mut_x, mut_data, 50, mut_color, 'filled', 'MarkerFaceAlpha', 0.6);
    
    % Add means as diamonds
    plot(1, mean(wt_data), 'd', 'MarkerSize', 12, 'MarkerFaceColor', 'blue', ...
        'MarkerEdgeColor', 'black', 'LineWidth', 2);
    plot(2, mean(mut_data), 'd', 'MarkerSize', 12, 'MarkerFaceColor', 'red', ...
        'MarkerEdgeColor', 'black', 'LineWidth', 2);
    
    ylabel(ylabel_text);
    grid on;
    
    % Statistical test
    if n_wt > 1 && n_mut > 1
        [p_value, ~] = ranksum(wt_data, mut_data);
        y_max = max(all_data);
        y_min = min(all_data);
        y_range = y_max - y_min;
        y_text = y_max + y_range * 0.15;
        
        if p_value < 0.001
            p_text = 'p < 0.001 ***';
        elseif p_value < 0.01
            p_text = sprintf('p = %.4f **', p_value);
        elseif p_value < 0.05
            p_text = sprintf('p = %.4f *', p_value);
        else
            p_text = sprintf('p = %.4f n.s.', p_value);
        end
        
        text(1.5, y_text, p_text, 'HorizontalAlignment', 'center', ...
            'FontSize', 9, 'FontWeight', 'bold');
        
        % Significance bracket
        line([1, 2], [y_text * 0.95, y_text * 0.95], 'Color', 'k', 'LineWidth', 1);
        line([1, 1], [y_text * 0.95, y_text * 0.90], 'Color', 'k', 'LineWidth', 1);
        line([2, 2], [y_text * 0.95, y_text * 0.90], 'Color', 'k', 'LineWidth', 1);
        
        ylim([y_min - y_range * 0.05, y_text + y_range * 0.05]);
    end
end

function [metrics, filenames] = extract_file_metrics(file_results, file_paths)
    % Extract per-file summary metrics
    
    n_files = length(file_results);
    
    metrics = struct();
    metrics.mean_frequency = zeros(n_files, 1);
    metrics.mean_amplitude = zeros(n_files, 1);
    metrics.active_roi_fraction = zeros(n_files, 1);
    metrics.total_events = zeros(n_files, 1);
    metrics.active_rois = zeros(n_files, 1);
    metrics.total_rois = zeros(n_files, 1);
    
    filenames = cell(n_files, 1);
    
    for i = 1:n_files
        if isempty(file_results{i})
            continue;
        end
        
        file_result = file_results{i};
        
        % Extract filename
        [~, fname, ~] = fileparts(file_paths{i});
        filenames{i} = fname;
        
        % Event statistics
        if isfield(file_result, 'event_stats')
            es = file_result.event_stats;
            
            % Mean frequency (only active ROIs)
            if isfield(es, 'frequency_hz')
                active_freqs = es.frequency_hz(es.frequency_hz > 0);
                if ~isempty(active_freqs)
                    metrics.mean_frequency(i) = mean(active_freqs);
                end
            end
            
            % Mean amplitude
            if isfield(es, 'all_event_peaks') && ~isempty(es.all_event_peaks)
                metrics.mean_amplitude(i) = mean(es.all_event_peaks);
            end
            
            % Total events
            if isfield(es, 'total_events')
                metrics.total_events(i) = es.total_events;
            end
            
            % Active ROIs
            if isfield(es, 'rois_with_events')
                metrics.active_rois(i) = es.rois_with_events;
            end
        end
        
        % Total ROIs
        if isfield(file_result, 'numROIs')
            metrics.total_rois(i) = file_result.numROIs;
            metrics.active_roi_fraction(i) = metrics.active_rois(i) / metrics.total_rois(i);
        end
    end
end
