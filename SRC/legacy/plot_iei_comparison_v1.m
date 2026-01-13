function fig = plot_iei_comparison(results, options)
    % PLOT_IEI_COMPARISON - Inter-Event Interval (IEI) distribution comparison
    % Shows temporal patterning of events: regular vs irregular firing
    %
    % Args:
    %   results - Output from integrated_batch_analysis_final()
    %   options - Optional struct with:
    %             .frame_rate (default: 100 Hz for 10ms exposure)
    %             .max_iei_s (default: 10s, max IEI to plot)
    %             .bin_size_ms (default: 50ms)
    
    if nargin < 2
        options = struct();
    end
    
    % Set defaults
    if ~isfield(options, 'frame_rate'), options.frame_rate = 100; end  % 10ms exposure = 100 Hz
    if ~isfield(options, 'max_iei_s'), options.max_iei_s = 10; end
    if ~isfield(options, 'bin_size_ms'), options.bin_size_ms = 50; end
    
    fprintf('Creating inter-event interval comparison...\n');
    
    %% Extract IEIs from both conditions
    [wt_ieis, wt_stats] = extract_all_ieis(results.file_results.wt_results, options.frame_rate);
    [mut_ieis, mut_stats] = extract_all_ieis(results.file_results.mut_results, options.frame_rate);
    
    fprintf('  WT: %d IEIs from %d ROIs\n', length(wt_ieis), wt_stats.num_rois_with_multiple_events);
    fprintf('  R213W: %d IEIs from %d ROIs\n', length(mut_ieis), mut_stats.num_rois_with_multiple_events);
    
    %% Create figure
    fig = figure('Name', 'Inter-Event Interval Comparison: WT vs R213W', ...
        'Position', [100, 100, 1400, 900]);
    
    % Filter to reasonable range
    wt_ieis_plot = wt_ieis(wt_ieis <= options.max_iei_s);
    mut_ieis_plot = mut_ieis(mut_ieis <= options.max_iei_s);
    
    %% Subplot 1: IEI Distribution (log scale)
    subplot(2, 3, 1);
    edges = logspace(log10(0.01), log10(options.max_iei_s), 50);
    
    hold on;
    histogram(wt_ieis_plot, edges, 'FaceColor', [0.2, 0.6, 1.0], 'FaceAlpha', 0.7, ...
        'EdgeColor', 'none', 'DisplayName', sprintf('WT (n=%d)', length(wt_ieis_plot)));
    histogram(mut_ieis_plot, edges, 'FaceColor', [1.0, 0.4, 0.2], 'FaceAlpha', 0.7, ...
        'EdgeColor', 'none', 'DisplayName', sprintf('R213W (n=%d)', length(mut_ieis_plot)));
    
    set(gca, 'XScale', 'log');
    xlabel('Inter-Event Interval (s)');
    ylabel('Count');
    title('IEI Distribution (Log Scale)');
    legend('Location', 'best');
    grid on;
    
    %% Subplot 2: IEI Distribution (linear scale, zoomed)
    subplot(2, 3, 2);
    edges_linear = 0:options.bin_size_ms/1000:5;  % 0-5s with specified bin size
    
    hold on;
    histogram(wt_ieis_plot, edges_linear, 'FaceColor', [0.2, 0.6, 1.0], 'FaceAlpha', 0.7, ...
        'EdgeColor', 'none', 'DisplayName', 'WT');
    histogram(mut_ieis_plot, edges_linear, 'FaceColor', [1.0, 0.4, 0.2], 'FaceAlpha', 0.7, ...
        'EdgeColor', 'none', 'DisplayName', 'R213W');
    
    xlabel('Inter-Event Interval (s)');
    ylabel('Count');
    title('IEI Distribution (Linear Scale, 0-5s)');
    legend('Location', 'best');
    grid on;
    xlim([0, 5]);
    
    %% Subplot 3: Cumulative Distribution
    subplot(2, 3, 3);
    [wt_f, wt_x] = ecdf(wt_ieis_plot);
    [mut_f, mut_x] = ecdf(mut_ieis_plot);
    
    plot(wt_x, wt_f, 'b-', 'LineWidth', 2, 'DisplayName', 'WT');
    hold on;
    plot(mut_x, mut_f, 'r-', 'LineWidth', 2, 'DisplayName', 'R213W');
    
    xlabel('Inter-Event Interval (s)');
    ylabel('Cumulative Probability');
    title('Cumulative IEI Distribution');
    legend('Location', 'best');
    grid on;
    
    % Add median lines
    wt_median = median(wt_ieis_plot);
    mut_median = median(mut_ieis_plot);
    xline(wt_median, 'b--', sprintf('WT: %.2fs', wt_median), 'LineWidth', 1.5);
    xline(mut_median, 'r--', sprintf('R213W: %.2fs', mut_median), 'LineWidth', 1.5);
    
    %% Subplot 4: Box Plot Comparison
    subplot(2, 3, 4);
    all_ieis = [wt_ieis_plot, mut_ieis_plot];
    groups = [ones(1, length(wt_ieis_plot)), 2*ones(1, length(mut_ieis_plot))];
    
    boxplot(all_ieis, groups, 'Labels', {'WT', 'R213W'}, ...
        'Colors', [0.2, 0.6, 1.0; 1.0, 0.4, 0.2]);
    ylabel('Inter-Event Interval (s)');
    title('IEI Box Plot Comparison');
    grid on;
    
    % Add statistical test
    if length(wt_ieis_plot) > 0 && length(mut_ieis_plot) > 0
        [p_value, ~] = ranksum(wt_ieis_plot, mut_ieis_plot);
        y_max = max(all_ieis);
        y_range = y_max - min(all_ieis);
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
            'FontSize', 11, 'FontWeight', 'bold');
        
        % Significance bracket
        line([1, 2], [y_text*0.95, y_text*0.95], 'Color', 'k', 'LineWidth', 1);
        line([1, 1], [y_text*0.95, y_text*0.9], 'Color', 'k', 'LineWidth', 1);
        line([2, 2], [y_text*0.95, y_text*0.9], 'Color', 'k', 'LineWidth', 1);
        
        ylim([min(all_ieis) - y_range*0.05, y_text + y_range*0.05]);
    end
    
    %% Subplot 5: Coefficient of Variation
    subplot(2, 3, 5);
    
    wt_cv = wt_stats.cv_per_roi;
    mut_cv = mut_stats.cv_per_roi;
    
    wt_cv_valid = wt_cv(~isnan(wt_cv) & ~isinf(wt_cv));
    mut_cv_valid = mut_cv(~isnan(mut_cv) & ~isinf(mut_cv));
    
    if ~isempty(wt_cv_valid) && ~isempty(mut_cv_valid)
        edges_cv = 0:0.1:3;
        hold on;
        histogram(wt_cv_valid, edges_cv, 'FaceColor', [0.2, 0.6, 1.0], 'FaceAlpha', 0.7, ...
            'EdgeColor', 'none', 'DisplayName', 'WT', 'Normalization', 'probability');
        histogram(mut_cv_valid, edges_cv, 'FaceColor', [1.0, 0.4, 0.2], 'FaceAlpha', 0.7, ...
            'EdgeColor', 'none', 'DisplayName', 'R213W', 'Normalization', 'probability');
        
        xlabel('Coefficient of Variation (IEI)');
        ylabel('Probability');
        title('IEI Regularity (CV)');
        legend('Location', 'best');
        grid on;
        
        % Add reference line at CV=1 (Poisson process)
        xline(1, 'k--', 'Random (Poisson)', 'LineWidth', 1.5, 'Alpha', 0.7);
    else
        text(0.5, 0.5, 'Insufficient data for CV analysis', ...
            'HorizontalAlignment', 'center', 'FontSize', 12);
    end
    
    %% Subplot 6: Summary Statistics
    subplot(2, 3, 6); axis off;
    
    summary_text = {
        'INTER-EVENT INTERVAL SUMMARY';
        '';
        sprintf('WT (n = %d IEIs, %d ROIs):', length(wt_ieis), wt_stats.num_rois_with_multiple_events);
        sprintf('  Mean IEI: %.3f ± %.3f s', mean(wt_ieis_plot), std(wt_ieis_plot));
        sprintf('  Median IEI: %.3f s', median(wt_ieis_plot));
        sprintf('  Mean CV: %.3f ± %.3f', mean(wt_cv_valid), std(wt_cv_valid));
        sprintf('  Range: %.3f - %.3f s', min(wt_ieis_plot), max(wt_ieis_plot));
        '';
        sprintf('R213W (n = %d IEIs, %d ROIs):', length(mut_ieis), mut_stats.num_rois_with_multiple_events);
        sprintf('  Mean IEI: %.3f ± %.3f s', mean(mut_ieis_plot), std(mut_ieis_plot));
        sprintf('  Median IEI: %.3f s', median(mut_ieis_plot));
        sprintf('  Mean CV: %.3f ± %.3f', mean(mut_cv_valid), std(mut_cv_valid));
        sprintf('  Range: %.3f - %.3f s', min(mut_ieis_plot), max(mut_ieis_plot));
        '';
        'STATISTICAL TEST:';
        'Mann-Whitney U test';
    };
    
    if length(wt_ieis_plot) > 0 && length(mut_ieis_plot) > 0
        [p_value, ~] = ranksum(wt_ieis_plot, mut_ieis_plot);
        if p_value < 0.001
            summary_text{end+1} = 'p < 0.001 ***';
        else
            summary_text{end+1} = sprintf('p = %.4f', p_value);
        end
    end
    
    text(0.05, 0.95, summary_text, 'FontSize', 10, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', 'Units', 'normalized', 'FontName', 'FixedWidth');
    
    sgtitle('Inter-Event Interval Comparison: WT vs R213W', 'FontSize', 16, 'FontWeight', 'bold');
    
    fprintf('  Created IEI comparison plot\n');
end

function [all_ieis, stats] = extract_all_ieis(file_results, frame_rate)
    % Extract all inter-event intervals from file results
    
    all_ieis = [];
    cv_per_roi = [];
    num_rois_with_multiple_events = 0;
    
    for file_idx = 1:length(file_results)
        if isempty(file_results{file_idx})
            continue;
        end
        
        file_result = file_results{file_idx};
        
        if ~isfield(file_result, 'event_mask') || isempty(file_result.event_mask)
            continue;
        end
        
        event_mask = file_result.event_mask;
        [num_frames, num_rois] = size(event_mask);
        
        % Process each ROI
        for roi = 1:num_rois
            roi_events = event_mask(:, roi);
            
            % Find event start times (transitions from false to true)
            event_starts = find(diff([false; roi_events]) == 1);
            
            if length(event_starts) >= 2
                % Calculate IEIs in seconds
                ieis_frames = diff(event_starts);
                ieis_seconds = ieis_frames / frame_rate;
                
                all_ieis = [all_ieis; ieis_seconds];
                
                % Calculate CV for this ROI
                if length(ieis_seconds) >= 2
                    roi_cv = std(ieis_seconds) / mean(ieis_seconds);
                    cv_per_roi = [cv_per_roi; roi_cv];
                    num_rois_with_multiple_events = num_rois_with_multiple_events + 1;
                end
            end
        end
    end
    
    stats = struct();
    stats.num_rois_with_multiple_events = num_rois_with_multiple_events;
    stats.cv_per_roi = cv_per_roi;
end
