function fig = plot_cumulative_events_over_time(results, options)
    % PLOT_CUMULATIVE_EVENTS_OVER_TIME - Shows cumulative event counts over recording
    % Detects if there's rundown, buildup, or constant activity
    %
    % Args:
    %   results - Output from integrated_batch_analysis_final()
    %   options - Optional struct with:
    %             .frame_rate (default: 100 Hz)
    %             .recording_duration_s (default: 30s)
    
    if nargin < 2
        options = struct();
    end
    
    % Set defaults
    if ~isfield(options, 'frame_rate'), options.frame_rate = 100; end
    if ~isfield(options, 'recording_duration_s'), options.recording_duration_s = 30; end
    
    fprintf('Creating cumulative event count analysis...\n');
    
    %% Extract cumulative data
    [wt_cum_mean, wt_cum_files, wt_stats] = extract_cumulative_data(...
        results.file_results.wt_results, options.frame_rate);
    [mut_cum_mean, mut_cum_files, mut_stats] = extract_cumulative_data(...
        results.file_results.mut_results, options.frame_rate);
    
    fprintf('  WT: %d files\n', wt_stats.n_files);
    fprintf('  R213W: %d files\n', mut_stats.n_files);
    
    %% Create figure
    fig = figure('Name', 'Cumulative Event Count: WT vs R213W', ...
        'Position', [100, 100, 1400, 900]);
    
    %% Subplot 1: WT Cumulative Curves
    subplot(2, 3, 1);
    plot_cumulative_curves(wt_cum_files, wt_cum_mean, 'WT', [0.2, 0.6, 1.0], options);
    
    %% Subplot 2: R213W Cumulative Curves
    subplot(2, 3, 2);
    plot_cumulative_curves(mut_cum_files, mut_cum_mean, 'R213W', [1.0, 0.4, 0.2], options);
    
    %% Subplot 3: Mean Comparison
    subplot(2, 3, 3);
    plot(wt_cum_mean.time, wt_cum_mean.cumulative_normalized, 'b-', 'LineWidth', 3, ...
        'DisplayName', 'WT Mean');
    hold on;
    plot(mut_cum_mean.time, mut_cum_mean.cumulative_normalized, 'r-', 'LineWidth', 3, ...
        'DisplayName', 'R213W Mean');
    
    xlabel('Time (s)');
    ylabel('Normalized Cumulative Events');
    title('Mean Cumulative Curves (Normalized to 1.0)');
    legend('Location', 'northwest');
    grid on;
    xlim([0, options.recording_duration_s]);
    ylim([0, 1.1]);
    
    % Add linear reference
    plot([0, options.recording_duration_s], [0, 1], 'k--', 'LineWidth', 1.5, ...
        'DisplayName', 'Linear (constant rate)');
    
    %% Subplot 4: Instantaneous Rate (WT)
    subplot(2, 3, 4);
    plot_instantaneous_rate(wt_cum_mean, 'WT', [0.2, 0.6, 1.0], options);
    
    %% Subplot 5: Instantaneous Rate (R213W)
    subplot(2, 3, 5);
    plot_instantaneous_rate(mut_cum_mean, 'R213W', [1.0, 0.4, 0.2], options);
    
    %% Subplot 6: Summary Statistics
    subplot(2, 3, 6); axis off;
    
    % Calculate linearity metrics
    wt_linearity = calculate_linearity(wt_cum_mean);
    mut_linearity = calculate_linearity(mut_cum_mean);
    
    summary_text = {
        'CUMULATIVE EVENT ANALYSIS';
        '';
        'WT:';
        sprintf('  Linearity: %.4f', wt_linearity.r_squared);
        sprintf('  Early rate: %.1f events/s', wt_stats.early_rate);
        sprintf('  Late rate: %.1f events/s', wt_stats.late_rate);
        sprintf('  Rate change: %+.1f%%', wt_stats.rate_change_percent);
        '';
        'R213W:';
        sprintf('  Linearity: %.4f', mut_linearity.r_squared);
        sprintf('  Early rate: %.1f events/s', mut_stats.early_rate);
        sprintf('  Late rate: %.1f events/s', mut_stats.late_rate);
        sprintf('  Rate change: %+.1f%%', mut_stats.rate_change_percent);
        '';
        'Interpretation:';
        sprintf('  R² near 1.0 = constant rate');
        sprintf('  Positive rate change = buildup');
        sprintf('  Negative rate change = rundown');
    };
    
    text(0.05, 0.95, summary_text, 'FontSize', 10, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', 'Units', 'normalized', 'FontName', 'FixedWidth');
    
    sgtitle('Cumulative Event Count Over Time: WT vs R213W', 'FontSize', 16, 'FontWeight', 'bold');
    
    fprintf('  Created cumulative event count analysis\n');
end

function plot_cumulative_curves(file_curves, mean_curve, condition_name, color, options)
    % Plot individual file curves with mean overlay
    
    % Plot individual files in lighter color
    for i = 1:length(file_curves)
        if ~isempty(file_curves{i})
            plot(file_curves{i}.time, file_curves{i}.cumulative_normalized, '-', ...
                'Color', [color, 0.2], 'LineWidth', 0.8);
            hold on;
        end
    end
    
    % Plot mean in bold
    plot(mean_curve.time, mean_curve.cumulative_normalized, '-', ...
        'Color', color, 'LineWidth', 3);
    
    % Add linear reference
    plot([0, options.recording_duration_s], [0, 1], 'k--', 'LineWidth', 1.5);
    
    xlabel('Time (s)');
    ylabel('Normalized Cumulative Events');
    title(sprintf('%s Cumulative Curves (n=%d files)', condition_name, length(file_curves)));
    grid on;
    xlim([0, options.recording_duration_s]);
    ylim([0, 1.1]);
end

function plot_instantaneous_rate(cum_data, condition_name, color, options)
    % Plot instantaneous event rate over time
    
    % Calculate rate as derivative of cumulative curve
    time = cum_data.time;
    cumulative = cum_data.cumulative;
    
    % Smooth and differentiate
    window_size = 10;  % Smooth over 10 time points
    if length(cumulative) > window_size
        cumulative_smooth = movmean(cumulative, window_size);
        rate = diff(cumulative_smooth) ./ diff(time);
        time_rate = time(1:end-1) + diff(time) / 2;
        
        % Plot as area
        area(time_rate, rate, 'FaceColor', color, 'FaceAlpha', 0.5, 'EdgeColor', color, 'LineWidth', 2);
        hold on;
        
        % Add mean line
        mean_rate = mean(rate);
        yline(mean_rate, 'k--', sprintf('Mean: %.1f events/s', mean_rate), ...
            'LineWidth', 2, 'LabelHorizontalAlignment', 'left');
        
        xlabel('Time (s)');
        ylabel('Event Rate (events/s)');
        title(sprintf('%s Instantaneous Rate', condition_name));
        grid on;
        xlim([0, options.recording_duration_s]);
    else
        text(0.5, 0.5, 'Insufficient data', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'FontSize', 12);
        title(sprintf('%s Instantaneous Rate', condition_name));
    end
end

function [cum_mean, cum_files, stats] = extract_cumulative_data(file_results, frame_rate)
    % Extract cumulative event counts over time
    
    n_files = length(file_results);
    cum_files = cell(n_files, 1);
    
    all_event_times = [];
    max_time = 0;
    
    file_idx_valid = 0;
    
    for i = 1:n_files
        if isempty(file_results{i})
            continue;
        end
        
        file_result = file_results{i};
        
        if ~isfield(file_result, 'event_mask') || isempty(file_result.event_mask)
            continue;
        end
        
        event_mask = file_result.event_mask;
        [num_frames, num_rois] = size(event_mask);
        
        % Extract all event times for this file
        file_event_times = [];
        for roi = 1:num_rois
            roi_events = event_mask(:, roi);
            event_starts = find(diff([false; roi_events]) == 1);
            if ~isempty(event_starts)
                event_times_s = event_starts / frame_rate;
                file_event_times = [file_event_times; event_times_s];
            end
        end
        
        if ~isempty(file_event_times)
            % Sort event times
            file_event_times = sort(file_event_times);
            
            % Create cumulative curve
            n_events = length(file_event_times);
            cumulative = (1:n_events)';
            cumulative_normalized = cumulative / n_events;
            
            file_idx_valid = file_idx_valid + 1;
            cum_files{file_idx_valid} = struct();
            cum_files{file_idx_valid}.time = file_event_times;
            cum_files{file_idx_valid}.cumulative = cumulative;
            cum_files{file_idx_valid}.cumulative_normalized = cumulative_normalized;
            
            % Accumulate for mean
            all_event_times = [all_event_times; file_event_times];
            max_time = max(max_time, max(file_event_times));
        end
    end
    
    % Remove empty cells
    cum_files = cum_files(1:file_idx_valid);
    
    % Calculate mean cumulative curve
    if ~isempty(all_event_times)
        all_event_times = sort(all_event_times);
        n_events_total = length(all_event_times);
        
        cum_mean = struct();
        cum_mean.time = all_event_times;
        cum_mean.cumulative = (1:n_events_total)';
        cum_mean.cumulative_normalized = cum_mean.cumulative / n_events_total;
    else
        cum_mean = struct('time', [], 'cumulative', [], 'cumulative_normalized', []);
    end
    
    % Calculate statistics
    stats = struct();
    stats.n_files = file_idx_valid;
    
    if ~isempty(all_event_times) && length(all_event_times) > 10
        % Early vs late rate
        recording_duration = max_time;
        early_events = sum(all_event_times <= recording_duration / 3);
        late_events = sum(all_event_times > 2 * recording_duration / 3);
        
        stats.early_rate = early_events / (recording_duration / 3);
        stats.late_rate = late_events / (recording_duration / 3);
        stats.rate_change_percent = 100 * (stats.late_rate - stats.early_rate) / stats.early_rate;
    else
        stats.early_rate = 0;
        stats.late_rate = 0;
        stats.rate_change_percent = 0;
    end
end

function linearity = calculate_linearity(cum_data)
    % Calculate how linear the cumulative curve is (R²)
    
    if isempty(cum_data.time) || length(cum_data.time) < 3
        linearity.r_squared = NaN;
        linearity.slope = NaN;
        return;
    end
    
    % Fit linear model
    p = polyfit(cum_data.time, cum_data.cumulative_normalized, 1);
    y_fit = polyval(p, cum_data.time);
    
    % Calculate R²
    ss_res = sum((cum_data.cumulative_normalized - y_fit).^2);
    ss_tot = sum((cum_data.cumulative_normalized - mean(cum_data.cumulative_normalized)).^2);
    
    linearity.r_squared = 1 - (ss_res / ss_tot);
    linearity.slope = p(1);
end
