function plot_handles = baseline_plotter_clean(raw_data, baseline, dfof_data, outlier_mask, baseline_stats, dfof_stats, metadata, config, event_stats)
    % BASELINE_PLOTTER_CLEAN - Version without peak markers
    % Red highlighted events are sufficient for visualization
    % Maintains 2x4 layout as requested
    
    % Handle backward compatibility
    if nargin < 9
        event_stats = create_empty_event_stats(size(raw_data));
        fprintf('Warning: No event_stats provided, using empty event data\n');
    end
    
    % Extract event_mask consistently
    if isfield(event_stats, 'event_mask')
        event_mask = event_stats.event_mask;
    else
        fprintf('Warning: event_mask not found in event_stats, creating empty mask\n');
        event_mask = false(size(raw_data));
    end
    
    [numFrames, numROIs] = size(raw_data);
    time_vector = (1:numFrames) / config.frame_rate;
    
    plot_handles = struct();
    
    fprintf('Creating clean event plots (no peak markers) for %s...\n', metadata.filename);
    
    %% === Plot 1: Clean dF/F Traces (2x4 layout) ===
    plot_handles.event_traces = create_clean_event_traces(dfof_data, event_mask, event_stats, ...
        time_vector, metadata, config);
    
    %% === Plot 2: Event Detection Summary ===
    plot_handles.event_summary = create_event_summary(event_stats, metadata, config);
    
    %% === Plot 3: Most Active ROIs (2x4 layout) ===
    plot_handles.active_rois = create_clean_active_roi_examples(dfof_data, event_mask, event_stats, ...
        time_vector, metadata, config);
    
    fprintf('  Created %d clean plots (2x4 layout, no peak markers)\n', length(fieldnames(plot_handles)));
end

function fig = create_clean_event_traces(dfof_data, event_mask, event_stats, time_vector, metadata, config)
    % CLEAN: 2x4 layout with red events, NO peak markers
    
    [numFrames, numROIs] = size(dfof_data);
    
    % Select 8 representative ROIs for 2x4 layout
    sample_rois = select_diverse_rois(event_stats, 8);
    
    fig = figure('Name', sprintf('dF/F Traces with Events (Clean) - %s', metadata.filename), ...
        'Position', [100, 100, 1000, 1200]);  % Optimized for 2x4 layout
    
    % STANDARDIZED: 2 columns x 4 rows = 8 subplots
    rows = 4;
    cols = 2;
    
    for i = 1:length(sample_rois)
        roi_idx = sample_rois(i);
        
        subplot(rows, cols, i);
        
        % Plot dF/F trace
        trace = dfof_data(:, roi_idx);
        plot(time_vector, trace, 'k-', 'LineWidth', 1.0);
        hold on;
        
        % Highlight detected events in red (NO PEAK MARKERS)
        if size(event_mask, 2) >= roi_idx
            events = event_mask(:, roi_idx);
            if any(events)
                % Only highlight event periods in red - no additional markers
                event_trace = trace;
                event_trace(~events) = NaN;
                plot(time_vector, event_trace, 'r-', 'LineWidth', 2.5);
            end
        end
        
        % Add Schmitt trigger thresholds
        if isfield(event_stats, 'upper_thresholds') && length(event_stats.upper_thresholds) >= roi_idx
            upper_thresh = event_stats.upper_thresholds(roi_idx);
            lower_thresh = event_stats.lower_thresholds(roi_idx);
            
            yline(upper_thresh, 'g--', 'LineWidth', 1, 'Alpha', 0.7);
            yline(lower_thresh, 'c--', 'LineWidth', 1, 'Alpha', 0.7);
        end
        
        % Add zero reference
        yline(0, 'k:', 'Alpha', 0.5, 'LineWidth', 0.8);
        
        % Get event count and noise info
        if isfield(event_stats, 'events_per_roi') && length(event_stats.events_per_roi) >= roi_idx
            num_events = event_stats.events_per_roi(roi_idx);
        else
            num_events = 0;
        end
        
        % Get noise level if available
        if isfield(event_stats, 'noise_metrics') && isfield(event_stats.noise_metrics, 'noise_std')
            noise_level = event_stats.noise_metrics.noise_std(roi_idx);
            method_used = event_stats.noise_metrics.noise_method{roi_idx};
            snr = event_stats.noise_metrics.signal_quality(roi_idx);
            
            title_text = sprintf('ROI %d: %d events, σ=%.5f\nSNR=%.1f (%s)', ...
                roi_idx, num_events, noise_level, snr, method_used);
        else
            max_dfof = max(trace, [], 'omitnan');
            title_text = sprintf('ROI %d: %d events, max=%.3f', roi_idx, num_events, max_dfof);
        end
        
        xlabel('Time (s)', 'FontSize', 10);
        ylabel('dF/F', 'FontSize', 10);
        title(title_text, 'FontSize', 10, 'FontWeight', 'bold');
        grid on; grid minor;
        
        % Set consistent y-limits (less margin since no peak markers)
        y_range = [min(trace, [], 'omitnan'), max(trace, [], 'omitnan')];
        if diff(y_range) > 0
            y_margin = diff(y_range) * 0.15;  % Reduced margin
            ylim([y_range(1) - y_margin, y_range(2) + y_margin]);
        end
        
        % Add legend only to first subplot
        if i == 1 && any(event_mask(:, roi_idx))
            legend('dF/F', 'Events', 'Upper (3.0σ)', 'Lower (1.5σ)', ...
                'Location', 'best', 'FontSize', 9);
        elseif i == 1
            legend('dF/F', 'Upper (3.0σ)', 'Lower (1.5σ)', ...
                'Location', 'best', 'FontSize', 9);
        end
    end
    
    total_events = 0;
    if isfield(event_stats, 'total_events')
        total_events = event_stats.total_events;
    end
    
    sgtitle(sprintf('Clean Event Detection - Red Highlights Only\n%s (%d total events)', ...
        metadata.filename, total_events), 'FontSize', 14, 'FontWeight', 'bold');
end

function fig = create_clean_active_roi_examples(dfof_data, event_mask, event_stats, time_vector, metadata, config)
    % CLEAN: Most active ROIs in 2x4 layout without peak markers
    
    if ~isfield(event_stats, 'events_per_roi') || isempty(event_stats.events_per_roi)
        fig = create_empty_plot(sprintf('Most Active ROIs - %s', metadata.filename), 'No event data available');
        return;
    end
    
    % Select top 8 most active ROIs for 2x4 layout
    [~, active_idx] = sort(event_stats.events_per_roi, 'descend');
    top_active = active_idx(1:min(8, length(active_idx)));
    
    fig = figure('Name', sprintf('Most Active ROIs (Clean) - %s', metadata.filename), ...
        'Position', [300, 100, 1000, 1200]);  % Optimized for 2x4
    
    % STANDARDIZED: 2 columns x 4 rows
    rows = 4;
    cols = 2;
    
    for i = 1:length(top_active)
        roi_idx = top_active(i);
        
        subplot(rows, cols, i);
        
        % Plot dF/F trace
        trace = dfof_data(:, roi_idx);
        plot(time_vector, trace, 'k-', 'LineWidth', 0.8);
        hold on;
        
        % Highlight events WITHOUT peak markers
        if size(event_mask, 2) >= roi_idx
            events = event_mask(:, roi_idx);
            if any(events)
                % Color individual event episodes differently for variety
                event_starts = find(diff([false; events]) == 1);
                event_ends = find(diff([events; false]) == -1);
                
                colors = lines(length(event_starts));
                
                for e = 1:length(event_starts)
                    event_frames = event_starts(e):event_ends(e);
                    event_trace = NaN(size(trace));
                    event_trace(event_frames) = trace(event_frames);
                    
                    plot(time_vector, event_trace, 'Color', colors(e, :), 'LineWidth', 2.5);
                end
                
                % NO PEAK MARKERS - just the colored event traces
            end
        end
        
        % Add thresholds
        if isfield(event_stats, 'upper_thresholds') && length(event_stats.upper_thresholds) >= roi_idx
            upper_thresh = event_stats.upper_thresholds(roi_idx);
            lower_thresh = event_stats.lower_thresholds(roi_idx);
            
            yline(upper_thresh, 'g--', 'LineWidth', 1.2, 'Alpha', 0.7);
            yline(lower_thresh, 'c--', 'LineWidth', 1.2, 'Alpha', 0.7);
        end
        
        yline(0, 'k:', 'Alpha', 0.5);
        
        num_events = event_stats.events_per_roi(roi_idx);
        max_response = max(trace, [], 'omitnan');
        
        % Enhanced title with noise info if available
        if isfield(event_stats, 'noise_metrics')
            noise_level = event_stats.noise_metrics.noise_std(roi_idx);
            snr = event_stats.noise_metrics.signal_quality(roi_idx);
            title_text = sprintf('ROI %d: %d events, σ=%.5f\nSNR=%.1f, max=%.3f', ...
                roi_idx, num_events, noise_level, snr, max_response);
        else
            title_text = sprintf('ROI %d: %d events, max=%.3f', roi_idx, num_events, max_response);
        end
        
        xlabel('Time (s)', 'FontSize', 10);
        ylabel('dF/F', 'FontSize', 10);
        title(title_text, 'FontSize', 10);
        grid on;
        
        % Set y-limits (reduced margin since no peak markers)
        y_range = [min(trace, [], 'omitnan'), max(trace, [], 'omitnan')];
        if diff(y_range) > 0
            y_margin = diff(y_range) * 0.15;  % Reduced margin
            ylim([y_range(1) - y_margin, y_range(2) + y_margin]);
        end
    end
    
    sgtitle(sprintf('Most Active ROIs - Clean Event Visualization\n%s', metadata.filename), ...
        'FontSize', 14, 'FontWeight', 'bold');
end

function fig = create_event_summary(event_stats, metadata, config)
    % Event detection summary (unchanged - no peak markers here anyway)
    
    fig = figure('Name', sprintf('Event Detection Summary - %s', metadata.filename), ...
        'Position', [200, 100, 1000, 600]);
    
    if ~isfield(event_stats, 'events_per_roi') || isempty(event_stats.events_per_roi)
        text(0.5, 0.5, 'No event data available for summary', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'FontSize', 14);
        title('Event Summary - NO DATA');
        return;
    end
    
    events_per_roi = event_stats.events_per_roi;
    
    % Create 2x2 summary layout
    subplot(2, 2, 1);
    histogram(events_per_roi, 0:max(events_per_roi), 'EdgeColor', 'none', 'FaceColor', [0.3, 0.6, 0.9]);
    mean_events = mean(events_per_roi);
    median_events = median(events_per_roi);
    xline(mean_events, 'r-', sprintf('Mean: %.1f', mean_events), 'LineWidth', 2);
    xline(median_events, 'g-', sprintf('Median: %.1f', median_events), 'LineWidth', 2);
    xlabel('Events per ROI'); ylabel('ROI Count'); title('Event Distribution'); grid on;
    
    subplot(2, 2, 2);
    if isfield(event_stats, 'upper_thresholds')
        histogram(event_stats.upper_thresholds, 30, 'EdgeColor', 'none', 'FaceColor', [0.6, 0.9, 0.4], 'FaceAlpha', 0.7);
        hold on;
        histogram(event_stats.lower_thresholds, 30, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.9, 0.6], 'FaceAlpha', 0.7);
        xlabel('Threshold Value (dF/F)'); ylabel('ROI Count'); title('Schmitt Trigger Thresholds');
        legend('Upper (3.0σ)', 'Lower (1.5σ)', 'Location', 'best'); grid on;
    else
        text(0.5, 0.5, 'Threshold data not available', 'HorizontalAlignment', 'center');
        title('Thresholds - NO DATA');
    end
    
    subplot(2, 2, 3);
    no_events = sum(events_per_roi == 0);
    low_activity = sum(events_per_roi >= 1 & events_per_roi <= 5);
    medium_activity = sum(events_per_roi >= 6 & events_per_roi <= 15);
    high_activity = sum(events_per_roi > 15);
    
    activity_labels = {'No Events', 'Low (1-5)', 'Medium (6-15)', 'High (>15)'};
    activity_counts = [no_events, low_activity, medium_activity, high_activity];
    pie(activity_counts, activity_labels); title('ROI Activity Classification');
    
    subplot(2, 2, 4); axis off;
    summary_text = {
        sprintf('Dataset: %s', metadata.filename);
        sprintf('Total ROIs: %d', length(events_per_roi));
        sprintf('Active ROIs: %d (%.1f%%)', event_stats.rois_with_events, ...
            100 * event_stats.rois_with_events / length(events_per_roi));
        sprintf('Total Events: %d', event_stats.total_events);
        '';
        'Detection Settings:';
        sprintf('Upper: %.1fσ, Lower: %.1fσ', config.event_detection.upper_threshold_sigma, ...
            config.event_detection.lower_threshold_sigma);
        sprintf('Decay Extension: %d frames', config.event_detection.decay_extension_frames);
    };
    
    text(0.05, 0.95, summary_text, 'FontSize', 10, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', 'Units', 'normalized');
    
    sgtitle('Event Detection Summary', 'FontSize', 14, 'FontWeight', 'bold');
end

% === HELPER FUNCTIONS (unchanged) ===

function sample_rois = select_diverse_rois(event_stats, num_samples)
    % Select ROIs representing different activity levels for display
    
    if ~isfield(event_stats, 'events_per_roi') || isempty(event_stats.events_per_roi)
        sample_rois = 1:min(num_samples, 100);  % Default fallback
        return;
    end
    
    events_per_roi = event_stats.events_per_roi;
    numROIs = length(events_per_roi);
    
    if numROIs <= num_samples
        sample_rois = 1:numROIs;
        return;
    end
    
    % Sort ROIs by activity level
    [~, sort_idx] = sort(events_per_roi, 'descend');
    
    % Select diverse ROIs across activity spectrum
    high_activity = sort_idx(1:ceil(length(sort_idx)*0.2));      % Top 20%
    med_activity = sort_idx(ceil(length(sort_idx)*0.3):ceil(length(sort_idx)*0.7));  % Middle 40%
    low_activity = sort_idx(ceil(length(sort_idx)*0.8):end);     % Bottom 20%
    
    % Distribute samples across categories
    n_high = ceil(num_samples * 0.5);  % 50% high activity
    n_med = ceil(num_samples * 0.3);   % 30% medium activity  
    n_low = num_samples - n_high - n_med;  % 20% low activity
    
    sample_rois = [high_activity(1:min(n_high, length(high_activity))), ...
                   med_activity(1:min(n_med, length(med_activity))), ...
                   low_activity(1:min(n_low, length(low_activity)))];
    
    % Fill with any remaining if needed
    if length(sample_rois) < num_samples
        remaining = setdiff(1:numROIs, sample_rois);
        needed = num_samples - length(sample_rois);
        sample_rois = [sample_rois, remaining(1:min(needed, length(remaining)))];
    end
    
    sample_rois = sort(sample_rois(1:num_samples));
end

function event_stats = create_empty_event_stats(data_size)
    % Create empty event stats for backward compatibility
    
    [numFrames, numROIs] = data_size;
    
    event_stats = struct();
    event_stats.events_per_roi = zeros(1, numROIs);
    event_stats.total_events = 0;
    event_stats.rois_with_events = 0;
    event_stats.event_mask = false(numFrames, numROIs);
    event_stats.upper_thresholds = zeros(1, numROIs);
    event_stats.lower_thresholds = zeros(1, numROIs);
    
    event_stats.event_summary = struct();
    event_stats.event_summary.total_events = 0;
    event_stats.event_summary.rois_with_events = 0;
    event_stats.event_summary.fraction_active_rois = 0;
end

function fig = create_empty_plot(fig_name, message)
    % Create empty plot with message
    
    fig = figure('Name', fig_name);
    text(0.5, 0.5, message, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontSize', 14);
    title(fig_name);
end