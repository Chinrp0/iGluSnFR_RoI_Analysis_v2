function plot_handles = simple_event_trace_plotter(raw_data, baseline, dfof_data, event_mask, event_stats, metadata, config)
    % SIMPLE_EVENT_TRACE_PLOTTER - Clean visualization of dF/F traces with event markings
    % Focus: Show baseline-corrected traces with clear event identification
    %
    % Inputs:
    %   raw_data    - [frames x ROIs] original fluorescence data
    %   baseline    - [frames x ROIs] calculated baseline
    %   dfof_data   - [frames x ROIs] dF/F traces
    %   event_mask  - [frames x ROIs] logical mask of detected events
    %   event_stats - Event detection statistics
    %   metadata    - File metadata
    %   config      - Configuration parameters
    
    [numFrames, numROIs] = size(dfof_data);
    time_vector = (1:numFrames) / config.frame_rate;  % Time in seconds
    
    plot_handles = struct();
    
    fprintf('Creating focused event trace plots for %s...\n', metadata.filename);
    
    %% === Plot 1: Sample dF/F Traces with Events ===
    plot_handles.event_traces = create_event_trace_overview(dfof_data, event_mask, event_stats, ...
        time_vector, metadata, config);
    
    %% === Plot 2: Event Detection Summary ===
    plot_handles.event_summary = create_event_detection_summary(event_stats, metadata, config);
    
    %% === Plot 3: High-Activity ROI Examples ===
    plot_handles.active_rois = create_active_roi_examples(dfof_data, event_mask, event_stats, ...
        time_vector, metadata, config);
    
    fprintf('  Created %d focused event plots\n', length(fieldnames(plot_handles)));
end

function fig = create_event_trace_overview(dfof_data, event_mask, event_stats, time_vector, metadata, config)
    % Show sample dF/F traces with events clearly marked
    
    [numFrames, numROIs] = size(dfof_data);
    
    % Select diverse ROIs for display
    sample_rois = select_representative_rois(event_stats, 12);  % Show 12 traces
    
    fig = figure('Name', sprintf('dF/F Traces with Events - %s', metadata.filename), ...
        'Position', [100, 100, 1600, 1000]);
    
    rows = 3;
    cols = 4;
    
    for i = 1:length(sample_rois)
        roi_idx = sample_rois(i);
        
        subplot(rows, cols, i);
        
        % Plot dF/F trace
        trace = dfof_data(:, roi_idx);
        plot(time_vector, trace, 'k-', 'LineWidth', 1.0);
        hold on;
        
        % Highlight detected events in red
        events = event_mask(:, roi_idx);
        if any(events)
            event_trace = trace;
            event_trace(~events) = NaN;  % Only show event portions
            plot(time_vector, event_trace, 'r-', 'LineWidth', 2.5);
        end
        
        % Add Schmitt trigger thresholds if available
        if isfield(event_stats, 'upper_thresholds') && length(event_stats.upper_thresholds) >= roi_idx
            upper_thresh = event_stats.upper_thresholds(roi_idx);
            lower_thresh = event_stats.lower_thresholds(roi_idx);
            
            yline(upper_thresh, 'g--', 'LineWidth', 1, 'Alpha', 0.6);
            yline(lower_thresh, 'c--', 'LineWidth', 1, 'Alpha', 0.6);
        end
        
        % Add zero reference line
        yline(0, 'k:', 'Alpha', 0.4, 'LineWidth', 0.5);
        
        % Get event count for this ROI
        num_events = event_stats.events_per_roi(roi_idx);
        max_dfof = max(trace, [], 'omitnan');
        
        xlabel('Time (s)', 'FontSize', 9);
        ylabel('dF/F', 'FontSize', 9);
        title(sprintf('ROI %d: %d events, max=%.3f', roi_idx, num_events, max_dfof), 'FontSize', 10);
        grid on;
        
        % Set reasonable y-limits
        y_range = [min(trace, [], 'omitnan'), max(trace, [], 'omitnan')];
        if diff(y_range) > 0
            y_margin = diff(y_range) * 0.1;
            ylim([y_range(1) - y_margin, y_range(2) + y_margin]);
        end
        
        % Add legend only to first subplot
        if i == 1
            legend('dF/F', 'Events', 'Upper (3.0σ)', 'Lower (1.5σ)', 'Location', 'best', 'FontSize', 8);
        end
    end
    
    sgtitle(sprintf('Baseline-Corrected Traces with Schmitt Trigger Events\n%s (%d total events)', ...
        metadata.filename, event_stats.total_events), 'FontSize', 14, 'FontWeight', 'bold');
end

function fig = create_event_detection_summary(event_stats, metadata, config)
    % Summary of event detection performance
    
    fig = figure('Name', sprintf('Event Detection Summary - %s', metadata.filename), ...
        'Position', [200, 100, 1000, 600]);
    
    % Subplot 1: Events per ROI histogram
    subplot(2, 2, 1);
    events_per_roi = event_stats.events_per_roi;
    histogram(events_per_roi, 0:max(events_per_roi), 'EdgeColor', 'none', 'FaceColor', [0.3, 0.6, 0.9]);
    
    % Add statistical lines
    mean_events = mean(events_per_roi);
    median_events = median(events_per_roi);
    xline(mean_events, 'r-', sprintf('Mean: %.1f', mean_events), 'LineWidth', 2);
    xline(median_events, 'g-', sprintf('Median: %.1f', median_events), 'LineWidth', 2);
    
    xlabel('Events per ROI');
    ylabel('ROI Count');
    title('Event Distribution');
    grid on;
    
    % Subplot 2: Threshold distribution
    subplot(2, 2, 2);
    if isfield(event_stats, 'upper_thresholds')
        histogram(event_stats.upper_thresholds, 30, 'EdgeColor', 'none', 'FaceColor', [0.6, 0.9, 0.4], 'FaceAlpha', 0.7);
        hold on;
        histogram(event_stats.lower_thresholds, 30, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.9, 0.6], 'FaceAlpha', 0.7);
        
        xlabel('Threshold Value (dF/F)');
        ylabel('ROI Count');
        title('Schmitt Trigger Thresholds');
        legend('Upper (3.0σ)', 'Lower (1.5σ)', 'Location', 'best');
        grid on;
    else
        text(0.5, 0.5, 'Threshold data not available', 'HorizontalAlignment', 'center');
        title('Thresholds - NO DATA');
    end
    
    % Subplot 3: Activity classification
    subplot(2, 2, 3);
    
    % Classify ROIs by activity level
    no_events = sum(events_per_roi == 0);
    low_activity = sum(events_per_roi >= 1 & events_per_roi <= 5);
    medium_activity = sum(events_per_roi >= 6 & events_per_roi <= 15);
    high_activity = sum(events_per_roi > 15);
    
    activity_labels = {'No Events', 'Low (1-5)', 'Medium (6-15)', 'High (>15)'};
    activity_counts = [no_events, low_activity, medium_activity, high_activity];
    
    pie(activity_counts, activity_labels);
    title('ROI Activity Classification');
    
    % Subplot 4: Summary statistics
    subplot(2, 2, 4);
    axis off;
    
    summary_text = {
        sprintf('Dataset: %s', metadata.filename),
        '',
        sprintf('Total ROIs: %d', length(events_per_roi)),
        sprintf('Active ROIs: %d (%.1f%%)', event_stats.rois_with_events, ...
            100 * event_stats.rois_with_events / length(events_per_roi)),
        sprintf('Total Events: %d', event_stats.total_events),
        sprintf('Events/Active ROI: %.1f', event_stats.total_events / max(1, event_stats.rois_with_events)),
        '',
        'Schmitt Trigger Settings:',
        sprintf('Upper Threshold: %.1fσ', config.event_detection.upper_threshold_sigma),
        sprintf('Lower Threshold: %.1fσ', config.event_detection.lower_threshold_sigma),
        sprintf('Decay Extension: %d frames (%.0f ms)', ...
            config.event_detection.decay_extension_frames, ...
            config.event_detection.decay_extension_frames * 1000 / config.frame_rate),
        '',
        'Performance:',
        sprintf('Mean Duration: %.1f frames (%.0f ms)', ...
            event_stats.mean_event_duration, ...
            event_stats.mean_event_duration * 1000 / config.frame_rate),
        sprintf('Event Coverage: %.1f%% of total frames', ...
            100 * event_stats.fraction_frames_in_events)
    };
    
    text(0.05, 0.95, summary_text, 'FontSize', 10, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', 'Units', 'normalized');
    
    sgtitle('Schmitt Trigger Event Detection Summary', 'FontSize', 14, 'FontWeight', 'bold');
end

function fig = create_active_roi_examples(dfof_data, event_mask, event_stats, time_vector, metadata, config)
    % Show examples of most active ROIs with detailed event markings
    
    % Select top 8 most active ROIs
    [~, active_idx] = sort(event_stats.events_per_roi, 'descend');
    top_active = active_idx(1:min(8, length(active_idx)));
    
    fig = figure('Name', sprintf('Most Active ROIs - %s', metadata.filename), ...
        'Position', [300, 100, 1400, 900]);
    
    rows = 2;
    cols = 4;
    
    for i = 1:length(top_active)
        roi_idx = top_active(i);
        
        subplot(rows, cols, i);
        
        % Plot dF/F trace
        trace = dfof_data(:, roi_idx);
        plot(time_vector, trace, 'k-', 'LineWidth', 0.8);
        hold on;
        
        % Highlight events with different colors for different event episodes
        events = event_mask(:, roi_idx);
        if any(events)
            % Find individual event episodes
            event_starts = find(diff([false; events]) == 1);
            event_ends = find(diff([events; false]) == -1);
            
            % Color each event episode differently
            colors = lines(length(event_starts));
            
            for e = 1:length(event_starts)
                event_frames = event_starts(e):event_ends(e);
                event_trace = NaN(size(trace));
                event_trace(event_frames) = trace(event_frames);
                
                plot(time_vector, event_trace, 'Color', colors(e, :), 'LineWidth', 3);
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
        
        xlabel('Time (s)', 'FontSize', 9);
        ylabel('dF/F', 'FontSize', 9);
        title(sprintf('ROI %d: %d events, max=%.3f', roi_idx, num_events, max_response), 'FontSize', 10);
        grid on;
        
        % Set y-limits to show events clearly
        y_range = [min(trace, [], 'omitnan'), max(trace, [], 'omitnan')];
        if diff(y_range) > 0
            y_margin = diff(y_range) * 0.15;
            ylim([y_range(1) - y_margin, y_range(2) + y_margin]);
        end
    end
    
    sgtitle(sprintf('Most Active ROIs - Individual Event Episodes\n%s', metadata.filename), ...
        'FontSize', 14, 'FontWeight', 'bold');
end

function sample_rois = select_representative_rois(event_stats, num_samples)
    % Select ROIs representing different activity levels
    
    events_per_roi = event_stats.events_per_roi;
    numROIs = length(events_per_roi);
    
    if numROIs <= num_samples
        sample_rois = 1:numROIs;
        return;
    end
    
    % Sort ROIs by activity level
    [sorted_events, sort_idx] = sort(events_per_roi, 'descend');
    
    % Select ROIs from different activity quartiles
    sample_rois = [];
    
    % High activity (top 25%)
    high_end = ceil(length(sort_idx) * 0.25);
    n_high = ceil(num_samples * 0.4);
    sample_rois = [sample_rois, sort_idx(1:min(n_high, high_end))];
    
    % Medium activity (25-75%)
    med_start = high_end + 1;
    med_end = ceil(length(sort_idx) * 0.75);
    n_med = ceil(num_samples * 0.4);
    med_candidates = sort_idx(med_start:med_end);
    if ~isempty(med_candidates)
        step = max(1, floor(length(med_candidates) / n_med));
        sample_rois = [sample_rois, med_candidates(1:step:end)];
    end
    
    % Low activity (bottom 25%)
    low_start = med_end + 1;
    n_low = num_samples - length(sample_rois);
    low_candidates = sort_idx(low_start:end);
    if ~isempty(low_candidates) && n_low > 0
        sample_rois = [sample_rois, low_candidates(1:min(n_low, length(low_candidates)))];
    end
    
    % Ensure we have exactly num_samples
    sample_rois = sample_rois(1:min(num_samples, length(sample_rois)));
    
    % Fill with random ROIs if needed
    if length(sample_rois) < num_samples
        remaining = setdiff(1:numROIs, sample_rois);
        n_fill = num_samples - length(sample_rois);
        sample_rois = [sample_rois, remaining(1:min(n_fill, length(remaining)))];
    end
    
    sample_rois = sort(sample_rois);
end