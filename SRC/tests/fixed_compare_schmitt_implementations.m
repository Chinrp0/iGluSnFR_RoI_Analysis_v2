%% FIXED COMPARE PURE SCHMITT TRIGGER
% Compare current vs pure Schmitt trigger (no artificial extensions)
% FIXED: MATLAB compatibility issues and threshold logic

clear; clc; close all;

fprintf('=== Comparing FIXED Pure Schmitt Trigger Implementation ===\n');

%% === Load Data ===
addpath(genpath(pwd));
test_file = 'D:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglusnfr4f_NGR\Spont\GPU_SNR_Processed\5_raw_mean\CP_Snfr4-NGR_Doc2b-R213W_Cs1-c1_spont-01_mean.csv';

loader = csv_loader_v2();
[data, metadata] = loader.loadSingleFile(test_file);
config = tracenorm_config();

[baseline, outlier_mask, baseline_stats] = baseline_detector(data, config);
[dfof_data, dfof_stats] = dfof_calculator(data, baseline, config);

% Use the same ROI that showed the problem
demo_roi = 610;
roi_trace = dfof_data(:, demo_roi);
time_vector = (1:length(roi_trace)) / config.frame_rate;

fprintf('\nAnalyzing ROI %d:\n', demo_roi);

%% === Test Current Implementation ===
fprintf('\n1. CURRENT implementation (with extensions):\n');
try
    [event_mask_current, event_stats_current] = optimal_schmitt_event_detector(dfof_data, config, baseline_stats);
    
    current_events = event_mask_current(:, demo_roi);
    current_starts = find(diff([false; current_events]) == 1);
    current_ends = find(diff([current_events; false]) == -1);
    
    fprintf('   ROI %d events with CURRENT method:\n', demo_roi);
    for i = 1:min(5, length(current_starts))
        duration = current_ends(i) - current_starts(i) + 1;
        fprintf('     Event %d: [%d-%d] = %d frames (%.0f ms)\n', ...
            i, current_starts(i), current_ends(i), duration, duration * 1000 / config.frame_rate);
    end
catch ME
    fprintf('   Error with current method: %s\n', ME.message);
    current_events = false(size(roi_trace));
    current_starts = [];
    current_ends = [];
    event_stats_current = struct();
end

%% === Test FIXED Pure Schmitt Implementation ===
fprintf('\n2. FIXED PURE SCHMITT implementation:\n');

% Calculate thresholds using research paper method
valid_data = roi_trace(~isnan(roi_trace));
if length(valid_data) > 20
    baseline_70th = prctile(valid_data, 70);
    baseline_data = valid_data(valid_data <= baseline_70th);
    noise_std = mad(baseline_data, 1) * 1.4826;
else
    noise_std = mad(valid_data, 1) * 1.4826;
end

upper_threshold = config.event_detection.upper_threshold_sigma * noise_std;
lower_threshold = config.event_detection.lower_threshold_sigma * noise_std;

% Apply manual Schmitt trigger
pure_events = apply_manual_schmitt_trigger(roi_trace, upper_threshold, lower_threshold);

pure_starts = find(diff([false; pure_events]) == 1);
pure_ends = find(diff([pure_events; false]) == -1);

fprintf('   ROI %d events with FIXED PURE SCHMITT method:\n', demo_roi);
for i = 1:min(5, length(pure_starts))
    duration = pure_ends(i) - pure_starts(i) + 1;
    fprintf('     Event %d: [%d-%d] = %d frames (%.0f ms)\n', ...
        i, pure_starts(i), pure_ends(i), duration, duration * 1000 / config.frame_rate);
end

%% === Debug Threshold Calculation ===
fprintf('\n3. Debug threshold calculation for ROI %d:\n', demo_roi);

fprintf('   Noise std: %.6f\n', noise_std);
fprintf('   Upper threshold (%.1fσ): %.6f\n', config.event_detection.upper_threshold_sigma, upper_threshold);
fprintf('   Lower threshold (%.1fσ): %.6f\n', config.event_detection.lower_threshold_sigma, lower_threshold);
fprintf('   Hysteresis gap: %.6f\n', upper_threshold - lower_threshold);

% Show trace statistics
fprintf('   Trace range: [%.6f, %.6f]\n', min(valid_data), max(valid_data));
fprintf('   Frames above upper: %d\n', sum(valid_data > upper_threshold));
fprintf('   Frames above lower: %d\n', sum(valid_data > lower_threshold));

%% === Visual Comparison ===
fprintf('\n4. Creating visual comparison...\n');

fig = figure('Name', sprintf('FIXED Pure Schmitt Trigger - ROI %d', demo_roi), ...
    'Position', [100, 100, 1600, 1200]);

% Plot 1: Full trace comparison
subplot(3, 1, 1);
plot(time_vector, roi_trace, 'k-', 'LineWidth', 1.0);
hold on;

% Highlight current events (with extensions)
if any(current_events)
    current_trace = roi_trace;
    current_trace(~current_events) = NaN;
    plot(time_vector, current_trace, 'r-', 'LineWidth', 2.5);
end

% Highlight pure events (no extensions)
if any(pure_events)
    pure_trace = roi_trace;
    pure_trace(~pure_events) = NaN;
    plot(time_vector, pure_trace, 'b-', 'LineWidth', 2.0);
end

% Add threshold lines
line([time_vector(1), time_vector(end)], [upper_threshold, upper_threshold], ...
    'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.5);
line([time_vector(1), time_vector(end)], [lower_threshold, lower_threshold], ...
    'Color', 'c', 'LineStyle', '--', 'LineWidth', 1.5);
line([time_vector(1), time_vector(end)], [0, 0], ...
    'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 1);

xlabel('Time (s)'); ylabel('dF/F');
title(sprintf('ROI %d: Current (Red, %d events) vs Fixed Pure Schmitt (Blue, %d events)', ...
    demo_roi, length(current_starts), length(pure_starts)));
legend('dF/F', 'Current (Extended)', 'Fixed Pure Schmitt', 'Upper', 'Lower', 'Zero', 'Location', 'best');
grid on;

% Plot 2: Zoomed view
subplot(3, 1, 2);
if ~isempty(pure_starts) && length(pure_starts) >= 1
    % Zoom around first pure event
    first_event = pure_starts(1);
    zoom_start = max(1, first_event - 20);
    zoom_end = min(length(roi_trace), first_event + 50);
    zoom_indices = zoom_start:zoom_end;
    
    plot(time_vector(zoom_indices), roi_trace(zoom_indices), 'k-', 'LineWidth', 1.5);
    hold on;
    
    % Current events in zoom
    if any(current_events)
        current_zoom = current_events(zoom_indices);
        if any(current_zoom)
            current_trace_zoom = roi_trace(zoom_indices);
            current_trace_zoom(~current_zoom) = NaN;
            plot(time_vector(zoom_indices), current_trace_zoom, 'r-', 'LineWidth', 3);
        end
    end
    
    % Pure events in zoom
    pure_zoom = pure_events(zoom_indices);
    if any(pure_zoom)
        pure_trace_zoom = roi_trace(zoom_indices);
        pure_trace_zoom(~pure_zoom) = NaN;
        plot(time_vector(zoom_indices), pure_trace_zoom, 'b-', 'LineWidth', 2);
    end
    
    % Add threshold lines to zoom
    line([time_vector(zoom_start), time_vector(zoom_end)], [upper_threshold, upper_threshold], ...
        'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.5);
    line([time_vector(zoom_start), time_vector(zoom_end)], [lower_threshold, lower_threshold], ...
        'Color', 'c', 'LineStyle', '--', 'LineWidth', 1.5);
    line([time_vector(zoom_start), time_vector(zoom_end)], [0, 0], ...
        'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 1);
    
    xlabel('Time (s)'); ylabel('dF/F');
    title('Zoomed View: Fixed Schmitt Trigger Behavior');
    legend('dF/F', 'Current (Extended)', 'Fixed Pure Schmitt', 'Upper', 'Lower', 'Zero', 'Location', 'best');
    grid on;
else
    text(0.5, 0.5, 'No pure events to zoom', 'HorizontalAlignment', 'center');
    title('No Pure Events Detected');
end

% Plot 3: Event duration comparison (FIXED MATLAB ERROR)
subplot(3, 1, 3);

% Get all event durations - FIXED: Handle empty arrays properly
current_durations = [];
if ~isempty(current_starts)
    current_durations = current_ends - current_starts + 1;
end

pure_durations = [];
if ~isempty(pure_starts)
    pure_durations = pure_ends - pure_starts + 1;
end

% Create histogram comparison - FIXED: Handle dimension mismatch
if ~isempty(current_durations) || ~isempty(pure_durations)
    % Determine max duration safely
    max_duration = 20;  % Default minimum
    if ~isempty(current_durations)
        max_duration = max(max_duration, max(current_durations));
    end
    if ~isempty(pure_durations)
        max_duration = max(max_duration, max(pure_durations));
    end
    
    bins = 1:max_duration;
    
    % Create histograms manually to avoid compatibility issues
    if ~isempty(current_durations)
        [counts1, ~] = hist(current_durations, bins);
        bar(bins, counts1, 'FaceAlpha', 0.5, 'FaceColor', 'r', 'EdgeColor', 'none');
        hold on;
    end
    
    if ~isempty(pure_durations)
        [counts2, ~] = hist(pure_durations, bins);
        bar(bins, counts2, 'FaceAlpha', 0.5, 'FaceColor', 'b', 'EdgeColor', 'none');
    end
    
    xlabel('Event Duration (frames)'); ylabel('Count');
    title(sprintf('Event Duration Distribution - ROI %d', demo_roi));
    
    legend_labels = {};
    if ~isempty(current_durations)
        legend_labels{end+1} = 'Current (Extended)';
    end
    if ~isempty(pure_durations)
        legend_labels{end+1} = 'Fixed Pure Schmitt';
    end
    if ~isempty(legend_labels)
        legend(legend_labels, 'Location', 'best');
    end
    
    grid on;
else
    text(0.5, 0.5, 'No events detected', 'HorizontalAlignment', 'center');
    title('No Events for Duration Analysis');
end

%% === Compare All Methods ===
fprintf('\n=== COMPARISON SUMMARY ===\n');
fprintf('Method                     Events    Avg Duration    Total Frames\n');
fprintf('-------------------------  --------  --------------  ------------\n');

current_avg_dur = 0;
if ~isempty(current_starts)
    current_avg_dur = mean(current_ends - current_starts + 1);
end

pure_avg_dur = 0;
if ~isempty(pure_starts)
    pure_avg_dur = mean(pure_ends - pure_starts + 1);
end

fprintf('Current (with extensions)  %8d      %8.1f        %8d\n', ...
    length(current_starts), current_avg_dur, sum(current_events));

fprintf('Fixed Pure Schmitt         %8d      %8.1f        %8d\n', ...
    length(pure_starts), pure_avg_dur, sum(pure_events));

fprintf('\n✅ Fixed comparison complete!\n');

%% === HELPER FUNCTION ===
function event_mask = apply_manual_schmitt_trigger(trace, upper_thresh, lower_thresh)
    % Manual Schmitt trigger implementation for single ROI
    
    event_mask = false(size(trace));
    state = 'baseline';
    event_start = 0;
    raw_events = [];
    
    % Phase 1: Pure Schmitt trigger
    for frame = 1:length(trace)
        signal = trace(frame);
        
        if isnan(signal)
            continue;
        end
        
        if strcmp(state, 'baseline')
            if signal > upper_thresh
                state = 'in_event';
                event_start = frame;
            end
        elseif strcmp(state, 'in_event')
            if signal < lower_thresh
                event_end = frame - 1;  % Last frame above lower threshold
                if event_end >= event_start
                    raw_events(end+1, :) = [event_start, event_end];
                end
                state = 'baseline';
                event_start = 0;
            end
        end
    end
    
    % Handle event extending to end
    if strcmp(state, 'in_event') && event_start > 0
        raw_events(end+1, :) = [event_start, length(trace)];
    end
    
    if isempty(raw_events)
        return;
    end
    
    % Phase 2: Merge nearby events (≤2 frame gap)
    merged_events = [];
    current_start = raw_events(1, 1);
    current_end = raw_events(1, 2);
    
    for i = 2:size(raw_events, 1)
        next_start = raw_events(i, 1);
        next_end = raw_events(i, 2);
        
        gap = next_start - current_end - 1;
        if gap <= 2  % Merge nearby events
            current_end = next_end;
        else
            merged_events(end+1, :) = [current_start, current_end];
            current_start = next_start;
            current_end = next_end;
        end
    end
    merged_events(end+1, :) = [current_start, current_end];
    
    % Phase 3: Apply minimum duration filter (≥3 frames)
    for i = 1:size(merged_events, 1)
        start_frame = merged_events(i, 1);
        end_frame = merged_events(i, 2);
        duration = end_frame - start_frame + 1;
        
        if duration >= 3  % Minimum duration filter
            event_mask(start_frame:end_frame) = true;
        end
    end
end
