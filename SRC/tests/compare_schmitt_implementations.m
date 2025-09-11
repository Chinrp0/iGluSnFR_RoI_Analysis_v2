%% COMPARE PURE SCHMITT TRIGGER
% Compare current vs pure Schmitt trigger (no artificial extensions)
% Fixed MATLAB compatibility issues

clear; clc; close all;

fprintf('=== Comparing Pure Schmitt Trigger Implementation ===\n');

%% === Load Data ===
addpath(genpath(pwd));
test_file = 'D:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglusnfr4f_NGR\Spont\GPU_SNR_Processed\5_raw_mean\CP_Snfr4-NGR_Doc2b-R213W_Cs1-c1_spont-01_mean.csv';

loader = csv_loader_v2();
[data, metadata] = loader.loadSingleFile(test_file);
config = tracenorm_config();

[baseline, outlier_mask, baseline_stats] = baseline_detector(data, config);
[dfof_data, dfof_stats] = dfof_calculator(data, baseline, config);

% Use the same ROI that showed the problem
demo_roi = 613;
roi_trace = dfof_data(:, demo_roi);
time_vector = (1:length(roi_trace)) / config.frame_rate;

fprintf('\nAnalyzing ROI %d:\n', demo_roi);

%% === Test Current Implementation ===
fprintf('\n1. CURRENT implementation (with extensions):\n');
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

%% === Test Pure Schmitt Implementation ===
fprintf('\n2. PURE SCHMITT implementation (no extensions):\n');
[event_mask_pure, event_stats_pure] = pure_schmitt_trigger_detector(dfof_data, config, baseline_stats);

pure_events = event_mask_pure(:, demo_roi);
pure_starts = find(diff([false; pure_events]) == 1);
pure_ends = find(diff([pure_events; false]) == -1);

fprintf('   ROI %d events with PURE SCHMITT method:\n', demo_roi);
for i = 1:min(5, length(pure_starts))
    duration = pure_ends(i) - pure_starts(i) + 1;
    fprintf('     Event %d: [%d-%d] = %d frames (%.0f ms)\n', ...
        i, pure_starts(i), pure_ends(i), duration, duration * 1000 / config.frame_rate);
end

%% === Visual Comparison ===
fprintf('\n3. Creating visual comparison...\n');

fig = figure('Name', sprintf('Pure Schmitt Trigger - ROI %d', demo_roi), ...
    'Position', [100, 100, 1600, 1200]);

% Get thresholds for this ROI
if demo_roi <= length(event_stats_current.upper_thresholds)
    upper_thresh = event_stats_current.upper_thresholds(demo_roi);
    lower_thresh = event_stats_current.lower_thresholds(demo_roi);
else
    upper_thresh = 0.1;
    lower_thresh = 0.05;
end

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

yline(upper_thresh, 'g--', sprintf('Upper (%.4f)', upper_thresh), 'LineWidth', 1.5);
yline(lower_thresh, 'c--', sprintf('Lower (%.4f)', lower_thresh), 'LineWidth', 1.5);
yline(0, 'k:', 'Color', [0.5 0.5 0.5]);

xlabel('Time (s)'); ylabel('dF/F');
title(sprintf('ROI %d: Current (Red, %d events) vs Pure Schmitt (Blue, %d events)', ...
    demo_roi, length(current_starts), length(pure_starts)));
legend('dF/F', 'Current (Extended)', 'Pure Schmitt', 'Upper', 'Lower', 'Location', 'best');
grid on;

% Plot 2: Zoomed view showing exact threshold crossings
subplot(3, 1, 2);
if ~isempty(pure_starts) && length(pure_starts) >= 1
    % Zoom around first pure event
    first_event = pure_starts(1);
    zoom_start = max(1, first_event - 10);
    zoom_end = min(length(roi_trace), first_event + 30);
    zoom_indices = zoom_start:zoom_end;
    
    plot(time_vector(zoom_indices), roi_trace(zoom_indices), 'k-', 'LineWidth', 1.5);
    hold on;
    
    % Current events in zoom
    current_zoom = current_events(zoom_indices);
    if any(current_zoom)
        current_trace_zoom = roi_trace(zoom_indices);
        current_trace_zoom(~current_zoom) = NaN;
        plot(time_vector(zoom_indices), current_trace_zoom, 'r-', 'LineWidth', 3);
    end
    
    % Pure events in zoom
    pure_zoom = pure_events(zoom_indices);
    if any(pure_zoom)
        pure_trace_zoom = roi_trace(zoom_indices);
        pure_trace_zoom(~pure_zoom) = NaN;
        plot(time_vector(zoom_indices), pure_trace_zoom, 'b-', 'LineWidth', 2);
    end
    
    yline(upper_thresh, 'g--', 'LineWidth', 1.5);
    yline(lower_thresh, 'c--', 'LineWidth', 1.5);
    yline(0, 'k:', 'Color', [0.5 0.5 0.5]);
    
    % Mark exact threshold crossings for pure events
    for i = zoom_start:zoom_end-1
        if ~pure_events(i) && pure_events(i+1)  % Event start
            plot(time_vector(i+1), roi_trace(i+1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
        end
        if pure_events(i) && ~pure_events(i+1)  % Event end
            plot(time_vector(i), roi_trace(i), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        end
    end
    
    xlabel('Time (s)'); ylabel('dF/F');
    title('Zoomed View: Pure Schmitt Ends Exactly at Lower Threshold');
    legend('dF/F', 'Current (Extended)', 'Pure Schmitt', 'Upper', 'Lower', 'Zero', ...
        'Event Start', 'Event End', 'Location', 'best');
    grid on;
else
    text(0.5, 0.5, 'No pure events to zoom', 'HorizontalAlignment', 'center');
    title('No Pure Events Detected');
end

% Plot 3: Event duration comparison
subplot(3, 1, 3);

% Get all event durations
current_durations = [];
if ~isempty(current_starts)
    current_durations = current_ends - current_starts + 1;
end

pure_durations = [];
if ~isempty(pure_starts)
    pure_durations = pure_ends - pure_starts + 1;
end

% Create histogram comparison
if ~isempty(current_durations) || ~isempty(pure_durations)
    max_duration = max([current_durations, pure_durations, 20]);  % At least up to 20 frames
    bins = 1:max_duration;
    
    if ~isempty(current_durations)
        h1 = histogram(current_durations, bins, 'FaceAlpha', 0.5, 'FaceColor', 'r', 'EdgeColor', 'none');
        hold on;
    end
    
    if ~isempty(pure_durations)
        h2 = histogram(pure_durations, bins, 'FaceAlpha', 0.5, 'FaceColor', 'b', 'EdgeColor', 'none');
        hold on;
    end
    
    xlabel('Event Duration (frames)'); ylabel('Count');
    title(sprintf('Event Duration Distribution - ROI %d', demo_roi));
    
    if ~isempty(current_durations) && ~isempty(pure_durations)
        legend('Current (Extended)', 'Pure Schmitt', 'Location', 'best');
    elseif ~isempty(current_durations)
        legend('Current (Extended)', 'Location', 'best');
    elseif ~isempty(pure_durations)
        legend('Pure Schmitt', 'Location', 'best');
    end
    
    grid on;
else
    text(0.5, 0.5, 'No events detected', 'HorizontalAlignment', 'center');
    title('No Events for Duration Analysis');
end

sgtitle(sprintf('Pure Schmitt Trigger: Events End at Exact Threshold Crossings'), ...
    'FontSize', 14, 'FontWeight', 'bold');

%% === Dataset-Wide Comparison ===
fprintf('\n=== DATASET-WIDE COMPARISON ===\n');
fprintf('Metric                        Current    Pure Schmitt    Difference\n');
fprintf('----------------------------  ---------  --------------  ----------\n');
fprintf('Total events                  %8d      %8d        %+d\n', ...
    event_stats_current.total_events, event_stats_pure.total_events, ...
    event_stats_pure.total_events - event_stats_current.total_events);

fprintf('Active ROIs                   %8d      %8d        %+d\n', ...
    event_stats_current.rois_with_events, event_stats_pure.rois_with_events, ...
    event_stats_pure.rois_with_events - event_stats_current.rois_with_events);

fprintf('Mean event duration (frames)  %8.1f      %8.1f        %+.1f\n', ...
    event_stats_current.mean_event_duration, event_stats_pure.mean_event_duration, ...
    event_stats_pure.mean_event_duration - event_stats_current.mean_event_duration);

fprintf('Mean event duration (ms)      %8.0f      %8.0f        %+.0f\n', ...
    event_stats_current.mean_event_duration * 1000 / config.frame_rate, ...
    event_stats_pure.mean_event_duration * 1000 / config.frame_rate, ...
    (event_stats_pure.mean_event_duration - event_stats_current.mean_event_duration) * 1000 / config.frame_rate);

%% === Key Insights ===
fprintf('\n=== KEY INSIGHTS ===\n');
fprintf('Pure Schmitt trigger behavior:\n');
fprintf('1. Events START when signal crosses ABOVE upper threshold\n');
fprintf('2. Events END when signal drops BELOW lower threshold\n');
fprintf('3. NO artificial extensions - shows true biological event duration\n');
fprintf('4. Still filters out <3 frame noise spikes\n');
fprintf('5. Still merges events separated by ≤2 frames\n\n');

if event_stats_pure.mean_event_duration < event_stats_current.mean_event_duration
    fprintf('✅ IMPROVEMENT: Pure method shows shorter, more accurate event durations\n');
    fprintf('   Extensions were artificially lengthening events beyond biological reality\n');
else
    fprintf('📊 RESULT: Pure method shows similar durations to current method\n');
    fprintf('   Extensions were not significantly affecting this dataset\n');
end

fprintf('\n=== RECOMMENDATION ===\n');
fprintf('Use the PURE SCHMITT TRIGGER for accurate event timing:\n');
fprintf('1. Events end exactly where signal biology dictates\n');
fprintf('2. No artificial kinetic extensions masking real duration\n');
fprintf('3. Better for downstream analysis requiring precise timing\n');
fprintf('4. More honest representation of glutamate release kinetics\n');

%% === Integration Instructions ===
fprintf('\n=== INTEGRATION INTO YOUR PIPELINE ===\n');
fprintf('Replace your event detector call with:\n');
fprintf('   [event_mask, event_stats] = pure_schmitt_trigger_detector(dfof_data, config, baseline_stats);\n\n');

fprintf('This will give you:\n');
fprintf('• Events that end exactly at lower threshold crossings\n');
fprintf('• No artificial extensions beyond biological reality\n');
fprintf('• More accurate event timing for analysis\n');
fprintf('• Same noise filtering (min 3 frames, merge nearby)\n');

fprintf('\n✅ Pure Schmitt trigger analysis complete!\n');