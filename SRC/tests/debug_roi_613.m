%% DEBUG ROI 613 SCHMITT TRIGGER
% Detailed analysis of why events are detected differently
% FIXED: MATLAB compatibility issues

clear; clc; close all;

fprintf('=== Debug ROI 613 Schmitt Trigger ===\n');

%% === Load and Prepare Data ===
test_file = 'D:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglusnfr4f_NGR\Spont\GPU_SNR_Processed\5_raw_mean\CP_Snfr4-NGR_Doc2b-R213W_Cs1-c1_spont-01_mean.csv';

try
    loader = csv_loader_v2();
    [data, metadata] = loader.loadSingleFile(test_file);
    config = tracenorm_config();
    [baseline, outlier_mask, baseline_stats] = baseline_detector(data, config);
    [dfof_data, dfof_stats] = dfof_calculator(data, baseline, config);
    
    fprintf('✓ Data loaded: %d frames × %d ROIs\n', size(data, 1), size(data, 2));
catch ME
    fprintf('✗ Error loading data: %s\n', ME.message);
    return;
end

%% === Focus on ROI 613 ===
demo_roi = 613;
roi_trace = dfof_data(:, demo_roi);
time_vector = (1:length(roi_trace)) / config.frame_rate;

fprintf('\n=== ROI %d Analysis ===\n', demo_roi);
fprintf('Trace range: [%.6f, %.6f]\n', min(roi_trace), max(roi_trace));
fprintf('Valid frames: %d/%d\n', sum(~isnan(roi_trace)), length(roi_trace));

%% === Calculate Thresholds (FIXED MATLAB compatibility) ===
fprintf('\n--- Threshold Calculation Methods ---\n');

% Remove NaN values first
valid_data = roi_trace(~isnan(roi_trace));

% Method 1: Current implementation - use full trace
full_trace_mad = mad(valid_data, 1);  % Remove 'omitnan' flag
full_trace_std = full_trace_mad * 1.4826;
upper_1 = config.event_detection.upper_threshold_sigma * full_trace_std;
lower_1 = config.event_detection.lower_threshold_sigma * full_trace_std;

fprintf('Method 1 (full trace MAD): σ=%.6f, upper=%.6f, lower=%.6f\n', ...
    full_trace_std, upper_1, lower_1);

% Method 2: Research paper approach - exclude likely events
baseline_70th = prctile(valid_data, 70);
baseline_data = valid_data(valid_data <= baseline_70th);
if length(baseline_data) > 10
    baseline_mad = mad(baseline_data, 1);
    baseline_std = baseline_mad * 1.4826;
else
    baseline_std = full_trace_std;  % Fallback
end
upper_2 = config.event_detection.upper_threshold_sigma * baseline_std;
lower_2 = config.event_detection.lower_threshold_sigma * baseline_std;

fprintf('Method 2 (70th percentile): σ=%.6f, upper=%.6f, lower=%.6f\n', ...
    baseline_std, upper_2, lower_2);

% Show what fraction of data is above each threshold
fprintf('\nFraction of frames above thresholds:\n');
fprintf('  Method 1 upper: %.1f%% | lower: %.1f%%\n', ...
    100*sum(valid_data > upper_1)/length(valid_data), ...
    100*sum(valid_data > lower_1)/length(valid_data));
fprintf('  Method 2 upper: %.1f%% | lower: %.1f%%\n', ...
    100*sum(valid_data > upper_2)/length(valid_data), ...
    100*sum(valid_data > lower_2)/length(valid_data));

%% === Manual Schmitt Trigger Implementation ===
fprintf('\n--- Manual Schmitt Trigger with Method 2 Thresholds ---\n');

% Use method 2 thresholds (research paper approach)
upper_thresh = upper_2;
lower_thresh = lower_2;

events_manual = false(size(roi_trace));
state = 'baseline';
event_start = 0;
event_count = 0;

for frame = 1:length(roi_trace)
    signal = roi_trace(frame);
    
    if isnan(signal)
        continue;
    end
    
    if strcmp(state, 'baseline')
        if signal > upper_thresh
            state = 'in_event';
            event_start = frame;
            event_count = event_count + 1;
            fprintf('Event %d START at frame %d (%.3fs): signal=%.6f > %.6f\n', ...
                event_count, frame, frame/config.frame_rate, signal, upper_thresh);
        end
    elseif strcmp(state, 'in_event')
        if signal < lower_thresh
            % Event ends
            event_end = frame - 1;  % Last frame above lower threshold
            events_manual(event_start:event_end) = true;
            
            duration = event_end - event_start + 1;
            fprintf('Event %d END at frame %d (%.3fs): signal=%.6f < %.6f, duration=%d frames\n', ...
                event_count, event_end, event_end/config.frame_rate, signal, lower_thresh, duration);
            
            state = 'baseline';
            event_start = 0;
        end
    end
end

% Handle event extending to end
if strcmp(state, 'in_event') && event_start > 0
    events_manual(event_start:end) = true;
    duration = length(roi_trace) - event_start + 1;
    fprintf('Event %d extends to END: duration=%d frames\n', event_count, duration);
end

%% === Apply Minimum Duration Filter ===
fprintf('\n--- Applying 3-frame minimum duration filter ---\n');

events_filtered = false(size(events_manual));
event_starts = find(diff([false; events_manual]) == 1);
event_ends = find(diff([events_manual; false]) == -1);

fprintf('Before filtering: %d raw events\n', length(event_starts));

for i = 1:length(event_starts)
    start_frame = event_starts(i);
    end_frame = event_ends(i);
    duration = end_frame - start_frame + 1;
    
    if duration >= 3  % Minimum duration filter
        events_filtered(start_frame:end_frame) = true;
        fprintf('  Keeping event %d: [%d-%d] duration=%d\n', i, start_frame, end_frame, duration);
    else
        fprintf('  Filtering event %d: [%d-%d] duration=%d (too short)\n', i, start_frame, end_frame, duration);
    end
end

filtered_starts = find(diff([false; events_filtered]) == 1);
filtered_ends = find(diff([events_filtered; false]) == -1);

fprintf('After filtering: %d events kept\n', length(filtered_starts));

%% === Visual Comparison ===
figure('Name', 'Debug ROI 613 Schmitt Trigger', 'Position', [100, 100, 1400, 800]);

% Plot the trace with events and thresholds
plot(time_vector, roi_trace, 'k-', 'LineWidth', 1);
hold on;

% Highlight detected events
if any(events_filtered)
    event_trace = roi_trace;
    event_trace(~events_filtered) = NaN;
    plot(time_vector, event_trace, 'r-', 'LineWidth', 3);
end

% Show thresholds
yline(upper_thresh, 'g--', sprintf('Upper (%.6f)', upper_thresh), 'LineWidth', 2);
yline(lower_thresh, 'c--', sprintf('Lower (%.6f)', lower_thresh), 'LineWidth', 2);
yline(0, 'k:', 'LineWidth', 1);

% Mark event start/end points
for i = 1:length(filtered_starts)
    start_time = filtered_starts(i) / config.frame_rate;
    end_time = filtered_ends(i) / config.frame_rate;
    
    plot(start_time, roi_trace(filtered_starts(i)), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
    plot(end_time, roi_trace(filtered_ends(i)), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
end

xlabel('Time (s)');
ylabel('dF/F');
title(sprintf('ROI %d: Manual Schmitt Trigger Debug (Method 2 Thresholds)', demo_roi));
legend('dF/F', 'Detected Events', 'Upper Threshold', 'Lower Threshold', 'Zero', ...
    'Event Start', 'Event End', 'Location', 'best');
grid on;

%% === Compare with Original Implementations ===
fprintf('\n--- Comparing with Original Implementations ---\n');

try
    [events_current, stats_current] = optimal_schmitt_event_detector(dfof_data, config, baseline_stats);
    current_events_613 = events_current(:, demo_roi);
    current_starts = find(diff([false; current_events_613]) == 1);
    current_ends = find(diff([current_events_613; false]) == -1);
    
    fprintf('Current implementation: %d events\n', length(current_starts));
    for i = 1:min(3, length(current_starts))
        duration = current_ends(i) - current_starts(i) + 1;
        fprintf('  Event %d: [%d-%d] duration=%d\n', i, current_starts(i), current_ends(i), duration);
    end
catch ME
    fprintf('Error with current implementation: %s\n', ME.message);
end

try
    [events_pure, stats_pure] = pure_schmitt_trigger_detector(dfof_data, config, baseline_stats);
    pure_events_613 = events_pure(:, demo_roi);
    pure_starts = find(diff([false; pure_events_613]) == 1);
    pure_ends = find(diff([pure_events_613; false]) == -1);
    
    fprintf('Original pure implementation: %d events\n', length(pure_starts));
    for i = 1:min(3, length(pure_starts))
        duration = pure_ends(i) - pure_starts(i) + 1;
        fprintf('  Event %d: [%d-%d] duration=%d\n', i, pure_starts(i), pure_ends(i), duration);
    end
catch ME
    fprintf('Error with pure implementation: %s\n', ME.message);
end

fprintf('Manual implementation: %d events\n', length(filtered_starts));
for i = 1:min(3, length(filtered_starts))
    duration = filtered_ends(i) - filtered_starts(i) + 1;
    fprintf('  Event %d: [%d-%d] duration=%d\n', i, filtered_starts(i), filtered_ends(i), duration);
end

%% === Detailed Trace Analysis ===
fprintf('\n--- Detailed Trace Analysis ---\n');

% Look for periods where signal is above upper threshold
above_upper = valid_data > upper_thresh;
above_lower = valid_data > lower_thresh;

fprintf('Frames above upper threshold: %d (%.1f%%)\n', sum(above_upper), 100*sum(above_upper)/length(valid_data));
fprintf('Frames above lower threshold: %d (%.1f%%)\n', sum(above_lower), 100*sum(above_lower)/length(valid_data));

% Find all potential event start points
potential_starts = find(diff([false; above_upper]) == 1);
fprintf('Potential event starts (crossing above upper): %d\n', length(potential_starts));

if length(potential_starts) > 0
    fprintf('First few potential starts at frames: ');
    fprintf('%d ', potential_starts(1:min(5, length(potential_starts))));
    fprintf('\n');
end

%% === Summary ===
fprintf('\n=== SUMMARY ===\n');
fprintf('Key findings for ROI %d:\n', demo_roi);
fprintf('• Method 1 (full trace): σ=%.6f → thresholds: [%.6f, %.6f]\n', full_trace_std, upper_1, lower_1);
fprintf('• Method 2 (baseline only): σ=%.6f → thresholds: [%.6f, %.6f]\n', baseline_std, upper_2, lower_2);
fprintf('• Manual Schmitt trigger found %d events\n', length(filtered_starts));
fprintf('• %d potential starts found by crossing upper threshold\n', length(potential_starts));

if baseline_std < full_trace_std
    fprintf('\nBaseline-only noise estimation gives LOWER thresholds → MORE sensitive detection\n');
    fprintf('This follows the research paper approach of excluding events from noise calculation.\n');
else
    fprintf('\nBaseline-only noise estimation gives HIGHER thresholds → LESS sensitive detection\n');
end

fprintf('\n✅ Debug complete!\n');
