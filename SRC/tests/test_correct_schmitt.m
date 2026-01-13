%% TEST CORRECT SCHMITT TRIGGER IMPLEMENTATION
% Compare the CORRECT validation-based implementation with previous versions

clear; clc; close all;

fprintf('=== Testing CORRECT Schmitt Trigger Implementation ===\n');

%% === Load Data ===
test_file = 'D:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglusnfr4f_NGR\Spont\GPU_SNR_Processed\5_raw_mean\CP_Snfr4-NGR_Doc2b-R213W_Cs1-c1_spont-01_mean.csv';

try
    loader = csv_loader_v2();
    [data, metadata] = loader.loadSingleFile(test_file);
    config = tracenorm_config();
    config.verbose = true;  % Enable debug output for ROI 613
    
    [baseline, outlier_mask, baseline_stats] = baseline_detector(data, config);
    [dfof_data, dfof_stats] = dfof_calculator(data, baseline, config);
    
    fprintf('✓ Data loaded: %d frames × %d ROIs\n', size(data, 1), size(data, 2));
catch ME
    fprintf('✗ Error loading data: %s\n', ME.message);
    return;
end

%% === Test ROI 613 ===
demo_roi = 613;
roi_trace = dfof_data(:, demo_roi);
time_vector = (1:length(roi_trace)) / config.frame_rate;

fprintf('\n=== Testing ROI %d ===\n', demo_roi);

%% === 1. Current Implementation ===
fprintf('\n1. CURRENT implementation (with extensions):\n');
try
    [events_current, stats_current] = optimal_schmitt_event_detector(dfof_data, config, baseline_stats);
    current_events = events_current(:, demo_roi);
    current_starts = find(diff([false; current_events]) == 1);
    current_ends = find(diff([current_events; false]) == -1);
    
    fprintf('   Current method: %d events\n', length(current_starts));
    for i = 1:min(3, length(current_starts))
        duration = current_ends(i) - current_starts(i) + 1;
        fprintf('     Event %d: [%d-%d] = %d frames\n', i, current_starts(i), current_ends(i), duration);
    end
catch ME
    fprintf('   Error: %s\n', ME.message);
    current_events = false(size(roi_trace));
    current_starts = [];
end

%% === 2. Previous "Fixed" Implementation ===
fprintf('\n2. Previous "fixed" implementation (incorrect):\n');
try
    % Calculate thresholds
    valid_data = roi_trace(~isnan(roi_trace));
    baseline_70th = prctile(valid_data, 70);
    baseline_data = valid_data(valid_data <= baseline_70th);
    noise_std = mad(baseline_data, 1) * 1.4826;
    upper_thresh = config.event_detection.upper_threshold_sigma * noise_std;
    lower_thresh = config.event_detection.lower_threshold_sigma * noise_std;
    
    % Apply simple Schmitt trigger (no validation)
    simple_events = apply_simple_schmitt(roi_trace, upper_thresh, lower_thresh);
    simple_starts = find(diff([false; simple_events]) == 1);
    simple_ends = find(diff([simple_events; false]) == -1);
    
    fprintf('   Simple Schmitt: %d events\n', length(simple_starts));
    for i = 1:min(3, length(simple_starts))
        duration = simple_ends(i) - simple_starts(i) + 1;
        fprintf('     Event %d: [%d-%d] = %d frames\n', i, simple_starts(i), simple_ends(i), duration);
    end
catch ME
    fprintf('   Error: %s\n', ME.message);
    simple_events = false(size(roi_trace));
    simple_starts = [];
end

%% === 3. CORRECT Implementation (with validation) ===
fprintf('\n3. CORRECT implementation (with validation period):\n');
try
    [events_correct, stats_correct] = correct_schmitt_trigger_detector(dfof_data, config, baseline_stats);
    correct_events = events_correct(:, demo_roi);
    correct_starts = find(diff([false; correct_events]) == 1);
    correct_ends = find(diff([correct_events; false]) == -1);
    
    fprintf('\n   Correct method: %d events\n', length(correct_starts));
    for i = 1:min(3, length(correct_starts))
        duration = correct_ends(i) - correct_starts(i) + 1;
        fprintf('     Event %d: [%d-%d] = %d frames\n', i, correct_starts(i), correct_ends(i), duration);
    end
catch ME
    fprintf('   Error: %s\n', ME.message);
    correct_events = false(size(roi_trace));
    correct_starts = [];
end

%% === Visual Comparison ===
fprintf('\n4. Creating visual comparison...\n');

% Calculate thresholds for plotting
valid_data = roi_trace(~isnan(roi_trace));
baseline_70th = prctile(valid_data, 70);
baseline_data = valid_data(valid_data <= baseline_70th);
noise_std = mad(baseline_data, 1) * 1.4826;
upper_threshold = config.event_detection.upper_threshold_sigma * noise_std;
lower_threshold = config.event_detection.lower_threshold_sigma * noise_std;

fig = figure('Name', sprintf('CORRECT Schmitt Trigger - ROI %d', demo_roi), ...
    'Position', [100, 100, 1600, 1200]);

% Plot 1: Full trace comparison
subplot(3, 1, 1);
plot(time_vector, roi_trace, 'k-', 'LineWidth', 1.0);
hold on;

% Highlight events from different methods
if any(current_events)
    current_trace = roi_trace;
    current_trace(~current_events) = NaN;
    plot(time_vector, current_trace, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Current');
end

if any(simple_events)
    simple_trace = roi_trace;
    simple_trace(~simple_events) = NaN;
    plot(time_vector, simple_trace, 'b-', 'LineWidth', 2.0, 'DisplayName', 'Simple Schmitt');
end

if any(correct_events)
    correct_trace = roi_trace;
    correct_trace(~correct_events) = NaN;
    plot(time_vector, correct_trace, 'm-', 'LineWidth', 3.0, 'DisplayName', 'CORRECT Schmitt');
end

% Add threshold lines
line([time_vector(1), time_vector(end)], [upper_threshold, upper_threshold], ...
    'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Upper');
line([time_vector(1), time_vector(end)], [lower_threshold, lower_threshold], ...
    'Color', 'c', 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Lower');
line([time_vector(1), time_vector(end)], [0, 0], ...
    'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 1);

xlabel('Time (s)'); ylabel('dF/F');
title(sprintf('ROI %d: Current (%d) vs Simple (%d) vs CORRECT (%d) events', ...
    demo_roi, length(current_starts), length(simple_starts), length(correct_starts)));
legend('dF/F', 'Current', 'Simple', 'CORRECT', 'Upper', 'Lower', 'Location', 'best');
grid on;

% Plot 2: Zoomed view of first correct event
subplot(3, 1, 2);
if ~isempty(correct_starts)
    first_event = correct_starts(1);
    zoom_start = max(1, first_event - 15);
    zoom_end = min(length(roi_trace), first_event + 30);
    zoom_indices = zoom_start:zoom_end;
    
    plot(time_vector(zoom_indices), roi_trace(zoom_indices), 'k-', 'LineWidth', 2);
    hold on;
    
    % Highlight validated event
    correct_zoom = correct_events(zoom_indices);
    if any(correct_zoom)
        correct_trace_zoom = roi_trace(zoom_indices);
        correct_trace_zoom(~correct_zoom) = NaN;
        plot(time_vector(zoom_indices), correct_trace_zoom, 'm-', 'LineWidth', 4);
    end
    
    % Add threshold lines
    line([time_vector(zoom_start), time_vector(zoom_end)], [upper_threshold, upper_threshold], ...
        'Color', 'g', 'LineStyle', '--', 'LineWidth', 2);
    line([time_vector(zoom_start), time_vector(zoom_end)], [lower_threshold, lower_threshold], ...
        'Color', 'c', 'LineStyle', '--', 'LineWidth', 2);
    
    % Mark validation period
    if length(correct_starts) > 0
        validation_end = min(zoom_end, correct_starts(1) + 2);  % 3 frames = indices 0,1,2
        validation_time = (correct_starts(1):validation_end) / config.frame_rate;
        plot(validation_time, roi_trace(correct_starts(1):validation_end), 'go', ...
            'MarkerSize', 8, 'MarkerFaceColor', 'g', 'DisplayName', 'Validation Period');
    end
    
    xlabel('Time (s)'); ylabel('dF/F');
    title('Zoomed: CORRECT Method Shows Validation Period (Green Circles)');
    legend('dF/F', 'CORRECT Event', 'Upper', 'Lower', 'Validation', 'Location', 'best');
    grid on;
else
    text(0.5, 0.5, 'No correct events to zoom', 'HorizontalAlignment', 'center');
    title('No CORRECT Events Detected');
end

% Plot 3: Event duration comparison
subplot(3, 1, 3);

% Get durations
current_durations = [];
if ~isempty(current_starts)
    current_durations = current_ends - current_starts + 1;
end

simple_durations = [];
if ~isempty(simple_starts)
    simple_durations = simple_ends - simple_starts + 1;
end

correct_durations = [];
if ~isempty(correct_starts)
    correct_durations = correct_ends - correct_starts + 1;
end

% Plot histogram
max_dur = max([current_durations, simple_durations, correct_durations, 20]);
bins = 1:max_dur;

if ~isempty(current_durations)
    [counts1, ~] = hist(current_durations, bins);
    bar(bins, counts1, 'FaceAlpha', 0.3, 'FaceColor', 'r');
    hold on;
end

if ~isempty(simple_durations)
    [counts2, ~] = hist(simple_durations, bins);
    bar(bins, counts2, 'FaceAlpha', 0.3, 'FaceColor', 'b');
end

if ~isempty(correct_durations)
    [counts3, ~] = hist(correct_durations, bins);
    bar(bins, counts3, 'FaceAlpha', 0.7, 'FaceColor', 'm');
end

xlabel('Event Duration (frames)'); ylabel('Count');
title('Event Duration Distribution');
legend('Current', 'Simple', 'CORRECT', 'Location', 'best');
grid on;

%% === Summary Comparison ===
fprintf('\n=== SUMMARY COMPARISON ===\n');
fprintf('Method                     Events    Avg Duration    Notes\n');
fprintf('-------------------------  --------  --------------  ------------------------\n');

if ~isempty(current_starts)
    current_avg = mean(current_ends - current_starts + 1);
else
    current_avg = 0;
end

if ~isempty(simple_starts)
    simple_avg = mean(simple_ends - simple_starts + 1);
else
    simple_avg = 0;
end

if ~isempty(correct_starts)
    correct_avg = mean(correct_ends - correct_starts + 1);
else
    correct_avg = 0;
end

fprintf('Current (extensions)       %8d      %8.1f        With kinetic extensions\n', ...
    length(current_starts), current_avg);
fprintf('Simple Schmitt             %8d      %8.1f        No validation period\n', ...
    length(simple_starts), simple_avg);
fprintf('CORRECT (validation)       %8d      %8.1f        ≥3 frame validation\n', ...
    length(correct_starts), correct_avg);

fprintf('\n=== KEY INSIGHT ===\n');
fprintf('The CORRECT method requires ≥3 frames above lower threshold\n');
fprintf('before accepting an event. This eliminates noise spikes that\n');
fprintf('briefly cross thresholds but aren''t sustained biological events.\n');

if length(correct_starts) < length(simple_starts)
    fprintf('\n✅ CORRECT method detected FEWER events → better noise rejection\n');
else
    fprintf('\n📊 CORRECT method shows similar detection to simple method\n');
end

fprintf('\n✅ CORRECT Schmitt trigger testing complete!\n');

%% === Helper Function ===
function event_mask = apply_simple_schmitt(trace, upper_thresh, lower_thresh)
    % Simple Schmitt trigger without validation (for comparison)
    
    event_mask = false(size(trace));
    state = 'baseline';
    event_start = 0;
    events = [];
    
    for frame = 1:length(trace)
        signal = trace(frame);
        if isnan(signal), continue; end
        
        if strcmp(state, 'baseline') && signal > upper_thresh
            state = 'in_event';
            event_start = frame;
        elseif strcmp(state, 'in_event') && signal < lower_thresh
            events(end+1, :) = [event_start, frame-1];
            state = 'baseline';
        end
    end
    
    if strcmp(state, 'in_event')
        events(end+1, :) = [event_start, length(trace)];
    end
    
    % Apply minimum duration filter
    for i = 1:size(events, 1)
        if events(i, 2) - events(i, 1) + 1 >= 3
            event_mask(events(i, 1):events(i, 2)) = true;
        end
    end
end
