%% CLARIFIED THRESHOLD PARAMETERS
% Separate noise calculation exclusion from event detection minimum duration

clear; clc; close all;

fprintf('=== CLARIFIED: Two Different Window Parameters ===\n\n');

%% === Load Data ===
test_file = 'D:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglusnfr4f_NGR\Spont\GPU_SNR_Processed\5_raw_mean\CP_Snfr4-NGR_Doc2b-R213W_Cs1-c1_spont-01_mean.csv';

try
    loader = csv_loader_v2();
    [data, metadata] = loader.loadSingleFile(test_file);
    config = tracenorm_config();
    
    [baseline, outlier_mask, baseline_stats] = baseline_detector(data, config);
    [dfof_data, dfof_stats] = dfof_calculator(data, baseline, config);
    
    fprintf('✓ Data loaded\n\n');
catch ME
    fprintf('✗ Error loading data: %s\n', ME.message);
    return;
end

%% === Two Separate Parameters ===
fprintf('=== TWO DIFFERENT PARAMETERS ===\n');
fprintf('1. NOISE EXCLUSION WINDOW: How long must a high period be to exclude from noise calculation?\n');
fprintf('2. EVENT MINIMUM DURATION: How long must an event be to count as real?\n\n');

% Test ROI
test_roi = 744;
roi_trace = dfof_data(:, test_roi);
roi_outliers = outlier_mask(:, test_roi);

%% === Parameter 1: Noise Exclusion Window ===
fprintf('PARAMETER 1: Noise Exclusion Window\n');
fprintf('Purpose: Exclude sustained biological events from noise calculation\n');
fprintf('Recommendation: Use longer window (7-10 frames) to exclude only clear events\n\n');

noise_exclusion_windows = [5, 7, 10];  % Test different exclusion windows
upper_sigma = 3.5;
lower_sigma = 1.5;

noise_results = struct();

for w = 1:length(noise_exclusion_windows)
    exclusion_window = noise_exclusion_windows(w);
    
    fprintf('Testing %d-frame noise exclusion window:\n', exclusion_window);
    
    % Calculate noise excluding sustained periods
    [noise_std, excluded_frames] = calculate_noise_with_exclusion(roi_trace, roi_outliers, exclusion_window);
    
    % Calculate thresholds
    upper_thresh = upper_sigma * noise_std;
    lower_thresh = lower_sigma * noise_std;
    
    fprintf('  • Excluded %d frames as "sustained events"\n', excluded_frames);
    fprintf('  • Noise std: %.6f\n', noise_std);
    fprintf('  • Upper threshold (%.1fσ): %.6f\n', upper_sigma, upper_thresh);
    fprintf('  • Lower threshold (%.1fσ): %.6f\n\n', lower_sigma, lower_thresh);
    
    noise_results(w).exclusion_window = exclusion_window;
    noise_results(w).noise_std = noise_std;
    noise_results(w).excluded_frames = excluded_frames;
    noise_results(w).upper_thresh = upper_thresh;
    noise_results(w).lower_thresh = lower_thresh;
end

%% === Parameter 2: Event Minimum Duration ===
fprintf('PARAMETER 2: Event Minimum Duration\n');
fprintf('Purpose: Filter out brief spikes in Schmitt trigger detection\n');
fprintf('Recommendation: Use 3-4 frames for glutamate release kinetics\n\n');

% Use the middle noise exclusion result (7-frame)
best_noise_idx = 2;  % 7-frame exclusion
upper_thresh = noise_results(best_noise_idx).upper_thresh;
lower_thresh = noise_results(best_noise_idx).lower_thresh;

event_min_durations = [3, 4, 5];  % Test different minimum event durations
event_results = struct();

for d = 1:length(event_min_durations)
    min_duration = event_min_durations(d);
    
    fprintf('Testing %d-frame minimum event duration:\n', min_duration);
    
    % Test event detection with this minimum duration
    events = test_schmitt_with_duration(roi_trace, upper_thresh, lower_thresh, min_duration);
    event_starts = find(diff([false; events]) == 1);
    event_ends = find(diff([events; false]) == -1);
    
    if ~isempty(event_starts)
        durations = event_ends - event_starts + 1;
        avg_duration = mean(durations);
        total_frames = sum(durations);
    else
        avg_duration = 0;
        total_frames = 0;
    end
    
    fprintf('  • Detected events: %d\n', length(event_starts));
    fprintf('  • Average duration: %.1f frames\n', avg_duration);
    fprintf('  • Total event frames: %d\n\n', total_frames);
    
    event_results(d).min_duration = min_duration;
    event_results(d).events = events;
    event_results(d).count = length(event_starts);
    event_results(d).avg_duration = avg_duration;
    event_results(d).total_frames = total_frames;
end

%% === VISUALIZATION ===
fig = figure('Name', sprintf('Parameter Clarification - ROI %d', test_roi), ...
    'Position', [100, 100, 1600, 1200]);

% Plot 1: Noise exclusion window effects
subplot(2, 3, 1);
exclusion_windows = [noise_results.exclusion_window];
noise_stds = [noise_results.noise_std];
excluded_counts = [noise_results.excluded_frames];

yyaxis left;
bar(exclusion_windows, noise_stds, 'FaceColor', [0.7 0.3 0.3], 'EdgeColor', 'k');
ylabel('Noise Standard Deviation');
xlabel('Noise Exclusion Window (frames)');

yyaxis right;
plot(exclusion_windows, excluded_counts, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
ylabel('Excluded Frames');

title('Noise Exclusion Window Effect');
grid on;

% Plot 2: Event minimum duration effects
subplot(2, 3, 2);
min_durations = [event_results.min_duration];
event_counts = [event_results.count];
avg_durations = [event_results.avg_duration];

yyaxis left;
bar(min_durations, event_counts, 'FaceColor', [0.3 0.7 0.3], 'EdgeColor', 'k');
ylabel('Event Count');
xlabel('Minimum Event Duration (frames)');

yyaxis right;
plot(min_durations, avg_durations, 'ro-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
ylabel('Average Event Duration');

title('Event Duration Filter Effect');
grid on;

% Plot 3: Trace with recommended parameters
subplot(2, 3, 3);
time_vector = (1:length(roi_trace)) / config.frame_rate;
plot(time_vector, roi_trace, 'k-', 'LineWidth', 1);
hold on;

% Use recommended parameters: 7-frame exclusion, 3-frame minimum
recommended_events = event_results(1).events;  % 3-frame minimum
if any(recommended_events)
    event_trace = roi_trace;
    event_trace(~recommended_events) = NaN;
    plot(time_vector, event_trace, 'b-', 'LineWidth', 2);
end

line([time_vector(1), time_vector(end)], [upper_thresh, upper_thresh], ...
    'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
line([time_vector(1), time_vector(end)], [lower_thresh, lower_thresh], ...
    'Color', 'c', 'LineStyle', '--', 'LineWidth', 2);

title(sprintf('Recommended: 7-frame exclusion, 3-frame minimum (%d events)', event_results(1).count));
xlabel('Time (s)'); ylabel('dF/F');
legend('dF/F', 'Events', 'Upper', 'Lower', 'Location', 'best');
grid on;

% Plot 4-6: Detailed comparison of event detection with different parameters
for d = 1:3
    subplot(2, 3, 3+d);
    plot(time_vector, roi_trace, 'k-', 'LineWidth', 1);
    hold on;
    
    events = event_results(d).events;
    if any(events)
        event_trace = roi_trace;
        event_trace(~events) = NaN;
        plot(time_vector, event_trace, 'b-', 'LineWidth', 2);
    end
    
    line([time_vector(1), time_vector(end)], [upper_thresh, upper_thresh], ...
        'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
    line([time_vector(1), time_vector(end)], [lower_thresh, lower_thresh], ...
        'Color', 'c', 'LineStyle', '--', 'LineWidth', 1.5);
    
    title(sprintf('%d-frame minimum: %d events', event_min_durations(d), event_results(d).count));
    xlabel('Time (s)'); ylabel('dF/F');
    grid on;
end

sgtitle(sprintf('Parameter Clarification - ROI %d', test_roi), 'FontSize', 14, 'FontWeight', 'bold');

%% === SUMMARY TABLE ===
fprintf('=== PARAMETER SUMMARY ===\n');
fprintf('\nNOISE EXCLUSION WINDOW EFFECTS:\n');
fprintf('Window    Noise Std    Excluded    Upper Thresh    Lower Thresh\n');
fprintf('------    ---------    --------    ------------    ------------\n');
for w = 1:length(noise_results)
    fprintf('%6d    %9.6f    %8d    %12.6f    %12.6f\n', ...
        noise_results(w).exclusion_window, noise_results(w).noise_std, ...
        noise_results(w).excluded_frames, noise_results(w).upper_thresh, ...
        noise_results(w).lower_thresh);
end

fprintf('\nEVENT MINIMUM DURATION EFFECTS (using 7-frame exclusion):\n');
fprintf('Min Dur    Events    Avg Duration    Total Frames\n');
fprintf('-------    ------    ------------    ------------\n');
for d = 1:length(event_results)
    fprintf('%7d    %6d    %12.1f    %12d\n', ...
        event_results(d).min_duration, event_results(d).count, ...
        event_results(d).avg_duration, event_results(d).total_frames);
end

%% === RECOMMENDATIONS ===
fprintf('\n=== RECOMMENDATIONS ===\n');
fprintf('1. NOISE EXCLUSION WINDOW:\n');
fprintf('   • Use 7-10 frames to exclude only clear sustained events\n');
fprintf('   • Longer window = more conservative noise estimate\n');
fprintf('   • Prevents real events from contaminating noise calculation\n\n');

fprintf('2. EVENT MINIMUM DURATION:\n');
fprintf('   • Use 3-4 frames for glutamate release kinetics\n');
fprintf('   • 3 frames = 75ms at 40Hz (biological minimum)\n');
fprintf('   • 4 frames = 100ms (more conservative)\n\n');

fprintf('3. RECOMMENDED COMBINATION:\n');
fprintf('   • 7-frame exclusion window\n');
fprintf('   • 3-frame minimum event duration\n');
fprintf('   • This gives %d events for ROI %d\n', event_results(1).count, test_roi);

fprintf('\n✅ Parameter clarification complete!\n');

%% === Helper Functions ===
function [noise_std, excluded_frames] = calculate_noise_with_exclusion(trace, outliers, exclusion_window)
    % Calculate noise excluding sustained high periods
    
    valid_frames = ~isnan(trace);
    valid_data = trace(valid_frames);
    
    % Find sustained high periods to exclude
    sustained_threshold = prctile(valid_data, 85);
    high_values = trace > sustained_threshold;
    
    % Find runs of consecutive high values
    run_starts = find(diff([false; high_values]) == 1);
    run_ends = find(diff([high_values; false]) == -1);
    
    sustained_mask = false(size(trace));
    for i = 1:length(run_starts)
        run_length = run_ends(i) - run_starts(i) + 1;
        if run_length >= exclusion_window
            sustained_mask(run_starts(i):run_ends(i)) = true;
        end
    end
    
    % Calculate noise from everything except sustained events
    noise_calculation_mask = valid_frames & ~sustained_mask;
    noise_data = trace(noise_calculation_mask);
    
    noise_std = mad(noise_data, 1) * 1.4826;
    excluded_frames = sum(sustained_mask);
end

function events = test_schmitt_with_duration(trace, upper_thresh, lower_thresh, min_duration)
    % Schmitt trigger with specific minimum duration
    
    events = false(size(trace));
    state = 'baseline';
    event_start = 0;
    validated_events = [];
    
    for frame = 1:length(trace)
        signal = trace(frame);
        if isnan(signal), continue; end
        
        if strcmp(state, 'baseline')
            if signal > upper_thresh
                state = 'validating';
                event_start = frame;
            end
        elseif strcmp(state, 'validating')
            frames_since_start = frame - event_start + 1;
            if signal < lower_thresh
                state = 'baseline';
                event_start = 0;
            elseif frames_since_start >= min_duration
                state = 'validated_event';
            end
        elseif strcmp(state, 'validated_event')
            if signal < lower_thresh
                validated_events(end+1, :) = [event_start, frame-1];
                state = 'baseline';
                event_start = 0;
            end
        end
    end
    
    if strcmp(state, 'validated_event') && event_start > 0
        validated_events(end+1, :) = [event_start, length(trace)];
    end
    
    for i = 1:size(validated_events, 1)
        events(validated_events(i, 1):validated_events(i, 2)) = true;
    end
end
