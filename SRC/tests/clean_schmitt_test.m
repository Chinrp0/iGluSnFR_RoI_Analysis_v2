%% CLEAN SCHMITT TRIGGER PARAMETER TESTING
% Focus on parameter refinement with clean visualization

clear; clc; close all;

fprintf('=== Schmitt Trigger Parameter Testing ===\n');

%% === Load Data ===
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

%% === Test ROI 627 ===
test_roi = 627;
roi_trace = dfof_data(:, test_roi);
time_vector = (1:length(roi_trace)) / config.frame_rate;

fprintf('\n=== Testing ROI %d ===\n', test_roi);

%% === Calculate Thresholds ===
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

fprintf('Threshold Analysis:\n');
fprintf('  Noise std: %.6f\n', noise_std);
fprintf('  Upper threshold (%.1fσ): %.6f\n', config.event_detection.upper_threshold_sigma, upper_threshold);
fprintf('  Lower threshold (%.1fσ): %.6f\n', config.event_detection.lower_threshold_sigma, lower_threshold);
fprintf('  Trace range: [%.6f, %.6f]\n', min(valid_data), max(valid_data));
fprintf('  Frames above upper: %d (%.1f%%)\n', sum(valid_data > upper_threshold), 100*sum(valid_data > upper_threshold)/length(valid_data));
fprintf('  Frames above lower: %d (%.1f%%)\n', sum(valid_data > lower_threshold), 100*sum(valid_data > lower_threshold)/length(valid_data));

%% === Test Current Parameters ===
fprintf('\n--- Testing Current Parameters (3.0σ/1.5σ) ---\n');

% Enable debug for this ROI
config_debug = config;
config_debug.verbose = true;

% Temporarily modify the detector to debug our ROI
events_current = test_schmitt_trigger_single_roi(roi_trace, upper_threshold, lower_threshold, test_roi);

current_starts = find(diff([false; events_current]) == 1);
current_ends = find(diff([events_current; false]) == -1);

fprintf('Current parameters detected %d events:\n', length(current_starts));
for i = 1:length(current_starts)
    duration = current_ends(i) - current_starts(i) + 1;
    fprintf('  Event %d: [%d-%d] = %d frames (%.0f ms)\n', ...
        i, current_starts(i), current_ends(i), duration, duration * 1000 / config.frame_rate);
end

%% === Test Alternative Parameters ===
fprintf('\n--- Testing Alternative Parameters ---\n');

% Test different parameter combinations
param_sets = [
    2.5, 1.0;  % More sensitive
    2.5, 1.5;  % More sensitive start, same end
    3.0, 1.0;  % Same start, more sensitive end
    3.5, 1.5;  % Less sensitive start
    3.0, 2.0;  % Same start, less sensitive end
];

param_names = {
    '2.5σ/1.0σ (more sensitive)';
    '2.5σ/1.5σ (sensitive start)';
    '3.0σ/1.0σ (sensitive end)';
    '3.5σ/1.5σ (less sensitive start)';
    '3.0σ/2.0σ (less sensitive end)';
};

alternative_results = cell(size(param_sets, 1), 1);

for p = 1:size(param_sets, 1)
    upper_sigma = param_sets(p, 1);
    lower_sigma = param_sets(p, 2);
    
    alt_upper = upper_sigma * noise_std;
    alt_lower = lower_sigma * noise_std;
    
    fprintf('\nTesting %s:\n', param_names{p});
    fprintf('  Thresholds: %.6f / %.6f\n', alt_upper, alt_lower);
    
    alt_events = test_schmitt_trigger_single_roi(roi_trace, alt_upper, alt_lower, 0);  % No debug
    alt_starts = find(diff([false; alt_events]) == 1);
    alt_ends = find(diff([alt_events; false]) == -1);
    
    fprintf('  Detected %d events\n', length(alt_starts));
    
    alternative_results{p} = struct();
    alternative_results{p}.events = alt_events;
    alternative_results{p}.starts = alt_starts;
    alternative_results{p}.ends = alt_ends;
    alternative_results{p}.upper_thresh = alt_upper;
    alternative_results{p}.lower_thresh = alt_lower;
    alternative_results{p}.name = param_names{p};
end

%% === Visualization ===
fprintf('\n--- Creating Parameter Comparison Visualization ---\n');

fig = figure('Name', sprintf('Schmitt Trigger Parameters - ROI %d', test_roi), ...
    'Position', [100, 100, 1800, 1000]);

% Plot 1: Full trace with current parameters
subplot(2, 1, 1);
plot(time_vector, roi_trace, 'k-', 'LineWidth', 1.0);
hold on;

% Highlight current events
if any(events_current)
    current_trace = roi_trace;
    current_trace(~events_current) = NaN;
    plot(time_vector, current_trace, 'r-', 'LineWidth', 3.0);
end

% Add threshold lines
line([time_vector(1), time_vector(end)], [upper_threshold, upper_threshold], ...
    'Color', 'g', 'LineStyle', '--', 'LineWidth', 2);
line([time_vector(1), time_vector(end)], [lower_threshold, lower_threshold], ...
    'Color', 'c', 'LineStyle', '--', 'LineWidth', 2);
line([time_vector(1), time_vector(end)], [0, 0], ...
    'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 1);

xlabel('Time (s)'); ylabel('dF/F');
title(sprintf('ROI %d: Current Parameters (%.1fσ/%.1fσ) - %d events detected', ...
    test_roi, config.event_detection.upper_threshold_sigma, ...
    config.event_detection.lower_threshold_sigma, length(current_starts)));
legend('dF/F', 'Detected Events', 'Upper Threshold', 'Lower Threshold', 'Zero', 'Location', 'best');
grid on;

% Plot 2: Parameter comparison matrix
subplot(2, 1, 2);

% Create a comparison of event counts
event_counts = zeros(size(param_sets, 1) + 1, 1);
param_labels = cell(size(param_sets, 1) + 1, 1);

% Current parameters
event_counts(1) = length(current_starts);
param_labels{1} = sprintf('Current (%.1fσ/%.1fσ)', ...
    config.event_detection.upper_threshold_sigma, ...
    config.event_detection.lower_threshold_sigma);

% Alternative parameters
for p = 1:size(param_sets, 1)
    event_counts(p+1) = length(alternative_results{p}.starts);
    param_labels{p+1} = sprintf('%.1fσ/%.1fσ', param_sets(p, 1), param_sets(p, 2));
end

% Create bar plot
bar(event_counts, 'FaceColor', [0.3 0.6 0.9], 'EdgeColor', 'k');
set(gca, 'XTickLabel', param_labels, 'XTickLabelRotation', 45);
ylabel('Number of Events Detected');
title(sprintf('Parameter Sensitivity Analysis - ROI %d', test_roi));
grid on;

% Add text labels on bars
for i = 1:length(event_counts)
    text(i, event_counts(i) + 0.1, sprintf('%d', event_counts(i)), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

%% === Summary Table ===
fprintf('\n=== PARAMETER SENSITIVITY SUMMARY ===\n');
fprintf('Parameters        Events    Avg Duration    Notes\n');
fprintf('----------------  --------  --------------  ---------------------\n');

% Current parameters
if ~isempty(current_starts)
    current_avg_dur = mean(current_ends - current_starts + 1);
    current_avg_dur_ms = current_avg_dur * 1000 / config.frame_rate;
else
    current_avg_dur = 0;
    current_avg_dur_ms = 0;
end

fprintf('Current %.1f/%.1f    %8d      %8.1f        Reference\n', ...
    config.event_detection.upper_threshold_sigma, ...
    config.event_detection.lower_threshold_sigma, ...
    length(current_starts), current_avg_dur_ms);

% Alternative parameters
for p = 1:size(param_sets, 1)
    if ~isempty(alternative_results{p}.starts)
        alt_avg_dur = mean(alternative_results{p}.ends - alternative_results{p}.starts + 1);
        alt_avg_dur_ms = alt_avg_dur * 1000 / config.frame_rate;
    else
        alt_avg_dur_ms = 0;
    end
    
    fprintf('Test    %.1f/%.1f    %8d      %8.1f        %s\n', ...
        param_sets(p, 1), param_sets(p, 2), ...
        length(alternative_results{p}.starts), alt_avg_dur_ms, ...
        get_sensitivity_note(param_sets(p, :), [config.event_detection.upper_threshold_sigma, config.event_detection.lower_threshold_sigma]));
end

%% === Recommendations ===
fprintf('\n=== PARAMETER RECOMMENDATIONS ===\n');

% Find the parameter set that gives reasonable event count
reasonable_range = [2, 15];  % Reasonable number of events for this ROI
good_params = [];

if length(current_starts) >= reasonable_range(1) && length(current_starts) <= reasonable_range(2)
    fprintf('✓ Current parameters (%.1fσ/%.1fσ) give reasonable event count (%d)\n', ...
        config.event_detection.upper_threshold_sigma, ...
        config.event_detection.lower_threshold_sigma, length(current_starts));
end

for p = 1:size(param_sets, 1)
    event_count = length(alternative_results{p}.starts);
    if event_count >= reasonable_range(1) && event_count <= reasonable_range(2)
        good_params(end+1) = p;
        fprintf('✓ Parameters %.1fσ/%.1fσ give reasonable event count (%d)\n', ...
            param_sets(p, 1), param_sets(p, 2), event_count);
    end
end

if isempty(good_params) && (length(current_starts) < reasonable_range(1) || length(current_starts) > reasonable_range(2))
    fprintf('⚠️  Consider adjusting parameters - current count (%d) may be too high/low\n', length(current_starts));
    fprintf('   Try more sensitive parameters if count is low, less sensitive if high\n');
end

fprintf('\n✅ Parameter testing complete for ROI %d!\n', test_roi);

%% === Helper Functions ===
function events = test_schmitt_trigger_single_roi(trace, upper_thresh, lower_thresh, debug_roi)
    % Test Schmitt trigger on single ROI with optional debug output
    
    events = false(size(trace));
    state = 'baseline';
    event_start = 0;
    validated_events = [];
    min_sustained_frames = 3;
    
    show_debug = (debug_roi > 0);
    if show_debug
        fprintf('\n--- Debug for ROI %d ---\n', debug_roi);
    end
    
    for frame = 1:length(trace)
        signal = trace(frame);
        if isnan(signal), continue; end
        
        if strcmp(state, 'baseline')
            if signal > upper_thresh
                state = 'validating';
                event_start = frame;
                if show_debug
                    fprintf('Potential START at frame %d: signal=%.6f > %.6f\n', ...
                        frame, signal, upper_thresh);
                end
            end
            
        elseif strcmp(state, 'validating')
            frames_since_start = frame - event_start + 1;
            
            if signal < lower_thresh
                if show_debug
                    fprintf('  VALIDATION FAILED at frame %d: signal=%.6f < %.6f after %d frames\n', ...
                        frame, signal, lower_thresh, frames_since_start);
                end
                state = 'baseline';
                event_start = 0;
                
            elseif frames_since_start >= min_sustained_frames
                state = 'validated_event';
                if show_debug
                    fprintf('  VALIDATION PASSED at frame %d: continuing event...\n', frame);
                end
            end
            
        elseif strcmp(state, 'validated_event')
            if signal < lower_thresh
                event_end = frame - 1;
                validated_events(end+1, :) = [event_start, event_end];
                duration = event_end - event_start + 1;
                if show_debug
                    fprintf('  EVENT END at frame %d: duration=%d frames\n', event_end, duration);
                end
                state = 'baseline';
                event_start = 0;
            end
        end
    end
    
    % Handle event extending to end
    if strcmp(state, 'validated_event') && event_start > 0
        validated_events(end+1, :) = [event_start, length(trace)];
    end
    
    % Create event mask
    for i = 1:size(validated_events, 1)
        events(validated_events(i, 1):validated_events(i, 2)) = true;
    end
    
    if show_debug
        fprintf('Total validated events: %d\n', size(validated_events, 1));
    end
end

function note = get_sensitivity_note(test_params, current_params)
    % Generate sensitivity note for parameter comparison
    
    upper_change = test_params(1) - current_params(1);
    lower_change = test_params(2) - current_params(2);
    
    if upper_change < 0 && lower_change <= 0
        note = 'More sensitive';
    elseif upper_change > 0 && lower_change >= 0
        note = 'Less sensitive';
    elseif upper_change < 0
        note = 'Easier start';
    elseif lower_change < 0
        note = 'Longer events';
    elseif upper_change > 0
        note = 'Harder start';
    elseif lower_change > 0
        note = 'Shorter events';
    else
        note = 'Mixed effect';
    end
end
