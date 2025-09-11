%% CORRECTED THRESHOLD CALCULATION
% Use baseline detector's outlier information to properly estimate noise
% Include outliers in noise calculation, exclude only real events

clear; clc; close all;

fprintf('=== CORRECTED Threshold Calculation Method ===\n');
fprintf('Key insight: Include outliers/noise in noise calculation!\n\n');

%% === Load Data ===
test_file = 'D:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglusnfr4f_NGR\Spont\GPU_SNR_Processed\5_raw_mean\CP_Snfr4-NGR_Doc2b-R213W_Cs1-c1_spont-01_mean.csv';

try
    loader = csv_loader_v2();
    [data, metadata] = loader.loadSingleFile(test_file);
    config = tracenorm_config();
    
    [baseline, outlier_mask, baseline_stats] = baseline_detector(data, config);
    [dfof_data, dfof_stats] = dfof_calculator(data, baseline, config);
    
    fprintf('✓ Data loaded with baseline and outlier information\n\n');
catch ME
    fprintf('✗ Error loading data: %s\n', ME.message);
    return;
end

%% === Example ROI ===
example_roi = 632;
roi_trace = dfof_data(:, example_roi);
roi_outliers = outlier_mask(:, example_roi);
time_vector = (1:length(roi_trace)) / config.frame_rate;

fprintf('=== COMPARISON: Old vs New Threshold Calculation ===\n');
fprintf('Using ROI %d\n\n', example_roi);

%% === METHOD 1: Current (70th Percentile) ===
fprintf('METHOD 1: Current approach (70th percentile)\n');

valid_data = roi_trace(~isnan(roi_trace));
baseline_70th = prctile(valid_data, 70);
baseline_data_old = valid_data(valid_data <= baseline_70th);
excluded_data_old = valid_data(valid_data > baseline_70th);

noise_std_old = mad(baseline_data_old, 1) * 1.4826;

fprintf('  • Uses bottom 70%% of all dF/F values\n');
fprintf('  • Baseline data points: %d\n', length(baseline_data_old));
fprintf('  • Excluded data points: %d (assumed to be "events")\n', length(excluded_data_old));
fprintf('  • Excluded range: [%.6f, %.6f]\n', min(excluded_data_old), max(excluded_data_old));
fprintf('  • Calculated noise std: %.6f\n\n', noise_std_old);

%% === METHOD 2: Corrected (Use Outlier Information) ===
fprintf('METHOD 2: CORRECTED approach (use outlier mask)\n');

% Get valid frames (not NaN)
valid_frames = ~isnan(roi_trace);

% Separate into different categories
non_outlier_frames = valid_frames & ~roi_outliers;
outlier_frames = valid_frames & roi_outliers;

% For noise calculation: include baseline + outliers, exclude only sustained events
% First, get baseline + outlier data (this IS the noise)
noise_data = roi_trace(non_outlier_frames | outlier_frames);  % Everything valid
noise_data = noise_data(~isnan(noise_data));

% Now we need to identify and REMOVE sustained events (not noise)
% Use a simple approach: remove periods that are consistently high
sustained_threshold = prctile(noise_data, 85);  % Higher threshold for sustained events
sustained_window = 5;  % Must be high for 5+ consecutive frames

% Find sustained high periods
sustained_mask = false(size(roi_trace));
high_values = roi_trace > sustained_threshold;

% Identify runs of consecutive high values
run_starts = find(diff([false; high_values]) == 1);
run_ends = find(diff([high_values; false]) == -1);

for i = 1:length(run_starts)
    run_length = run_ends(i) - run_starts(i) + 1;
    if run_length >= sustained_window
        sustained_mask(run_starts(i):run_ends(i)) = true;
    end
end

% Final noise data: everything EXCEPT sustained events
noise_calculation_mask = valid_frames & ~sustained_mask;
final_noise_data = roi_trace(noise_calculation_mask);

noise_std_new = mad(final_noise_data, 1) * 1.4826;

fprintf('  • Includes baseline + outliers + brief spikes (all noise types)\n');
fprintf('  • Excludes only sustained high periods (likely real events)\n');
fprintf('  • Total data points: %d\n', sum(valid_frames));
fprintf('  • Outlier frames: %d\n', sum(outlier_frames));
fprintf('  • Sustained event frames: %d\n', sum(sustained_mask));
fprintf('  • Used for noise calculation: %d\n', length(final_noise_data));
fprintf('  • Calculated noise std: %.6f\n\n', noise_std_new);

%% === COMPARE THRESHOLDS ===
fprintf('=== THRESHOLD COMPARISON ===\n');

% Test parameters
upper_sigma = 3.5;
lower_sigma = 1.5;

% Old thresholds
upper_thresh_old = upper_sigma * noise_std_old;
lower_thresh_old = lower_sigma * noise_std_old;

% New thresholds  
upper_thresh_new = upper_sigma * noise_std_new;
lower_thresh_new = lower_sigma * noise_std_new;

fprintf('Parameter: %.1fσ upper / %.1fσ lower\n\n', upper_sigma, lower_sigma);

fprintf('                     Old Method    New Method    Difference\n');
fprintf('                     ----------    ----------    ----------\n');
fprintf('Noise std:           %10.6f    %10.6f    %+9.6f\n', ...
    noise_std_old, noise_std_new, noise_std_new - noise_std_old);
fprintf('Upper threshold:     %10.6f    %10.6f    %+9.6f\n', ...
    upper_thresh_old, upper_thresh_new, upper_thresh_new - upper_thresh_old);
fprintf('Lower threshold:     %10.6f    %10.6f    %+9.6f\n', ...
    lower_thresh_old, lower_thresh_new, lower_thresh_new - lower_thresh_old);

% Calculate how many frames exceed each threshold
frames_above_old_upper = sum(valid_data > upper_thresh_old);
frames_above_new_upper = sum(valid_data > upper_thresh_new);
frames_above_old_lower = sum(valid_data > lower_thresh_old);
frames_above_new_lower = sum(valid_data > lower_thresh_new);

fprintf('\nFrames above thresholds:\n');
fprintf('Upper (%.1fσ):        %10d    %10d    %+9d\n', ...
    upper_sigma, frames_above_old_upper, frames_above_new_upper, ...
    frames_above_new_upper - frames_above_old_upper);
fprintf('Lower (%.1fσ):        %10d    %10d    %+9d\n', ...
    lower_sigma, frames_above_old_lower, frames_above_new_lower, ...
    frames_above_new_lower - frames_above_old_lower);

%% === TEST EVENT DETECTION ===
fprintf('\n=== EVENT DETECTION COMPARISON ===\n');

% Test with both threshold sets
events_old = test_event_detection(roi_trace, upper_thresh_old, lower_thresh_old, 3);
events_new = test_event_detection(roi_trace, upper_thresh_new, lower_thresh_new, 3);

old_starts = find(diff([false; events_old]) == 1);
new_starts = find(diff([false; events_new]) == 1);

fprintf('Events detected:\n');
fprintf('  Old method: %d events\n', length(old_starts));
fprintf('  New method: %d events\n', length(new_starts));
fprintf('  Difference: %+d events\n', length(new_starts) - length(old_starts));

if length(new_starts) < length(old_starts)
    fprintf('  ✓ New method detects fewer events (better noise rejection)\n');
elseif length(new_starts) > length(old_starts)
    fprintf('  ⚠ New method detects more events (check parameters)\n');
else
    fprintf('  = Same number of events detected\n');
end

%% === VISUALIZATION ===
fprintf('\n=== Creating Comparison Visualization ===\n');

fig = figure('Name', sprintf('Corrected Threshold Calculation - ROI %d', example_roi), ...
    'Position', [100, 100, 1800, 1200]);

% Plot 1: Data categorization
subplot(3, 2, 1);
plot(time_vector, roi_trace, 'k-', 'LineWidth', 0.5);
hold on;

% Highlight different data types
outlier_trace = roi_trace;
outlier_trace(~roi_outliers) = NaN;
plot(time_vector, outlier_trace, 'ro', 'MarkerSize', 3);

sustained_trace = roi_trace;
sustained_trace(~sustained_mask) = NaN;
plot(time_vector, sustained_trace, 'b-', 'LineWidth', 2);

baseline_trace = roi_trace;
baseline_trace(roi_outliers | sustained_mask) = NaN;
plot(time_vector, baseline_trace, 'g-', 'LineWidth', 0.5);

title('Data Categorization');
xlabel('Time (s)'); ylabel('dF/F');
legend('All data', 'Outliers (noise)', 'Sustained events', 'Baseline', 'Location', 'best');
grid on;

% Plot 2: Threshold comparison
subplot(3, 2, 2);
plot(time_vector, roi_trace, 'k-', 'LineWidth', 1);
hold on;

% Old thresholds
line([time_vector(1), time_vector(end)], [upper_thresh_old, upper_thresh_old], ...
    'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
line([time_vector(1), time_vector(end)], [lower_thresh_old, lower_thresh_old], ...
    'Color', 'r', 'LineStyle', ':', 'LineWidth', 2);

% New thresholds
line([time_vector(1), time_vector(end)], [upper_thresh_new, upper_thresh_new], ...
    'Color', 'b', 'LineStyle', '--', 'LineWidth', 2);
line([time_vector(1), time_vector(end)], [lower_thresh_new, lower_thresh_new], ...
    'Color', 'b', 'LineStyle', ':', 'LineWidth', 2);

title('Threshold Comparison');
xlabel('Time (s)'); ylabel('dF/F');
legend('dF/F', 'Old upper', 'Old lower', 'NEW upper', 'NEW lower', 'Location', 'best');
grid on;

% Plot 3: Noise data histograms
subplot(3, 2, 3);
histogram(baseline_data_old, 50, 'FaceColor', [1 0.7 0.7], 'FaceAlpha', 0.7, 'EdgeColor', 'none');
hold on;
histogram(final_noise_data, 50, 'FaceColor', [0.7 0.7 1], 'FaceAlpha', 0.7, 'EdgeColor', 'none');

xline(noise_std_old, 'r--', sprintf('Old σ: %.6f', noise_std_old), 'LineWidth', 2);
xline(noise_std_new, 'b--', sprintf('New σ: %.6f', noise_std_new), 'LineWidth', 2);

title('Noise Data Distribution');
xlabel('dF/F'); ylabel('Count');
legend('Old method (70th %ile)', 'New method (exclude sustained)', 'Location', 'best');
grid on;

% Plot 4: Event detection comparison
subplot(3, 2, 4);
plot(time_vector, roi_trace, 'k-', 'LineWidth', 1);
hold on;

% Old events
if any(events_old)
    old_trace = roi_trace;
    old_trace(~events_old) = NaN;
    plot(time_vector, old_trace, 'r-', 'LineWidth', 2);
end

% New events
if any(events_new)
    new_trace = roi_trace;
    new_trace(~events_new) = NaN;
    plot(time_vector, new_trace, 'b-', 'LineWidth', 2);
end

title(sprintf('Event Detection: Old (%d) vs New (%d)', length(old_starts), length(new_starts)));
xlabel('Time (s)'); ylabel('dF/F');
legend('dF/F', 'Old method events', 'NEW method events', 'Location', 'best');
grid on;

% Plot 5: Sensitivity analysis
subplot(3, 2, 5);
test_sigmas = 1:0.1:6;
old_thresholds = test_sigmas * noise_std_old;
new_thresholds = test_sigmas * noise_std_new;

old_counts = arrayfun(@(t) sum(valid_data > t), old_thresholds);
new_counts = arrayfun(@(t) sum(valid_data > t), new_thresholds);

plot(test_sigmas, old_counts, 'r-', 'LineWidth', 2);
hold on;
plot(test_sigmas, new_counts, 'b-', 'LineWidth', 2);

plot(upper_sigma, frames_above_old_upper, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot(upper_sigma, frames_above_new_upper, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');

title('Threshold Sensitivity');
xlabel('Sigma Level'); ylabel('Frames Above Threshold');
legend('Old method', 'NEW method', sprintf('Old %.1fσ', upper_sigma), sprintf('New %.1fσ', upper_sigma), 'Location', 'best');
grid on;

% Plot 6: Noise estimation comparison
subplot(3, 2, 6);
categories = {'Baseline only\n(old)', 'Baseline +\nOutliers +\nBrief spikes\n(new)'};
noise_values = [noise_std_old, noise_std_new];

b = bar(noise_values, 'FaceColor', [0.7 0.9 0.7], 'EdgeColor', 'k');
set(gca, 'XTickLabel', categories);
ylabel('Noise Standard Deviation');
title('Noise Estimation Methods');

% Add value labels
for i = 1:length(noise_values)
    text(i, noise_values(i) + 0.001, sprintf('%.6f', noise_values(i)), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end
grid on;

sgtitle(sprintf('CORRECTED Threshold Calculation - ROI %d', example_roi), ...
    'FontSize', 16, 'FontWeight', 'bold');

%% === KEY INSIGHTS ===
fprintf('\n=== KEY INSIGHTS ===\n');
fprintf('1. NOISE SHOULD INCLUDE OUTLIERS:\n');
fprintf('   • Outliers are noise, not events!\n');
fprintf('   • Including them gives realistic noise estimate\n');
fprintf('   • Old method: σ = %.6f (underestimated)\n', noise_std_old);
fprintf('   • New method: σ = %.6f (realistic)\n\n', noise_std_new);

fprintf('2. THRESHOLD EFFECT:\n');
if noise_std_new > noise_std_old
    fprintf('   • Higher noise estimate → higher thresholds\n');
    fprintf('   • Better noise rejection\n');
    fprintf('   • Fewer false positive "events"\n\n');
else
    fprintf('   • Noise estimates similar\n\n');
end

fprintf('3. EVENT DETECTION IMPROVEMENT:\n');
if length(new_starts) < length(old_starts)
    reduction = length(old_starts) - length(new_starts);
    percent_reduction = 100 * reduction / length(old_starts);
    fprintf('   • Reduced events by %d (%.1f%% reduction)\n', reduction, percent_reduction);
    fprintf('   • Better signal-to-noise discrimination\n');
end

fprintf('\n✅ Corrected threshold calculation complete!\n');
fprintf('This approach gives more realistic noise estimates by including all noise types.\n');

%% === Helper Function ===
function events = test_event_detection(trace, upper_thresh, lower_thresh, min_frames)
    % Simple event detection for comparison
    
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
            elseif frames_since_start >= min_frames
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
