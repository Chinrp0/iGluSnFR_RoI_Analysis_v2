%% NOISE ESTIMATION METHODS COMPARISON
% Compare different approaches to calculating sigma for Schmitt trigger
% Shows impact on event detection sensitivity

clear; clc; close all;

fprintf('=== Noise Estimation Methods Comparison ===\n');

%% === Load Data ===
addpath(genpath(pwd));
test_file = 'D:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglusnfr4f_NGR\Spont\GPU_SNR_Processed\5_raw_mean\CP_Snfr4-NGR_Doc2b-R213W_Cs1-c1_spont-01_mean.csv';

loader = csv_loader_v2();
[data, metadata] = loader.loadSingleFile(test_file);
config = tracenorm_config();

[baseline, outlier_mask, baseline_stats] = baseline_detector(data, config);
[dfof_data, dfof_stats] = dfof_calculator(data, baseline, config);

time_vector = (1:size(data, 1)) / config.frame_rate;

fprintf('Data loaded: %d frames × %d ROIs\n', size(data, 1), size(data, 2));

%% === Compare Different Noise Estimation Methods ===

% Select a representative ROI for detailed comparison
roi_idx = 700;  % Change this to test different ROIs
roi_trace = dfof_data(:, roi_idx);

fprintf('\nTesting noise estimation methods on ROI %d:\n', roi_idx);

%% Method 1: Current approach (entire trace)
noise_std_entire = mad(roi_trace, 1) * 1.4826;
fprintf('1. Entire trace MAD:     σ = %.4f\n', noise_std_entire);

%% Method 2: Bottom 70% percentile (exclude large events)
baseline_percentile = 70;
threshold_70 = prctile(roi_trace, baseline_percentile);
baseline_periods_70 = roi_trace <= threshold_70;
baseline_data_70 = roi_trace(baseline_periods_70);
noise_std_70 = mad(baseline_data_70, 1) * 1.4826;
fprintf('2. Bottom 70%% periods:   σ = %.4f (using %d/%d points)\n', ...
    noise_std_70, length(baseline_data_70), length(roi_trace));

%% Method 3: Bottom 50% percentile (more conservative)
threshold_50 = prctile(roi_trace, 50);
baseline_periods_50 = roi_trace <= threshold_50;
baseline_data_50 = roi_trace(baseline_periods_50);
noise_std_50 = mad(baseline_data_50, 1) * 1.4826;
fprintf('3. Bottom 50%% periods:   σ = %.4f (using %d/%d points)\n', ...
    noise_std_50, length(baseline_data_50), length(roi_trace));

%% Method 4: Initial quiet period (first 20%)
initial_frames = 1:round(0.2 * length(roi_trace));
initial_data = roi_trace(initial_frames);
noise_std_initial = std(initial_data, 'omitnan');
fprintf('4. Initial 20%% period:   σ = %.4f (using %d/%d points)\n', ...
    noise_std_initial, length(initial_data), length(roi_trace));

%% Method 5: Iterative exclusion
current_noise = mad(roi_trace, 1) * 1.4826;
for iter = 1:3
    threshold_iter = 3.0 * current_noise;  % 3-sigma threshold
    potential_events = roi_trace > threshold_iter;
    baseline_periods_iter = ~potential_events;
    
    if sum(baseline_periods_iter) > 20
        baseline_data_iter = roi_trace(baseline_periods_iter);
        current_noise = mad(baseline_data_iter, 1) * 1.4826;
    else
        break;
    end
end
noise_std_iterative = current_noise;
fprintf('5. Iterative exclusion:  σ = %.4f (using %d/%d points)\n', ...
    noise_std_iterative, sum(baseline_periods_iter), length(roi_trace));

%% === Calculate Thresholds for Each Method ===
upper_sigma = 3.0;
lower_sigma = 1.5;

thresholds = struct();
thresholds.entire = [upper_sigma * noise_std_entire, lower_sigma * noise_std_entire];
thresholds.percent70 = [upper_sigma * noise_std_70, lower_sigma * noise_std_70];
thresholds.percent50 = [upper_sigma * noise_std_50, lower_sigma * noise_std_50];
thresholds.initial = [upper_sigma * noise_std_initial, lower_sigma * noise_std_initial];
thresholds.iterative = [upper_sigma * noise_std_iterative, lower_sigma * noise_std_iterative];

fprintf('\nResulting Schmitt trigger thresholds (3.0σ/1.5σ):\n');
fprintf('Method                Upper     Lower\n');
fprintf('-------------------- --------- ---------\n');
fprintf('1. Entire trace      %.4f    %.4f\n', thresholds.entire);
fprintf('2. Bottom 70%%        %.4f    %.4f\n', thresholds.percent70);
fprintf('3. Bottom 50%%        %.4f    %.4f\n', thresholds.percent50);
fprintf('4. Initial period    %.4f    %.4f\n', thresholds.initial);
fprintf('5. Iterative         %.4f    %.4f\n', thresholds.iterative);

%% === Visual Comparison ===
fprintf('\nCreating visual comparison...\n');

fig1 = figure('Name', sprintf('Noise Estimation Methods - ROI %d', roi_idx), ...
    'Position', [100, 100, 1400, 800]);

methods = {'entire', 'percent70', 'percent50', 'initial', 'iterative'};
method_names = {'Entire Trace', 'Bottom 70%', 'Bottom 50%', 'Initial 20%', 'Iterative'};

for i = 1:length(methods)
    subplot(2, 3, i);
    
    % Plot the trace
    plot(time_vector, roi_trace, 'k-', 'LineWidth', 0.8);
    hold on;
    
    % Add threshold lines
    method = methods{i};
    upper_thresh = thresholds.(method)(1);
    lower_thresh = thresholds.(method)(2);
    
    yline(upper_thresh, 'g--', sprintf('Upper: %.4f', upper_thresh), 'LineWidth', 1.5);
    yline(lower_thresh, 'c--', sprintf('Lower: %.4f', lower_thresh), 'LineWidth', 1.5);
    yline(0, 'k:', 'Alpha', 0.5);
    
    % Highlight baseline periods used for this method
    switch method
        case 'percent70'
            baseline_mask = baseline_periods_70;
        case 'percent50'
            baseline_mask = baseline_periods_50;
        case 'initial'
            baseline_mask = false(size(roi_trace));
            baseline_mask(initial_frames) = true;
        case 'iterative'
            baseline_mask = baseline_periods_iter;
        otherwise
            baseline_mask = true(size(roi_trace));  % Entire trace
    end
    
    if ~strcmp(method, 'entire')
        % Shade baseline periods
        baseline_trace = roi_trace;
        baseline_trace(~baseline_mask) = NaN;
        plot(time_vector, baseline_trace, 'r.', 'MarkerSize', 4);
    end
    
    xlabel('Time (s)');
    ylabel('dF/F');
    title(sprintf('%s\nσ = %.4f', method_names{i}, ...
        thresholds.(method)(1) / upper_sigma));
    grid on;
    
    % Set consistent y-limits
    ylim([min(roi_trace)*1.2, max(roi_trace)*1.2]);
end

% Add comparison subplot
subplot(2, 3, 6);
noise_values = [noise_std_entire, noise_std_70, noise_std_50, noise_std_initial, noise_std_iterative];
bar(noise_values);
set(gca, 'XTickLabel', {'Entire', '70%', '50%', 'Initial', 'Iterative'});
ylabel('Noise σ (dF/F)');
title('Noise Estimates Comparison');
grid on;

% Add text annotation
for i = 1:length(noise_values)
    text(i, noise_values(i) + 0.001, sprintf('%.4f', noise_values(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
end

sgtitle(sprintf('Noise Estimation Methods Comparison - ROI %d', roi_idx), 'FontSize', 14);

%% === Event Detection Comparison ===
fprintf('\nComparing event detection with different noise estimates...\n');

% Apply simple threshold detection with each method
fig2 = figure('Name', 'Event Detection Sensitivity Comparison', ...
    'Position', [200, 100, 1200, 800]);

for i = 1:length(methods)
    subplot(2, 3, i);
    
    % Plot trace
    plot(time_vector, roi_trace, 'k-', 'LineWidth', 0.8);
    hold on;
    
    % Apply threshold
    method = methods{i};
    upper_thresh = thresholds.(method)(1);
    
    % Detect events above upper threshold
    events = roi_trace > upper_thresh;
    
    if any(events)
        event_trace = roi_trace;
        event_trace(~events) = NaN;
        plot(time_vector, event_trace, 'r-', 'LineWidth', 2);
    end
    
    yline(upper_thresh, 'g--', 'LineWidth', 1.5);
    yline(0, 'k:', 'Alpha', 0.5);
    
    xlabel('Time (s)');
    ylabel('dF/F');
    title(sprintf('%s: %d events', method_names{i}, sum(events)));
    grid on;
    
    ylim([min(roi_trace)*1.2, max(roi_trace)*1.2]);
end

% Summary comparison
subplot(2, 3, 6);
event_counts = zeros(size(methods));
for i = 1:length(methods)
    method = methods{i};
    upper_thresh = thresholds.(method)(1);
    events = roi_trace > upper_thresh;
    event_counts(i) = sum(events);
end

bar(event_counts);
set(gca, 'XTickLabel', {'Entire', '70%', '50%', 'Initial', 'Iterative'});
ylabel('Event Frames Detected');
title('Detection Sensitivity');
grid on;

for i = 1:length(event_counts)
    text(i, event_counts(i) + 5, sprintf('%d', event_counts(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
end

sgtitle('Event Detection Sensitivity by Noise Estimation Method', 'FontSize', 14);

%% === Recommendations ===
fprintf('\n=== RECOMMENDATIONS ===\n');
fprintf('1. CURRENT METHOD (Entire trace): σ = %.4f\n', noise_std_entire);
fprintf('   - Pros: Simple, robust to outliers\n');
fprintf('   - Cons: Includes events in noise calculation, reduces sensitivity\n\n');

fprintf('2. RECOMMENDED (Bottom 70%%): σ = %.4f\n', noise_std_70);
fprintf('   - Pros: Excludes large events, more accurate baseline noise\n');
fprintf('   - Cons: Slightly more complex\n');
fprintf('   - Sensitivity increase: %.1fx\n', noise_std_entire/noise_std_70);

fprintf('\n3. CONSERVATIVE (Bottom 50%%): σ = %.4f\n', noise_std_50);
fprintf('   - Pros: Very conservative baseline estimate\n');
fprintf('   - Cons: May be too sensitive, include noise as events\n');
fprintf('   - Sensitivity increase: %.1fx\n', noise_std_entire/noise_std_50);

fprintf('\n4. ITERATIVE EXCLUSION: σ = %.4f\n', noise_std_iterative);
fprintf('   - Pros: Adapts to data, excludes detected events from noise\n');
fprintf('   - Cons: More computationally complex\n');
fprintf('   - Sensitivity increase: %.1fx\n', noise_std_entire/noise_std_iterative);

fprintf('\n=== SUGGESTED IMPLEMENTATION ===\n');
fprintf('For your spontaneous release data, I recommend:\n');
fprintf('1. Use BOTTOM 70%% method as default\n');
fprintf('2. Add iterative refinement as an option\n');
fprintf('3. Keep current method as fallback for problematic ROIs\n');

fprintf('\nThis will make your Schmitt trigger more sensitive to real events\n');
fprintf('while avoiding false positives from baseline noise.\n');