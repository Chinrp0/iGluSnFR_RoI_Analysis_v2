%% SIGMA THRESHOLD CALCULATION EXPLANATION
% Detailed walkthrough of how thresholds are calculated

clear; clc; close all;

fprintf('=== How Sigma Thresholds Are Calculated ===\n\n');

%% === Load Example Data ===
test_file = 'D:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglusnfr4f_NGR\Spont\GPU_SNR_Processed\5_raw_mean\CP_Snfr4-NGR_Doc2b-R213W_Cs1-c1_spont-01_mean.csv';

try
    loader = csv_loader_v2();
    [data, metadata] = loader.loadSingleFile(test_file);
    config = tracenorm_config();
    
    [baseline, outlier_mask, baseline_stats] = baseline_detector(data, config);
    [dfof_data, dfof_stats] = dfof_calculator(data, baseline, config);
    
    fprintf('✓ Data loaded for threshold calculation demonstration\n\n');
catch ME
    fprintf('✗ Error loading data: %s\n', ME.message);
    return;
end

%% === Example ROI for Detailed Walkthrough ===
example_roi = 627;
roi_trace = dfof_data(:, example_roi);
time_vector = (1:length(roi_trace)) / config.frame_rate;

fprintf('=== STEP-BY-STEP THRESHOLD CALCULATION ===\n');
fprintf('Using ROI %d as example\n\n', example_roi);

%% === STEP 1: Start with dF/F Data ===
fprintf('STEP 1: Start with dF/F normalized data\n');
fprintf('  • Raw fluorescence has been baseline-corrected\n');
fprintf('  • Converted to dF/F = (F - F0) / F0\n');
fprintf('  • Trace length: %d frames (%.1f seconds)\n', length(roi_trace), length(roi_trace)/config.frame_rate);
fprintf('  • dF/F range: [%.6f, %.6f]\n', min(roi_trace), max(roi_trace));
fprintf('  • Zero line represents baseline fluorescence\n\n');

%% === STEP 2: Remove Invalid Data ===
fprintf('STEP 2: Clean the data\n');
valid_data = roi_trace(~isnan(roi_trace));
fprintf('  • Original frames: %d\n', length(roi_trace));
fprintf('  • Valid frames: %d\n', length(valid_data));
fprintf('  • Removed %d NaN values\n\n', length(roi_trace) - length(valid_data));

%% === STEP 3: Baseline-Only Noise Estimation (Research Paper Method) ===
fprintf('STEP 3: Estimate noise from baseline periods only\n');
fprintf('  • KEY INSIGHT: Don''t include event periods in noise calculation\n');
fprintf('  • Method: Use 70th percentile approach\n\n');

% Calculate 70th percentile threshold
baseline_70th = prctile(valid_data, 70);
fprintf('  70th percentile of dF/F: %.6f\n', baseline_70th);
fprintf('  → This value separates likely baseline (bottom 70%%) from likely events (top 30%%)\n\n');

% Separate baseline periods from event periods
baseline_data = valid_data(valid_data <= baseline_70th);
event_like_data = valid_data(valid_data > baseline_70th);

fprintf('  Data separation:\n');
fprintf('    • Baseline periods: %d frames (%.1f%% of trace)\n', ...
    length(baseline_data), 100*length(baseline_data)/length(valid_data));
fprintf('    • Event-like periods: %d frames (%.1f%% of trace)\n', ...
    length(event_like_data), 100*length(event_like_data)/length(valid_data));
fprintf('    • Baseline range: [%.6f, %.6f]\n', min(baseline_data), max(baseline_data));
fprintf('    • Event-like range: [%.6f, %.6f]\n\n', min(event_like_data), max(event_like_data));

%% === STEP 4: Calculate Noise Using MAD (Median Absolute Deviation) ===
fprintf('STEP 4: Calculate noise standard deviation\n');
fprintf('  • Use MAD (Median Absolute Deviation) for robust noise estimation\n');
fprintf('  • MAD is less sensitive to outliers than standard deviation\n\n');

% Calculate MAD of baseline data
baseline_median = median(baseline_data);
absolute_deviations = abs(baseline_data - baseline_median);
mad_value = median(absolute_deviations);

fprintf('  Baseline statistics:\n');
fprintf('    • Baseline median: %.6f\n', baseline_median);
fprintf('    • MAD (median absolute deviation): %.6f\n', mad_value);

% Convert MAD to standard deviation equivalent
mad_to_std_factor = 1.4826;  % For normally distributed data
noise_std = mad_value * mad_to_std_factor;

fprintf('    • MAD → std conversion factor: %.4f\n', mad_to_std_factor);
fprintf('    • Estimated noise std (σ): %.6f\n\n', noise_std);

fprintf('  Why MAD instead of regular std?\n');
fprintf('    • Regular std includes ALL baseline variation\n');
fprintf('    • MAD focuses on median, ignoring extreme values\n');
fprintf('    • More robust against small events missed by 70th percentile\n\n');

%% === STEP 5: Calculate Threshold Values ===
fprintf('STEP 5: Calculate actual threshold values\n');

% Example with common sigma values
sigma_examples = [1.5, 2.0, 3.0, 3.5, 4.0, 5.0];

fprintf('  Formula: Threshold = σ_multiplier × noise_std\n');
fprintf('  Noise std = %.6f\n\n', noise_std);

fprintf('  Sigma    Threshold     Above baseline    Biological meaning\n');
fprintf('  -----    ---------     --------------    ------------------\n');

for i = 1:length(sigma_examples)
    sigma = sigma_examples(i);
    threshold = sigma * noise_std;
    frames_above = sum(valid_data > threshold);
    percent_above = 100 * frames_above / length(valid_data);
    
    % Biological interpretation
    if sigma <= 2.0
        meaning = 'Very sensitive (lots of noise)';
    elseif sigma <= 2.5
        meaning = 'Sensitive (some noise)';
    elseif sigma <= 3.5
        meaning = 'Moderate (good balance)';
    elseif sigma <= 4.5
        meaning = 'Conservative (miss small events)';
    else
        meaning = 'Very conservative (only large events)';
    end
    
    fprintf('  %4.1f    %9.6f     %6d (%4.1f%%)    %s\n', ...
        sigma, threshold, frames_above, percent_above, meaning);
end

fprintf('\n');

%% === STEP 6: Example with Your Parameters (3.5σ/1.5σ) ===
fprintf('STEP 6: Your specific parameters (3.5σ upper, 1.5σ lower)\n');

upper_sigma = 3.5;
lower_sigma = 1.5;

upper_threshold = upper_sigma * noise_std;
lower_threshold = lower_sigma * noise_std;

fprintf('  Upper threshold (%.1fσ): %.6f\n', upper_sigma, upper_threshold);
fprintf('  Lower threshold (%.1fσ): %.6f\n', lower_sigma, lower_threshold);
fprintf('  Hysteresis gap: %.6f (%.1fx noise)\n', ...
    upper_threshold - lower_threshold, (upper_threshold - lower_threshold)/noise_std);

% Calculate frame statistics
frames_above_upper = sum(valid_data > upper_threshold);
frames_above_lower = sum(valid_data > lower_threshold);

fprintf('\n  Frame statistics:\n');
fprintf('    • Frames above upper (%.1fσ): %d (%.1f%% of trace)\n', ...
    upper_sigma, frames_above_upper, 100*frames_above_upper/length(valid_data));
fprintf('    • Frames above lower (%.1fσ): %d (%.1f%% of trace)\n', ...
    lower_sigma, frames_above_lower, 100*frames_above_lower/length(valid_data));
fprintf('    • Potential event frames: %d\n', frames_above_lower - frames_above_upper);

%% === VISUALIZATION ===
fprintf('\n=== Creating Threshold Calculation Visualization ===\n');

fig = figure('Name', sprintf('Threshold Calculation - ROI %d', example_roi), ...
    'Position', [100, 100, 1600, 1000]);

% Plot 1: Full trace with thresholds
subplot(2, 2, 1);
plot(time_vector, roi_trace, 'k-', 'LineWidth', 1);
hold on;

% Add threshold lines
line([time_vector(1), time_vector(end)], [upper_threshold, upper_threshold], ...
    'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
line([time_vector(1), time_vector(end)], [lower_threshold, lower_threshold], ...
    'Color', 'b', 'LineStyle', '--', 'LineWidth', 2);
line([time_vector(1), time_vector(end)], [baseline_70th, baseline_70th], ...
    'Color', 'g', 'LineStyle', ':', 'LineWidth', 1);
line([time_vector(1), time_vector(end)], [0, 0], ...
    'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 1);

title(sprintf('ROI %d: Calculated Thresholds', example_roi));
xlabel('Time (s)'); ylabel('dF/F');
legend('dF/F trace', sprintf('Upper (%.1fσ)', upper_sigma), sprintf('Lower (%.1fσ)', lower_sigma), ...
    '70th percentile', 'Zero', 'Location', 'best');
grid on;

% Plot 2: Histogram of dF/F values
subplot(2, 2, 2);
histogram(baseline_data, 50, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.7, 'EdgeColor', 'none');
hold on;
if ~isempty(event_like_data)
    histogram(event_like_data, 50, 'FaceColor', [1 0.5 0.5], 'FaceAlpha', 0.7, 'EdgeColor', 'none');
end

% Add threshold lines
xline(upper_threshold, 'r--', sprintf('%.1fσ', upper_sigma), 'LineWidth', 2, 'LabelOrientation', 'horizontal');
xline(lower_threshold, 'b--', sprintf('%.1fσ', lower_sigma), 'LineWidth', 2, 'LabelOrientation', 'horizontal');
xline(baseline_70th, 'g:', '70th %ile', 'LineWidth', 1, 'LabelOrientation', 'horizontal');

title('Distribution of dF/F Values');
xlabel('dF/F'); ylabel('Count');
legend('Baseline periods', 'Event-like periods', 'Location', 'best');
grid on;

% Plot 3: Noise estimation process
subplot(2, 2, 3);
plot(1:length(baseline_data), sort(baseline_data), 'b-', 'LineWidth', 1);
hold on;
plot(1:length(baseline_data), baseline_median * ones(size(baseline_data)), 'r--', 'LineWidth', 2);

title('Baseline Data for Noise Estimation');
xlabel('Sorted Frame Index'); ylabel('dF/F');
legend('Sorted baseline values', sprintf('Median: %.6f', baseline_median), 'Location', 'best');
grid on;

% Plot 4: Threshold sensitivity
subplot(2, 2, 4);
test_sigmas = 1:0.1:5;
test_thresholds = test_sigmas * noise_std;
frames_above = arrayfun(@(t) sum(valid_data > t), test_thresholds);
percent_above = 100 * frames_above / length(valid_data);

plot(test_sigmas, percent_above, 'k-', 'LineWidth', 2);
hold on;
plot(upper_sigma, 100*frames_above_upper/length(valid_data), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot(lower_sigma, 100*frames_above_lower/length(valid_data), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');

title('Threshold Sensitivity');
xlabel('Sigma Level'); ylabel('% Frames Above Threshold');
legend('Sensitivity curve', sprintf('Your upper (%.1fσ)', upper_sigma), ...
    sprintf('Your lower (%.1fσ)', lower_sigma), 'Location', 'best');
grid on;

sgtitle(sprintf('Threshold Calculation Process - ROI %d (σ = %.6f)', example_roi, noise_std), ...
    'FontSize', 14, 'FontWeight', 'bold');

%% === KEY CONCEPTS SUMMARY ===
fprintf('\n=== KEY CONCEPTS SUMMARY ===\n');
fprintf('1. BASELINE-ONLY NOISE ESTIMATION:\n');
fprintf('   • Uses 70th percentile to exclude likely events\n');
fprintf('   • Follows research paper methodology\n');
fprintf('   • More accurate than using full trace\n\n');

fprintf('2. MAD vs STANDARD DEVIATION:\n');
fprintf('   • MAD = Median Absolute Deviation\n');
fprintf('   • More robust against outliers\n');
fprintf('   • Converted to std using factor 1.4826\n\n');

fprintf('3. SIGMA MULTIPLIERS:\n');
fprintf('   • Higher σ = higher threshold = fewer events\n');
fprintf('   • Lower σ = lower threshold = more events\n');
fprintf('   • Balance between sensitivity and noise rejection\n\n');

fprintf('4. SCHMITT TRIGGER LOGIC:\n');
fprintf('   • Event STARTS when crossing ABOVE upper threshold\n');
fprintf('   • Event CONTINUES until dropping BELOW lower threshold\n');
fprintf('   • Hysteresis prevents false events from noise\n\n');

fprintf('5. YOUR PARAMETERS (%.1fσ/%.1fσ):\n', upper_sigma, lower_sigma);
fprintf('   • Upper: %.6f → %.1f%% of frames exceed\n', upper_threshold, 100*frames_above_upper/length(valid_data));
fprintf('   • Lower: %.6f → %.1f%% of frames exceed\n', lower_threshold, 100*frames_above_lower/length(valid_data));
fprintf('   • Gap: %.6f → good hysteresis for noise rejection\n', upper_threshold - lower_threshold);

fprintf('\n✅ Threshold calculation explanation complete!\n');
