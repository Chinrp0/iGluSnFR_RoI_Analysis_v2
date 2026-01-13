%% QUICK SNR FIX TEST - Test corrected SNR calculation
% This fixes the SNR calculation bug and shows realistic values

clear; clc;

fprintf('=== Testing Fixed SNR Calculation ===\n');

%% Load your processed data
addpath(genpath(pwd));
test_file = 'D:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglusnfr4f_NGR\Spont\GPU_SNR_Processed\5_raw_mean\CP_Snfr4-NGR_Doc2b-R213W_Cs1-c1_spont-01_mean.csv';

loader = csv_loader_v2();
[data, metadata] = loader.loadSingleFile(test_file);
config = tracenorm_config();

[baseline, outlier_mask, baseline_stats] = baseline_detector(data, config);
[dfof_data, dfof_stats] = dfof_calculator(data, baseline, config);

%% Calculate CORRECTED SNR using proper biological metrics
fprintf('Calculating corrected SNR values...\n');

[numFrames, numROIs] = size(dfof_data);

% Method 1: Signal = Peak response amplitude, Noise = baseline dF/F noise
peak_response = max(dfof_data, [], 1);  % Peak dF/F per ROI
baseline_noise = std(dfof_data(1:100, :), 0, 1);  % Noise in first 100 frames (should be near baseline)
corrected_snr_v1 = peak_response ./ baseline_noise;

% Method 2: Signal = Event amplitude above baseline, Noise = inter-event noise  
median_dfof = median(dfof_data, 1);  % Baseline dF/F level
signal_amplitude = peak_response - median_dfof;  % Signal above baseline
noise_estimate = mad(dfof_data, 1, 1) * 1.4826;  % Robust noise estimate using MAD
corrected_snr_v2 = signal_amplitude ./ noise_estimate;

% Method 3: Classic event detection SNR
event_threshold = 2 * std(dfof_data, 0, 1);  % 2-sigma threshold
above_threshold = max(dfof_data - event_threshold, [], 1);  % Signal above threshold
corrected_snr_v3 = above_threshold ./ std(dfof_data, 0, 1);

%% Compare original vs corrected SNR
fprintf('\n=== SNR Comparison ===\n');
fprintf('Original SNR (buggy):\n');
fprintf('  Range: %.3f - %.3f\n', min(dfof_stats.snr), max(dfof_stats.snr));
fprintf('  Mean: %.3f\n', mean(dfof_stats.snr));

fprintf('\nCorrected SNR v1 (Peak/Noise):\n');
fprintf('  Range: %.1f - %.1f\n', min(corrected_snr_v1), max(corrected_snr_v1));
fprintf('  Mean: %.1f\n', mean(corrected_snr_v1));
fprintf('  ROIs with SNR > 2: %d (%.1f%%)\n', sum(corrected_snr_v1 > 2), 100*mean(corrected_snr_v1 > 2));
fprintf('  ROIs with SNR > 5: %d (%.1f%%)\n', sum(corrected_snr_v1 > 5), 100*mean(corrected_snr_v1 > 5));

fprintf('\nCorrected SNR v2 (Signal/MAD):\n');
fprintf('  Range: %.1f - %.1f\n', min(corrected_snr_v2), max(corrected_snr_v2));
fprintf('  Mean: %.1f\n', mean(corrected_snr_v2));
fprintf('  ROIs with SNR > 2: %d (%.1f%%)\n', sum(corrected_snr_v2 > 2), 100*mean(corrected_snr_v2 > 2));

fprintf('\nCorrected SNR v3 (Event detection):\n');
fprintf('  Range: %.1f - %.1f\n', min(corrected_snr_v3), max(corrected_snr_v3));
fprintf('  Mean: %.1f\n', mean(corrected_snr_v3));
fprintf('  ROIs with SNR > 2: %d (%.1f%%)\n', sum(corrected_snr_v3 > 2), 100*mean(corrected_snr_v3 > 2));

%% Test filtering with corrected SNR
fprintf('\n=== Testing Filtering with Corrected SNR ===\n');

% Use the most reasonable SNR (v2 - signal amplitude over robust noise)
corrected_snr = corrected_snr_v2;

% Your original criteria with corrected SNR
events_per_roi = dfof_stats.events_per_roi;
max_dfof = max(dfof_data, [], 1);
baseline_quality = baseline_stats.valid_fraction;

% Test different SNR thresholds
snr_thresholds = [1.0, 2.0, 3.0, 5.0];

for snr_thresh = snr_thresholds
    good_rois = find(events_per_roi >= 2 & ...
                     corrected_snr >= snr_thresh & ...
                     max_dfof >= 0.020 & ...
                     baseline_quality >= 0.8);
    
    fprintf('SNR >= %.1f: %d ROIs pass all criteria (%.1f%%)\n', ...
        snr_thresh, length(good_rois), 100*length(good_rois)/numROIs);
end

%% Show sample ROIs with corrected metrics
fprintf('\n=== Sample ROIs with Corrected Metrics ===\n');
sample_rois = [28, 79, 213, 246, 248];  % First 5 from your random list

fprintf('ROI\tEvents\tOld SNR\tNew SNR\tMax dF/F\tPasses Filter?\n');
fprintf('---\t------\t-------\t-------\t--------\t--------------\n');

for roi_idx = sample_rois
    old_snr = dfof_stats.snr(roi_idx);
    new_snr = corrected_snr(roi_idx);
    passes = events_per_roi(roi_idx) >= 2 && new_snr >= 2.0 && max_dfof(roi_idx) >= 0.020;
    
    if passes
        status = 'Yes';
    else
        status = 'No';
    end
    
    fprintf('%d\t%d\t%.3f\t\t%.1f\t%.3f\t\t%s\n', ...
        roi_idx, events_per_roi(roi_idx), old_snr, new_snr, max_dfof(roi_idx), status);
end

%% Recommended filtering criteria
fprintf('\n=== RECOMMENDED FILTERING CRITERIA ===\n');
fprintf('Based on corrected SNR, use these criteria:\n');
fprintf('  MIN_EVENTS = 5        (above median: %d events)\n', round(median(events_per_roi)));
fprintf('  MIN_SNR = 3.0         (corrected SNR - good signal quality)\n');
fprintf('  MIN_DFOF = 0.030      (3%% - clear biological signal)\n');
fprintf('  MIN_BASELINE_QUALITY = 0.9  (high quality baseline)\n');

recommended_rois = find(events_per_roi >= 5 & ...
                       corrected_snr >= 3.0 & ...
                       max_dfof >= 0.030 & ...
                       baseline_quality >= 0.9);

fprintf('\nRecommended criteria would select: %d ROIs (%.1f%%)\n', ...
    length(recommended_rois), 100*length(recommended_rois)/numROIs);

if length(recommended_rois) >= 10
    fprintf('First 10 recommended ROIs: [%s]\n', sprintf('%d ', recommended_rois(1:10)));
end

%% Save corrected SNR for use in other scripts
fprintf('\n=== Saving Corrected SNR ===\n');
fprintf('Saving corrected_snr variable for use in other scripts...\n');
save('corrected_snr_data.mat', 'corrected_snr', 'corrected_snr_v1', 'corrected_snr_v2', 'corrected_snr_v3');
fprintf('Saved to: corrected_snr_data.mat\n');

fprintf('\n✅ SNR calculation fixed! Your filtering should work now.\n');