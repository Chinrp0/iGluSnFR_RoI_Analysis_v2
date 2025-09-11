%% RANDOM ROI VIEWER - Show random ROIs to guide filtering decisions
% No filtering - just show what the data actually looks like

clear; clc; close all;

%% === Configuration ===
NUM_TRACES = 20;  % How many random ROIs to show
SHOW_DFOF_PLOTS = true;
SHOW_STATISTICS = true;

fprintf('=== Random ROI Viewer ===\n');
fprintf('Showing %d random ROIs to guide filtering decisions\n\n', NUM_TRACES);

%% === Load Data ===
addpath(genpath(pwd));
test_file = 'D:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglusnfr4f_NGR\Spont\GPU_SNR_Processed\5_raw_mean\CP_Snfr4-NGR_Doc2b-R213W_Cs1-c1_spont-01_mean.csv';

loader = csv_loader_v2();
[data, metadata] = loader.loadSingleFile(test_file);
config = tracenorm_config();

[baseline, outlier_mask, baseline_stats] = baseline_detector(data, config);
[dfof_data, dfof_stats] = dfof_calculator(data, baseline, config);

time_vector = (1:size(data, 1)) / config.frame_rate;
[numFrames, numROIs] = size(data);

fprintf('Data loaded: %d frames × %d ROIs\n', numFrames, numROIs);

%% === Select Random ROIs ===
rng(42);  % Set seed for reproducible random selection
random_rois = randperm(numROIs, NUM_TRACES);
random_rois = sort(random_rois);  % Sort for easier reference

fprintf('Randomly selected ROIs: [%s]\n', sprintf('%d ', random_rois));

%% === Extract Metrics for All ROIs ===
events_per_roi = dfof_stats.events_per_roi;
snr_values = dfof_stats.snr;
max_dfof = max(dfof_data, [], 1);
min_dfof = min(dfof_data, [], 1);
baseline_quality = baseline_stats.valid_fraction;
amplitude_range = max_dfof - min_dfof;

%% === Show Data Distribution First ===
if SHOW_STATISTICS
    fprintf('\n=== Overall Data Statistics (All %d ROIs) ===\n', numROIs);
    
    % Events per ROI
    fprintf('Events per ROI:\n');
    fprintf('  Range: %d - %d events\n', min(events_per_roi), max(events_per_roi));
    fprintf('  Mean: %.1f events\n', mean(events_per_roi));
    fprintf('  Median: %.1f events\n', median(events_per_roi));
    fprintf('  ROIs with >0 events: %d (%.1f%%)\n', sum(events_per_roi > 0), 100*mean(events_per_roi > 0));
    fprintf('  ROIs with >5 events: %d (%.1f%%)\n', sum(events_per_roi > 5), 100*mean(events_per_roi > 5));
    fprintf('  ROIs with >10 events: %d (%.1f%%)\n', sum(events_per_roi > 10), 100*mean(events_per_roi > 10));
    
    % SNR values
    fprintf('\nSNR values:\n');
    fprintf('  Range: %.2f - %.2f\n', min(snr_values), max(snr_values));
    fprintf('  Mean: %.2f\n', mean(snr_values));
    fprintf('  Median: %.2f\n', median(snr_values));
    fprintf('  ROIs with SNR > 1.0: %d (%.1f%%)\n', sum(snr_values > 1.0), 100*mean(snr_values > 1.0));
    fprintf('  ROIs with SNR > 2.0: %d (%.1f%%)\n', sum(snr_values > 2.0), 100*mean(snr_values > 2.0));
    
    % dF/F values
    fprintf('\nMax dF/F values:\n');
    fprintf('  Range: %.4f - %.4f\n', min(max_dfof), max(max_dfof));
    fprintf('  Mean: %.4f\n', mean(max_dfof));
    fprintf('  Median: %.4f\n', median(max_dfof));
    fprintf('  ROIs with max dF/F > 0.01: %d (%.1f%%)\n', sum(max_dfof > 0.01), 100*mean(max_dfof > 0.01));
    fprintf('  ROIs with max dF/F > 0.02: %d (%.1f%%)\n', sum(max_dfof > 0.02), 100*mean(max_dfof > 0.02));
    fprintf('  ROIs with max dF/F > 0.05: %d (%.1f%%)\n', sum(max_dfof > 0.05), 100*mean(max_dfof > 0.05));
    
    % Baseline quality
    fprintf('\nBaseline quality:\n');
    fprintf('  Range: %.2f - %.2f\n', min(baseline_quality), max(baseline_quality));
    fprintf('  Mean: %.3f\n', mean(baseline_quality));
    fprintf('  ROIs with quality > 0.8: %d (%.1f%%)\n', sum(baseline_quality > 0.8), 100*mean(baseline_quality > 0.8));
    fprintf('  ROIs with quality > 0.9: %d (%.1f%%)\n', sum(baseline_quality > 0.9), 100*mean(baseline_quality > 0.9));
end

%% === Create Raw + Baseline Plots ===
fprintf('\n=== Creating plots for random ROIs ===\n');

rows = ceil(sqrt(NUM_TRACES));
cols = ceil(NUM_TRACES / rows);

fig1 = figure('Name', 'Random ROIs - Raw + Baseline', 'Position', [50, 50, 1600, 1200]);

for i = 1:NUM_TRACES
    roi_idx = random_rois(i);
    
    subplot(rows, cols, i);
    
    % Plot raw data
    plot(time_vector, data(:, roi_idx), 'k-', 'LineWidth', 0.6, 'DisplayName', 'Raw');
    hold on;
    
    % Plot baseline
    plot(time_vector, baseline(:, roi_idx), 'r-', 'LineWidth', 1.2, 'DisplayName', 'Baseline');
    
    % Mark detected events (outliers)
    events = outlier_mask(:, roi_idx);
    if any(events)
        scatter(time_vector(events), data(events, roi_idx), 15, 'o', ...
            'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor', [1 0.5 0], 'MarkerSize', 4);
    end
    
    xlabel('Time (s)', 'FontSize', 8);
    ylabel('Fluorescence', 'FontSize', 8);
    title(sprintf('ROI %d: %d events, SNR=%.2f', roi_idx, events_per_roi(roi_idx), snr_values(roi_idx)), 'FontSize', 9);
    grid on;
    
    % Add text with key metrics
    text(0.02, 0.98, sprintf('dF/F: %.3f', max_dfof(roi_idx)), 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'FontSize', 7, 'BackgroundColor', 'white');
    
    if i == 1
        legend('Raw', 'Baseline', 'Events', 'Location', 'best', 'FontSize', 7);
    end
end

sgtitle(sprintf('Random ROI Sample - Raw Traces (%s)', metadata.filename), 'FontSize', 14);

%% === Create dF/F Plots ===
if SHOW_DFOF_PLOTS
    fig2 = figure('Name', 'Random ROIs - dF/F Traces', 'Position', [100, 100, 1600, 1200]);
    
    for i = 1:NUM_TRACES
        roi_idx = random_rois(i);
        
        subplot(rows, cols, i);
        
        % Plot dF/F
        plot(time_vector, dfof_data(:, roi_idx), 'b-', 'LineWidth', 0.8);
        hold on;
        
        % Add reference lines
        yline(0, 'k--', 'Alpha', 0.3);
        yline(0.01, 'g:', 'Alpha', 0.5);  % 1%
        yline(0.05, 'r:', 'Alpha', 0.5);  % 5%
        
        % Mark significant events
        big_events = dfof_data(:, roi_idx) > 0.01;  % 1% threshold
        if any(big_events)
            scatter(time_vector(big_events), dfof_data(big_events, roi_idx), 10, ...
                'filled', 'MarkerFaceColor', [0 0.8 0]);
        end
        
        xlabel('Time (s)', 'FontSize', 8);
        ylabel('dF/F', 'FontSize', 8);
        title(sprintf('ROI %d dF/F (range: %.3f to %.3f)', roi_idx, min_dfof(roi_idx), max_dfof(roi_idx)), 'FontSize', 9);
        grid on;
        
        % Set reasonable y-axis limits
        y_range = max_dfof(roi_idx) - min_dfof(roi_idx);
        y_center = (max_dfof(roi_idx) + min_dfof(roi_idx)) / 2;
        y_margin = max(y_range * 0.1, 0.01);  % At least 1% margin
        ylim([y_center - y_range/2 - y_margin, y_center + y_range/2 + y_margin]);
    end
    
    sgtitle('Random ROI Sample - dF/F Traces', 'FontSize', 14);
end

%% === Show Statistics for Selected Random ROIs ===
fprintf('\n=== Statistics for Displayed Random ROIs ===\n');
fprintf('ROI\tEvents\tSNR\tMax dF/F\tMin dF/F\tRange\tBaseline Quality\n');
fprintf('---\t------\t----\t--------\t--------\t-----\t----------------\n');

for i = 1:NUM_TRACES
    roi_idx = random_rois(i);
    fprintf('%d\t%d\t%.2f\t%.4f\t\t%.4f\t\t%.4f\t%.1f%%\n', ...
        roi_idx, events_per_roi(roi_idx), snr_values(roi_idx), ...
        max_dfof(roi_idx), min_dfof(roi_idx), amplitude_range(roi_idx), ...
        100*baseline_quality(roi_idx));
end

%% === Suggested Filtering Criteria Based on Data ===
fprintf('\n=== SUGGESTED FILTERING CRITERIA ===\n');
fprintf('Based on the data distribution, consider these thresholds:\n\n');

% Conservative criteria (top 20%)
events_80th = prctile(events_per_roi, 80);
snr_80th = prctile(snr_values, 80);
dfof_80th = prctile(max_dfof, 80);

% Moderate criteria (top 50%)
events_50th = prctile(events_per_roi, 50);
snr_50th = prctile(snr_values, 50);
dfof_50th = prctile(max_dfof, 50);

% Liberal criteria (top 70%)
events_30th = prctile(events_per_roi, 30);
snr_30th = prctile(snr_values, 30);
dfof_30th = prctile(max_dfof, 30);

fprintf('CONSERVATIVE (top 20%% of ROIs):\n');
fprintf('  MIN_EVENTS = %d\n', round(events_80th));
fprintf('  MIN_SNR = %.2f\n', snr_80th);
fprintf('  MIN_DFOF = %.4f\n', dfof_80th);
conservative_count = sum(events_per_roi >= events_80th & snr_values >= snr_80th & max_dfof >= dfof_80th);
fprintf('  → Would select ~%d ROIs\n', conservative_count);

fprintf('\nMODERATE (median values):\n');
fprintf('  MIN_EVENTS = %d\n', round(events_50th));
fprintf('  MIN_SNR = %.2f\n', snr_50th);
fprintf('  MIN_DFOF = %.4f\n', dfof_50th);
moderate_count = sum(events_per_roi >= events_50th & snr_values >= snr_50th & max_dfof >= dfof_50th);
fprintf('  → Would select ~%d ROIs\n', moderate_count);

fprintf('\nLIBERAL (lower thresholds):\n');
fprintf('  MIN_EVENTS = %d\n', round(events_30th));
fprintf('  MIN_SNR = %.2f\n', snr_30th);
fprintf('  MIN_DFOF = %.4f\n', dfof_30th);
liberal_count = sum(events_per_roi >= events_30th & snr_values >= snr_30th & max_dfof >= dfof_30th);
fprintf('  → Would select ~%d ROIs\n', liberal_count);

fprintf('\n=== DISCUSSION QUESTIONS ===\n');
fprintf('1. Looking at the random ROIs, which ones show clear spontaneous release activity?\n');
fprintf('2. What level of "events" looks meaningful vs. noise?\n');
fprintf('3. What dF/F amplitude range represents real biological signals?\n');
fprintf('4. Should we prioritize many small events or fewer large events?\n');
fprintf('5. Which filtering level (conservative/moderate/liberal) seems appropriate?\n');

fprintf('\n=== NEXT STEPS ===\n');
fprintf('1. Review the plotted random ROIs\n');
fprintf('2. Identify which ones you consider "good" spontaneous release\n');
fprintf('3. Note their Events/SNR/dF/F values to guide filtering\n');
fprintf('4. Try one of the suggested criteria sets above\n');
fprintf('5. Iterate until you get ROIs that match your biological expectations\n');