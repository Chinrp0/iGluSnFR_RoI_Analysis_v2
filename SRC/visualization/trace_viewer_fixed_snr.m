%% TRACE VIEWER WITH FIXED SNR - Now with realistic filtering!
% Uses corrected SNR calculation for proper ROI selection

clear; clc; close all;

%% === EASY CONFIGURATION ===
NUM_TRACES = 16;                    % How many traces to show
SELECTION_MODE = 'most_events';     % 'most_events', 'highest_snr', 'largest_amplitude', 'mixed'

% UPDATED filtering criteria based on corrected SNR
MIN_EVENTS = 10;                    % Minimum events (median was 39)
MIN_SNR = 4.0;                      % Corrected SNR threshold (mean was 4.9)  
MIN_DFOF = 0.030;                   % 3% minimum dF/F
MIN_BASELINE_QUALITY = 0.9;         % 90% baseline quality

fprintf('=== Trace Viewer with Fixed SNR ===\n');
fprintf('Configuration:\n');
fprintf('  - Number of traces: %d\n', NUM_TRACES);
fprintf('  - Selection mode: %s\n', SELECTION_MODE);
fprintf('  - Min events: %d, Min SNR: %.1f, Min dF/F: %.3f\n', MIN_EVENTS, MIN_SNR, MIN_DFOF);

%% === Load Data ===
addpath(genpath(pwd));
test_file = 'D:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglusnfr4f_NGR\Spont\GPU_SNR_Processed\5_raw_mean\CP_Snfr4-NGR_Doc2b-R213W_Cs1-c1_spont-01_mean.csv';

loader = csv_loader_v2();
[data, metadata] = loader.loadSingleFile(test_file);
config = tracenorm_config();

[baseline, outlier_mask, baseline_stats] = baseline_detector(data, config);
[dfof_data, dfof_stats] = dfof_calculator(data, baseline, config);

time_vector = (1:size(data, 1)) / config.frame_rate;

%% === Calculate Corrected SNR ===
fprintf('Calculating corrected SNR...\n');

% Use the best method from our testing (Signal amplitude / robust noise)
peak_response = max(dfof_data, [], 1);
median_dfof = median(dfof_data, 1);
signal_amplitude = peak_response - median_dfof;
noise_estimate = mad(dfof_data, 1, 1) * 1.4826;
corrected_snr = signal_amplitude ./ noise_estimate;

fprintf('✓ Corrected SNR calculated (range: %.1f - %.1f, mean: %.1f)\n', ...
    min(corrected_snr), max(corrected_snr), mean(corrected_snr));

%% === Apply Filtering with Corrected SNR ===
events_per_roi = dfof_stats.events_per_roi;
max_dfof = max(dfof_data, [], 1);
baseline_quality = baseline_stats.valid_fraction;

% Find ROIs that meet criteria  
good_rois = find(events_per_roi >= MIN_EVENTS & ...
                 corrected_snr >= MIN_SNR & ...
                 max_dfof >= MIN_DFOF & ...
                 baseline_quality >= MIN_BASELINE_QUALITY);

fprintf('\nFiltering results:\n');
fprintf('  Total ROIs: %d\n', size(data, 2));
fprintf('  ROIs meeting criteria: %d (%.1f%%)\n', length(good_rois), 100*length(good_rois)/size(data, 2));

if length(good_rois) < NUM_TRACES
    fprintf('  Warning: Only %d ROIs meet criteria, showing all of them\n', length(good_rois));
    NUM_TRACES = length(good_rois);
end

%% === Select Best ROIs Based on Mode ===
switch SELECTION_MODE
    case 'most_events'
        [~, sort_idx] = sort(events_per_roi(good_rois), 'descend');
        selected_rois = good_rois(sort_idx(1:NUM_TRACES));
        
    case 'highest_snr'
        [~, sort_idx] = sort(corrected_snr(good_rois), 'descend');
        selected_rois = good_rois(sort_idx(1:NUM_TRACES));
        
    case 'largest_amplitude'
        [~, sort_idx] = sort(max_dfof(good_rois), 'descend');
        selected_rois = good_rois(sort_idx(1:NUM_TRACES));
        
    case 'mixed'
        % Balanced selection
        n_per_cat = ceil(NUM_TRACES / 3);
        
        [~, events_idx] = sort(events_per_roi(good_rois), 'descend');
        [~, snr_idx] = sort(corrected_snr(good_rois), 'descend');
        [~, amp_idx] = sort(max_dfof(good_rois), 'descend');
        
        mixed_selection = [good_rois(events_idx(1:min(n_per_cat, length(events_idx)))), ...
                          good_rois(snr_idx(1:min(n_per_cat, length(snr_idx)))), ...
                          good_rois(amp_idx(1:min(n_per_cat, length(amp_idx))))];
        
        selected_rois = unique(mixed_selection);
        selected_rois = selected_rois(1:min(NUM_TRACES, length(selected_rois)));
        
    otherwise
        selected_rois = good_rois(1:NUM_TRACES);
end

fprintf('Selected %d ROIs using %s criteria\n', length(selected_rois), SELECTION_MODE);

%% === Create Raw + Baseline Traces ===
fprintf('\nCreating trace plots...\n');

rows = ceil(sqrt(length(selected_rois)));
cols = ceil(length(selected_rois) / rows);

fig1 = figure('Name', 'High-Quality Spontaneous Release ROIs', 'Position', [50, 50, 1600, 1200]);

for i = 1:length(selected_rois)
    roi_idx = selected_rois(i);
    
    subplot(rows, cols, i);
    
    % Plot raw data
    plot(time_vector, data(:, roi_idx), 'k-', 'LineWidth', 0.7);
    hold on;
    
    % Plot baseline
    plot(time_vector, baseline(:, roi_idx), 'r-', 'LineWidth', 1.2);
    
    % Mark events
    events = outlier_mask(:, roi_idx);
    if any(events)
        scatter(time_vector(events), data(events, roi_idx), 15, 'o', ...
            'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor', [1 0.5 0]);
    end
    
    xlabel('Time (s)');
    ylabel('Fluorescence');
    title(sprintf('ROI %d: %d events, SNR=%.1f', roi_idx, events_per_roi(roi_idx), corrected_snr(roi_idx)));
    grid on;
    
    if i == 1
        legend('Raw', 'Baseline', 'Events', 'Location', 'best', 'FontSize', 8);
    end
end

sgtitle(sprintf('High-Quality Spontaneous Release ROIs - %s', metadata.filename), 'FontSize', 14);

%% === Create dF/F Traces ===
fig2 = figure('Name', 'dF/F Traces - Selected ROIs', 'Position', [100, 100, 1600, 1200]);

for i = 1:length(selected_rois)
    roi_idx = selected_rois(i);
    
    subplot(rows, cols, i);
    
    % Plot dF/F
    plot(time_vector, dfof_data(:, roi_idx), 'b-', 'LineWidth', 1.0);
    hold on;
    
    % Reference lines
    yline(0, 'k--');
    yline(0.05, 'g:');
    yline(0.10, 'r:');
    
    % Mark significant events
    big_events = dfof_data(:, roi_idx) > 0.05;  % 5% threshold
    if any(big_events)
        scatter(time_vector(big_events), dfof_data(big_events, roi_idx), 12, ...
            'filled', 'MarkerFaceColor', [0 0.8 0]);
    end
    
    xlabel('Time (s)');
    ylabel('dF/F');
    title(sprintf('ROI %d: max=%.3f', roi_idx, max_dfof(roi_idx)));
    grid on;
end

sgtitle('dF/F Traces - High Quality ROIs', 'FontSize', 14);

%% === Display Detailed Statistics ===
fprintf('\n=== Selected ROI Statistics ===\n');
fprintf('ROI\tEvents\tSNR\tMax dF/F\tBaseline Quality\tSelection Rank\n');
fprintf('---\t------\t----\t--------\t----------------\t--------------\n');

for i = 1:length(selected_rois)
    roi_idx = selected_rois(i);
    fprintf('%d\t%d\t%.1f\t%.3f\t\t%.1f%%\t\t\t%d\n', ...
        roi_idx, events_per_roi(roi_idx), corrected_snr(roi_idx), ...
        max_dfof(roi_idx), 100*baseline_quality(roi_idx), i);
end

fprintf('\nSummary for selected ROIs:\n');
fprintf('  Mean events: %.1f\n', mean(events_per_roi(selected_rois)));
fprintf('  Mean SNR: %.1f\n', mean(corrected_snr(selected_rois)));
fprintf('  Mean max dF/F: %.3f\n', mean(max_dfof(selected_rois)));
fprintf('  Mean baseline quality: %.1f%%\n', 100*mean(baseline_quality(selected_rois)));

%% === Access Information ===
fprintf('\n=== Selected High-Quality ROIs ===\n');
fprintf('ROI indices: [%s]\n', sprintf('%d ', selected_rois));
fprintf('\nThese ROIs represent the best spontaneous release activity in your dataset!\n');
fprintf('\nTo access data:\n');
fprintf('  high_quality_data = data(:, [%s]);\n', sprintf('%d ', selected_rois));
fprintf('  high_quality_dfof = dfof_data(:, [%s]);\n', sprintf('%d ', selected_rois));

%% === Filtering Summary ===
fprintf('\n=== Filtering Effectiveness ===\n');
fprintf('Filter\t\t\tPassed\tPercent\n');
fprintf('------\t\t\t------\t-------\n');
fprintf('Events >= %d\t\t%d\t%.1f%%\n', MIN_EVENTS, sum(events_per_roi >= MIN_EVENTS), 100*mean(events_per_roi >= MIN_EVENTS));
fprintf('SNR >= %.1f\t\t\t%d\t%.1f%%\n', MIN_SNR, sum(corrected_snr >= MIN_SNR), 100*mean(corrected_snr >= MIN_SNR));
fprintf('dF/F >= %.3f\t\t%d\t%.1f%%\n', MIN_DFOF, sum(max_dfof >= MIN_DFOF), 100*mean(max_dfof >= MIN_DFOF));
fprintf('Baseline >= %.1f\t\t%d\t%.1f%%\n', MIN_BASELINE_QUALITY, sum(baseline_quality >= MIN_BASELINE_QUALITY), 100*mean(baseline_quality >= MIN_BASELINE_QUALITY));
fprintf('ALL criteria\t\t%d\t%.1f%%\n', length(good_rois), 100*length(good_rois)/size(data, 2));

fprintf('\n✅ Success! Now showing ROIs with genuine spontaneous release activity!\n');