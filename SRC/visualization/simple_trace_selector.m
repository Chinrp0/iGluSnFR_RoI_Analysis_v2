%% SIMPLE TRACE SELECTOR - Easy configuration for showing spontaneous release ROIs
% Just change the parameters below to get different ROI selections

clear; clc; close all;

%% === EASY CONFIGURATION ===
% Change these parameters to customize your trace selection

% How many traces to show?
NUM_TRACES = 16;  % Increase this to see more ROIs

% What kind of ROIs do you want to see?
SELECTION_MODE = 'mixed';  % Options: 'most_events', 'highest_snr', 'largest_amplitude', 'best_baseline', 'mixed'

% Minimum criteria (set lower to be more inclusive)
MIN_EVENTS = 2;        % Minimum number of detected events
MIN_SNR = 1.0;         % Minimum signal-to-noise ratio  
MIN_DFOF = 0.02;       % Minimum peak dF/F (2%)
MIN_BASELINE_QUALITY = 0.8;  % Minimum baseline quality (80%)

% Display options
SHOW_DFOF_PLOTS = true;   % Show dF/F traces in addition to raw traces
SHOW_DETAILED_INFO = true;  % Print detailed ROI information

fprintf('=== Simple Trace Selector for Spontaneous Release ===\n');
fprintf('Configuration:\n');
fprintf('  - Number of traces: %d\n', NUM_TRACES);
fprintf('  - Selection mode: %s\n', SELECTION_MODE);
fprintf('  - Min events: %d, Min SNR: %.1f, Min dF/F: %.3f\n', MIN_EVENTS, MIN_SNR, MIN_DFOF);

%% === Load Data (same as before) ===
addpath(genpath(pwd));
test_file = 'D:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglusnfr4f_NGR\Spont\GPU_SNR_Processed\5_raw_mean\CP_Snfr4-NGR_Doc2b-R213W_Cs1-c1_spont-01_mean.csv';

loader = csv_loader_v2();
[data, metadata] = loader.loadSingleFile(test_file);
config = tracenorm_config();

[baseline, outlier_mask, baseline_stats] = baseline_detector(data, config);
[dfof_data, dfof_stats] = dfof_calculator(data, baseline, config);

time_vector = (1:size(data, 1)) / config.frame_rate;

%% === Simple ROI Selection ===
% Apply basic filters first
events_per_roi = dfof_stats.events_per_roi;
snr_values = dfof_stats.snr;
max_dfof = max(dfof_data, [], 1);
baseline_quality = baseline_stats.valid_fraction;

% Find ROIs that meet minimum criteria
good_rois = find(events_per_roi >= MIN_EVENTS & ...
                 snr_values >= MIN_SNR & ...
                 max_dfof >= MIN_DFOF & ...
                 baseline_quality >= MIN_BASELINE_QUALITY);

fprintf('\nFiltering results:\n');
fprintf('  Total ROIs: %d\n', size(data, 2));
fprintf('  ROIs meeting minimum criteria: %d\n', length(good_rois));

if length(good_rois) < NUM_TRACES
    fprintf('  Warning: Only %d ROIs meet criteria, showing all of them\n', length(good_rois));
    NUM_TRACES = length(good_rois);
end

% Select specific ROIs based on mode
switch SELECTION_MODE
    case 'most_events'
        [~, sort_idx] = sort(events_per_roi(good_rois), 'descend');
        selected_rois = good_rois(sort_idx(1:NUM_TRACES));
        fprintf('  Selection: ROIs with most detected events\n');
        
    case 'highest_snr'
        [~, sort_idx] = sort(snr_values(good_rois), 'descend');
        selected_rois = good_rois(sort_idx(1:NUM_TRACES));
        fprintf('  Selection: ROIs with highest signal-to-noise ratio\n');
        
    case 'largest_amplitude'
        [~, sort_idx] = sort(max_dfof(good_rois), 'descend');
        selected_rois = good_rois(sort_idx(1:NUM_TRACES));
        fprintf('  Selection: ROIs with largest dF/F responses\n');
        
    case 'best_baseline'
        [~, sort_idx] = sort(baseline_quality(good_rois), 'descend');
        selected_rois = good_rois(sort_idx(1:NUM_TRACES));
        fprintf('  Selection: ROIs with most stable baselines\n');
        
    case 'mixed'
        % Take some from each category
        n_per_category = ceil(NUM_TRACES / 4);
        
        [~, idx1] = sort(events_per_roi(good_rois), 'descend');
        [~, idx2] = sort(snr_values(good_rois), 'descend');
        [~, idx3] = sort(max_dfof(good_rois), 'descend');
        [~, idx4] = sort(baseline_quality(good_rois), 'descend');
        
        mixed_selection = [good_rois(idx1(1:min(n_per_category, length(idx1)))), ...
                          good_rois(idx2(1:min(n_per_category, length(idx2)))), ...
                          good_rois(idx3(1:min(n_per_category, length(idx3)))), ...
                          good_rois(idx4(1:min(n_per_category, length(idx4))))];
        
        selected_rois = unique(mixed_selection);
        selected_rois = selected_rois(1:min(NUM_TRACES, length(selected_rois)));
        fprintf('  Selection: Mixed criteria (events + SNR + amplitude + baseline)\n');
        
    otherwise
        % Default: just take first N good ROIs
        selected_rois = good_rois(1:NUM_TRACES);
        fprintf('  Selection: First %d ROIs meeting criteria\n', NUM_TRACES);
end

%% === Create Raw + Baseline Plots ===
fprintf('\nCreating trace plots for %d selected ROIs...\n', length(selected_rois));

% Calculate subplot layout
rows = ceil(sqrt(length(selected_rois)));
cols = ceil(length(selected_rois) / rows);

fig1 = figure('Name', 'Raw Traces + Baseline', 'Position', [50, 50, 1400, 1000]);

for i = 1:length(selected_rois)
    roi_idx = selected_rois(i);
    
    subplot(rows, cols, i);
    
    % Plot raw data
    plot(time_vector, data(:, roi_idx), 'k-', 'LineWidth', 0.7, 'DisplayName', 'Raw');
    hold on;
    
    % Plot baseline
    plot(time_vector, baseline(:, roi_idx), 'r-', 'LineWidth', 1.2, 'DisplayName', 'Baseline');
    
    % Mark events
    events = outlier_mask(:, roi_idx);
    if any(events)
        scatter(time_vector(events), data(events, roi_idx), 20, 'o', ...
            'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor', [1 0.5 0], ...
            'DisplayName', sprintf('%d events', sum(events)));
    end
    
    xlabel('Time (s)');
    ylabel('Fluorescence');
    title(sprintf('ROI %d (%d events, SNR=%.1f)', roi_idx, events_per_roi(roi_idx), snr_values(roi_idx)));
    grid on;
    
    if i == 1
        legend('Location', 'best', 'FontSize', 8);
    end
end

sgtitle(sprintf('Raw Traces with Baseline - %s (%s selection)', metadata.filename, SELECTION_MODE), 'FontSize', 14);

%% === Create dF/F Plots (if requested) ===
if SHOW_DFOF_PLOTS
    fig2 = figure('Name', 'dF/F Traces', 'Position', [100, 100, 1400, 1000]);
    
    for i = 1:length(selected_rois)
        roi_idx = selected_rois(i);
        
        subplot(rows, cols, i);
        
        % Plot dF/F
        plot(time_vector, dfof_data(:, roi_idx), 'b-', 'LineWidth', 1.0);
        hold on;
        
        % Add reference lines
        yline(0, 'k--', 'Alpha', 0.3, 'Baseline');
        yline(0.05, 'g--', 'Alpha', 0.5, '5%');
        
        % Mark significant events (>5% dF/F)
        big_events = dfof_data(:, roi_idx) > 0.05;
        if any(big_events)
            scatter(time_vector(big_events), dfof_data(big_events, roi_idx), 15, ...
                'filled', 'MarkerFaceColor', [0 0.8 0]);
        end
        
        xlabel('Time (s)');
        ylabel('dF/F');
        title(sprintf('ROI %d (max dF/F = %.3f)', roi_idx, max_dfof(roi_idx)));
        grid on;
    end
    
    sgtitle('dF/F Normalized Traces', 'FontSize', 14);
end

%% === Print Detailed Information ===
if SHOW_DETAILED_INFO
    fprintf('\n=== Selected ROI Details ===\n');
    fprintf('ROI\tEvents\tSNR\tMax dF/F\tBaseline Quality\n');
    fprintf('---\t------\t---\t--------\t----------------\n');
    
    for i = 1:length(selected_rois)
        roi_idx = selected_rois(i);
        fprintf('%d\t%d\t%.1f\t%.3f\t\t%.1f%%\n', ...
            roi_idx, events_per_roi(roi_idx), snr_values(roi_idx), ...
            max_dfof(roi_idx), 100*baseline_quality(roi_idx));
    end
    
    fprintf('\nSummary statistics for selected ROIs:\n');
    fprintf('  Mean events per ROI: %.1f\n', mean(events_per_roi(selected_rois)));
    fprintf('  Mean SNR: %.2f\n', mean(snr_values(selected_rois)));
    fprintf('  Mean max dF/F: %.3f\n', mean(max_dfof(selected_rois)));
    fprintf('  Mean baseline quality: %.1f%%\n', 100*mean(baseline_quality(selected_rois)));
end

%% === Easy Access to Results ===
fprintf('\n=== Access Selected Data ===\n');
fprintf('Selected ROI indices: [%s]\n', sprintf('%d ', selected_rois));
fprintf('\nTo access data for these ROIs:\n');
fprintf('  selected_data = data(:, [%s]);\n', sprintf('%d ', selected_rois));
fprintf('  selected_baseline = baseline(:, [%s]);\n', sprintf('%s', sprintf('%d ', selected_rois)));
fprintf('  selected_dfof = dfof_data(:, [%s]);\n', sprintf('%s', sprintf('%d ', selected_rois)));

%% === Configuration Tips ===
fprintf('\n=== Tips for Different Selections ===\n');
fprintf('To see more active ROIs: Set SELECTION_MODE = ''most_events''\n');
fprintf('To see cleaner signals: Set SELECTION_MODE = ''highest_snr''\n');
fprintf('To see big responses: Set SELECTION_MODE = ''largest_amplitude''\n');
fprintf('To see stable baselines: Set SELECTION_MODE = ''best_baseline''\n');
fprintf('To see variety: Set SELECTION_MODE = ''mixed''\n');
fprintf('\nTo be more inclusive: Lower MIN_EVENTS, MIN_SNR, MIN_DFOF values\n');
fprintf('To be more selective: Raise these minimum values\n');