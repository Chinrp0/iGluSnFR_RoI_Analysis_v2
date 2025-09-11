%% QUICK EVENT VISUALIZATION - Run this to see clean dF/F traces with events
% This script shows exactly what you want: baseline-corrected traces with event markings

clear; clc; close all;

%% === Load and Process Data ===
fprintf('=== Quick Event Visualization ===\n');

% Your data file
test_file = 'D:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglusnfr4f_NGR\Spont\GPU_SNR_Processed\5_raw_mean\CP_Snfr4-NGR_Doc2b-R213W_Cs1-c1_spont-01_mean.csv';

% Add paths
addpath(genpath(pwd));

% Load data
fprintf('Loading data...\n');
loader = csv_loader_v2();
[data, metadata] = loader.loadSingleFile(test_file);
config = tracenorm_config();

% Process data
fprintf('Processing baseline and events...\n');
[baseline, outlier_mask, baseline_stats] = baseline_detector(data, config);
[dfof_data, dfof_stats] = dfof_calculator(data, baseline, config);
[event_mask, event_stats] = schmitt_event_detector(dfof_data, config);

fprintf('✓ Processed %d ROIs with %d total events\n', size(data, 2), event_stats.total_events);

%% === Create Focused Event Plots ===
fprintf('Creating focused event visualizations...\n');

% Call the focused plotting function
plot_handles = simple_event_trace_plotter(data, baseline, dfof_data, event_mask, ...
    event_stats, metadata, config);

%% === Quick Data Access ===
fprintf('\n=== Quick Results Summary ===\n');
fprintf('Total Events: %d\n', event_stats.total_events);
fprintf('Active ROIs: %d (%.1f%%)\n', event_stats.rois_with_events, ...
    100 * event_stats.rois_with_events / size(data, 2));
fprintf('Events per active ROI: %.1f\n', event_stats.total_events / event_stats.rois_with_events);
fprintf('Mean event duration: %.1f frames (%.0f ms)\n', ...
    event_stats.mean_event_duration, event_stats.mean_event_duration * 1000 / config.frame_rate);

% Show top 10 most active ROIs
[~, top_idx] = sort(event_stats.events_per_roi, 'descend');
fprintf('\nTop 10 most active ROIs:\n');
for i = 1:min(10, length(top_idx))
    roi_idx = top_idx(i);
    fprintf('  ROI %d: %d events\n', roi_idx, event_stats.events_per_roi(roi_idx));
end

%% === Access Event Data for Specific ROIs ===
fprintf('\n=== Event Data Access Example ===\n');

% Example: Get events for ROI 100
roi_of_interest = min(100, size(data, 2));
roi_events = event_mask(:, roi_of_interest);
roi_trace = dfof_data(:, roi_of_interest);

if any(roi_events)
    event_frames = find(roi_events);
    event_times = event_frames / config.frame_rate;
    
    fprintf('ROI %d has %d total event frames:\n', roi_of_interest, sum(roi_events));
    fprintf('  Event times: %.1f, %.1f, %.1f... seconds\n', event_times(1:min(3, length(event_times))));
    fprintf('  Max dF/F during events: %.3f\n', max(roi_trace(roi_events)));
else
    fprintf('ROI %d has no detected events\n', roi_of_interest);
end

%% === Save Results (Optional) ===
% Uncomment these lines if you want to save the event data
% fprintf('\nSaving event detection results...\n');
% save('event_detection_results.mat', 'event_mask', 'event_stats', 'dfof_data', 'metadata');
% fprintf('Saved to: event_detection_results.mat\n');

fprintf('\n✅ Event visualization complete!\n');
fprintf('Look at the 3 figures showing:\n');
fprintf('  1. Sample dF/F traces with events highlighted in red\n');
fprintf('  2. Event detection summary and statistics\n');
fprintf('  3. Most active ROIs with individual event episodes\n');