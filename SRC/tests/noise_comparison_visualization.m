%% NOISE COMPARISON VISUALIZATION
% Compare highest and lowest noise ROIs using optimal detection
% Creates two figures: 8 highest noise ROIs and 8 lowest noise ROIs

clear; clc; close all;

fprintf('=== Optimal Event Detection with Noise Comparison ===\n');

%% === Setup and Data Loading ===
addpath(genpath(pwd));
test_file = 'D:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglusnfr4f_NGR\Spont\GPU_SNR_Processed\5_raw_mean\CP_Snfr4-NGR_Doc2b-R213W_Cs1-c1_spont-01_mean.csv';

% Load and process data
fprintf('Loading data and processing with optimal detection...\n');
loader = csv_loader_v2();
[data, metadata] = loader.loadSingleFile(test_file);
config = tracenorm_config();

% Enable verbose output for detailed timing
config.verbose = true;

% Baseline calculation
fprintf('\n1. Baseline calculation...\n');
[baseline, outlier_mask, baseline_stats] = baseline_detector(data, config);

% dF/F calculation  
fprintf('\n2. dF/F calculation...\n');
[dfof_data, dfof_stats] = dfof_calculator(data, baseline, config);

% Check parallel pool
fprintf('\n3. Parallel processing setup...\n');
parpool_obj = gcp('nocreate');
if isempty(parpool_obj)
    fprintf('   Starting parallel pool...\n');
    parpool('local', min(4, feature('numcores')));
else
    fprintf('   Using existing pool with %d workers\n', parpool_obj.NumWorkers);
end

% OPTIMAL EVENT DETECTION using all toolboxes
fprintf('\n4. Optimal event detection (Parallel + Signal Processing + Statistics)...\n');
tic;
[event_mask, event_stats] = optimal_schmitt_event_detector(dfof_data, config, baseline_stats);
optimal_time = toc;

fprintf('\n✅ Optimal detection completed in %.3f seconds\n', optimal_time);

%% === Noise Analysis and ROI Selection ===
noise_metrics = event_stats.noise_metrics;
noise_levels = noise_metrics.noise_std;
signal_quality = noise_metrics.signal_quality;

% Remove invalid ROIs
valid_rois = ~isnan(noise_levels) & noise_levels > 0 & ~isinf(noise_levels);
valid_indices = find(valid_rois);

fprintf('\nNoise analysis:\n');
fprintf('  Valid ROIs: %d/%d\n', sum(valid_rois), length(noise_levels));
fprintf('  Noise range: %.5f - %.4f dF/F\n', min(noise_levels(valid_rois)), max(noise_levels(valid_rois)));
fprintf('  Mean noise: %.5f ± %.5f dF/F\n', mean(noise_levels(valid_rois)), std(noise_levels(valid_rois)));

% Sort ROIs by noise level
[sorted_noise, sort_indices] = sort(noise_levels(valid_rois));
sorted_roi_indices = valid_indices(sort_indices);

% Select highest and lowest noise ROIs
num_to_show = 8;
highest_noise_rois = sorted_roi_indices((end-num_to_show+1):end);  % Top 8
lowest_noise_rois = sorted_roi_indices(1:num_to_show);            % Bottom 8

fprintf('\nSelected ROIs for comparison:\n');
fprintf('Highest noise ROIs: [%s]\n', sprintf('%d ', highest_noise_rois));
fprintf('  Noise levels: [%s]\n', sprintf('%.5f ', noise_levels(highest_noise_rois)));
fprintf('Lowest noise ROIs: [%s]\n', sprintf('%d ', lowest_noise_rois));
fprintf('  Noise levels: [%s]\n', sprintf('%.5f ', noise_levels(lowest_noise_rois)));

%% === Create Visualization Function ===
time_vector = (1:size(dfof_data, 1)) / config.frame_rate;

% Create highest noise ROI figure
create_noise_comparison_figure(dfof_data, event_mask, event_stats, time_vector, ...
    highest_noise_rois, noise_metrics, metadata, 'Highest Noise ROIs', 1);

% Create lowest noise ROI figure  
create_noise_comparison_figure(dfof_data, event_mask, event_stats, time_vector, ...
    lowest_noise_rois, noise_metrics, metadata, 'Lowest Noise ROIs', 2);

%% === Summary Statistics ===
fprintf('\n=== NOISE COMPARISON SUMMARY ===\n');

% Calculate statistics for each group
high_noise_stats = calculate_group_stats(highest_noise_rois, noise_metrics, event_stats, dfof_data);
low_noise_stats = calculate_group_stats(lowest_noise_rois, noise_metrics, event_stats, dfof_data);

fprintf('\nHIGHEST NOISE ROIs (n=%d):\n', length(highest_noise_rois));
print_group_stats(high_noise_stats);

fprintf('\nLOWEST NOISE ROIs (n=%d):\n', length(lowest_noise_rois));
print_group_stats(low_noise_stats);

fprintf('\nCOMPARISON:\n');
fprintf('  Noise ratio (high/low): %.1fx\n', high_noise_stats.mean_noise / low_noise_stats.mean_noise);
fprintf('  Events ratio (high/low): %.1fx\n', high_noise_stats.mean_events / max(low_noise_stats.mean_events, 1));
fprintf('  SNR ratio (high/low): %.1fx\n', high_noise_stats.mean_snr / max(low_noise_stats.mean_snr, 0.01));

%% === Method Effectiveness Analysis ===
fprintf('\n=== NOISE ESTIMATION METHOD EFFECTIVENESS ===\n');

methods_used = noise_metrics.methods_used;
fprintf('Methods used across all ROIs: %s\n', strjoin(methods_used, ', '));

for method = methods_used
    method_rois = strcmp(noise_metrics.noise_method, method);
    method_count = sum(method_rois);
    method_noise = mean(noise_metrics.noise_std(method_rois));
    method_snr = mean(noise_metrics.signal_quality(method_rois));
    
    fprintf('  %s: %d ROIs, avg noise=%.5f, avg SNR=%.1f\n', ...
        method{1}, method_count, method_noise, method_snr);
end

if isfield(noise_metrics, 'convergence_iterations')
    convergence_stats = noise_metrics.convergence_iterations(noise_metrics.convergence_iterations > 0);
    if ~isempty(convergence_stats)
        fprintf('  Iterative convergence: %.1f iterations avg (range: %d-%d)\n', ...
            mean(convergence_stats), min(convergence_stats), max(convergence_stats));
    end
end

fprintf('\n✅ Visualization complete! Check the two figures showing noise comparison.\n');

%% === HELPER FUNCTIONS ===

function create_noise_comparison_figure(dfof_data, event_mask, event_stats, time_vector, ...
    roi_indices, noise_metrics, metadata, fig_title, fig_number)
    % Create 2x4 subplot figure comparing ROIs
    
    fig = figure(fig_number);
    fig.Name = sprintf('%s - %s', fig_title, metadata.filename);
    fig.Position = [50 + (fig_number-1)*50, 100, 1200, 1000];
    
    % 2x4 layout as requested
    rows = 4;
    cols = 2;
    
    for i = 1:length(roi_indices)
        roi_idx = roi_indices(i);
        
        subplot(rows, cols, i);
        
        % Plot dF/F trace
        trace = dfof_data(:, roi_idx);
        plot(time_vector, trace, 'k-', 'LineWidth', 1.0);
        hold on;
        
        % Highlight detected events
        if size(event_mask, 2) >= roi_idx
            events = event_mask(:, roi_idx);
            if any(events)
                % Highlight event periods in red
                event_trace = trace;
                event_trace(~events) = NaN;
                plot(time_vector, event_trace, 'r-', 'LineWidth', 2.5);
                
                % Add peak markers above events
                event_peaks = find_event_peaks_in_events(trace, events);
                if ~isempty(event_peaks.frames)
                    marker_height = event_peaks.amplitudes + 0.03 * range(trace, 'omitnan');
                    scatter(time_vector(event_peaks.frames), marker_height, 60, '^', ...
                        'MarkerFaceColor', [0 0.8 0], 'MarkerEdgeColor', [0 0.6 0], ...
                        'LineWidth', 1.5);
                end
            end
        end
        
        % Add Schmitt trigger thresholds
        if length(noise_metrics.upper_thresholds) >= roi_idx
            upper_thresh = noise_metrics.upper_thresholds(roi_idx);
            lower_thresh = noise_metrics.lower_thresholds(roi_idx);
            
            yline(upper_thresh, 'g--', 'LineWidth', 1.2, 'Alpha', 0.7);
            yline(lower_thresh, 'c--', 'LineWidth', 1.2, 'Alpha', 0.7);
        end
        
        % Add zero reference
        yline(0, 'k:', 'Alpha', 0.5, 'LineWidth', 0.8);
        
        % Get metrics for this ROI
        num_events = event_stats.events_per_roi(roi_idx);
        noise_level = noise_metrics.noise_std(roi_idx);
        snr = noise_metrics.signal_quality(roi_idx);
        method_used = noise_metrics.noise_method{roi_idx};
        max_dfof = max(trace, [], 'omitnan');
        
        xlabel('Time (s)', 'FontSize', 10);
        ylabel('dF/F', 'FontSize', 10);
        title(sprintf('ROI %d: %d events, σ=%.5f\nSNR=%.1f, max=%.3f (%s)', ...
            roi_idx, num_events, noise_level, snr, max_dfof, method_used), ...
            'FontSize', 10, 'FontWeight', 'bold');
        grid on; grid minor;
        
        % Set consistent y-limits with room for markers
        y_range = [min(trace, [], 'omitnan'), max(trace, [], 'omitnan')];
        if diff(y_range) > 0
            y_margin = diff(y_range) * 0.25;
            ylim([y_range(1) - y_margin, y_range(2) + y_margin]);
        end
        
        % Add legend only to first subplot
        if i == 1 && any(event_mask(:, roi_idx))
            legend('dF/F', 'Events', 'Peak Markers', 'Upper (3.0σ)', 'Lower (1.5σ)', ...
                'Location', 'best', 'FontSize', 8);
        elseif i == 1
            legend('dF/F', 'Upper (3.0σ)', 'Lower (1.5σ)', 'Location', 'best', 'FontSize', 8);
        end
    end
    
    % Calculate summary for this group
    group_noise = noise_metrics.noise_std(roi_indices);
    group_events = event_stats.events_per_roi(roi_indices);
    group_snr = noise_metrics.signal_quality(roi_indices);
    
    sgtitle(sprintf('%s - Optimal Detection Results\n%s\nNoise: %.5f±%.5f, Events: %.1f±%.1f, SNR: %.1f±%.1f', ...
        fig_title, metadata.filename, ...
        mean(group_noise), std(group_noise), ...
        mean(group_events), std(group_events), ...
        mean(group_snr), std(group_snr)), 'FontSize', 14, 'FontWeight', 'bold');
end

function event_peaks = find_event_peaks_in_events(trace, event_mask)
    % Find peak locations within event periods
    
    event_peaks = struct('frames', [], 'amplitudes', []);
    
    if ~any(event_mask)
        return;
    end
    
    % Find discrete event episodes
    event_starts = find(diff([false; event_mask]) == 1);
    event_ends = find(diff([event_mask; false]) == -1);
    
    peak_frames = [];
    peak_amplitudes = [];
    
    % For each event episode, find the peak
    for i = 1:length(event_starts)
        start_frame = event_starts(i);
        end_frame = event_ends(i);
        
        % Find peak within this event
        event_segment = trace(start_frame:end_frame);
        [max_val, max_idx] = max(event_segment);
        
        if ~isnan(max_val)
            % Convert back to original frame index
            peak_frame = start_frame + max_idx - 1;
            
            peak_frames(end+1) = peak_frame;
            peak_amplitudes(end+1) = max_val;
        end
    end
    
    event_peaks.frames = peak_frames;
    event_peaks.amplitudes = peak_amplitudes;
end

function stats = calculate_group_stats(roi_indices, noise_metrics, event_stats, dfof_data)
    % Calculate statistics for a group of ROIs
    
    stats = struct();
    
    % Noise statistics
    group_noise = noise_metrics.noise_std(roi_indices);
    stats.mean_noise = mean(group_noise);
    stats.std_noise = std(group_noise);
    stats.min_noise = min(group_noise);
    stats.max_noise = max(group_noise);
    
    % Event statistics
    group_events = event_stats.events_per_roi(roi_indices);
    stats.mean_events = mean(group_events);
    stats.std_events = std(group_events);
    stats.total_events = sum(group_events);
    stats.active_rois = sum(group_events > 0);
    
    % Signal quality
    group_snr = noise_metrics.signal_quality(roi_indices);
    stats.mean_snr = mean(group_snr);
    stats.std_snr = std(group_snr);
    
    % dF/F range
    group_max_dfof = max(dfof_data(:, roi_indices), [], 1);
    stats.mean_max_dfof = mean(group_max_dfof);
    stats.std_max_dfof = std(group_max_dfof);
    
    % Methods used
    group_methods = noise_metrics.noise_method(roi_indices);
    stats.methods_used = unique(group_methods);
    stats.method_counts = arrayfun(@(m) sum(strcmp(group_methods, m)), stats.methods_used);
end

function print_group_stats(stats)
    % Print formatted statistics for a group
    
    fprintf('  Noise: %.5f ± %.5f (range: %.5f - %.5f)\n', ...
        stats.mean_noise, stats.std_noise, stats.min_noise, stats.max_noise);
    fprintf('  Events: %.1f ± %.1f total (%d total, %d/%d active)\n', ...
        stats.mean_events, stats.std_events, stats.total_events, stats.active_rois, length(stats.active_rois));
    fprintf('  SNR: %.1f ± %.1f\n', stats.mean_snr, stats.std_snr);
    fprintf('  Max dF/F: %.3f ± %.3f\n', stats.mean_max_dfof, stats.std_max_dfof);
    fprintf('  Methods: %s\n', strjoin(stats.methods_used, ', '));
end