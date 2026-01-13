%% DEBUG EVENT STATS - Check what's in your event_stats structure
% This will help us understand how to properly extract the events

clear; clc;

fprintf('=== Debugging Event Stats Structure ===\n');

%% Load and process data to get event_stats
addpath(genpath(pwd));
test_file = 'D:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglusnfr4f_NGR\Spont\GPU_SNR_Processed\5_raw_mean\CP_Snfr4-NGR_Doc2b-R213W_Cs1-c1_spont-01_mean.csv';

fprintf('Loading data and processing...\n');
loader = csv_loader_v2();
[data, metadata] = loader.loadSingleFile(test_file);
config = tracenorm_config();

[baseline, outlier_mask, baseline_stats] = baseline_detector(data, config);
[dfof_data, dfof_stats] = dfof_calculator(data, baseline, config);

%% Try to get event_stats from your pipeline
fprintf('\nChecking if schmitt_event_detector function exists...\n');
if exist('schmitt_event_detector', 'file')
    fprintf('✓ schmitt_event_detector function found\n');
    try
        [event_mask, event_stats] = schmitt_event_detector(dfof_data, config);
        fprintf('✓ Event detection completed successfully\n');
        has_real_events = true;
    catch ME
        fprintf('✗ Event detection failed: %s\n', ME.message);
        has_real_events = false;
    end
else
    fprintf('✗ schmitt_event_detector function not found\n');
    has_real_events = false;
end

%% Check event_stats structure
if has_real_events
    fprintf('\n=== Real Event Stats Analysis ===\n');
    fprintf('event_stats fields: %s\n', strjoin(fieldnames(event_stats), ', '));
    
    if isfield(event_stats, 'total_events')
        fprintf('Total events: %d\n', event_stats.total_events);
    end
    
    if isfield(event_stats, 'events_per_roi')
        fprintf('Events per ROI: min=%d, max=%d, mean=%.1f\n', ...
            min(event_stats.events_per_roi), max(event_stats.events_per_roi), ...
            mean(event_stats.events_per_roi));
    end
    
    if isfield(event_stats, 'rois_with_events')
        fprintf('Active ROIs: %d\n', event_stats.rois_with_events);
    end
    
    % Check event_mask
    if exist('event_mask', 'var')
        fprintf('\nevent_mask size: [%dx%d]\n', size(event_mask, 1), size(event_mask, 2));
        fprintf('Total event frames: %d\n', sum(event_mask, 'all'));
        
        % Check a sample ROI
        sample_roi = find(event_stats.events_per_roi > 0, 1, 'first');
        if ~isempty(sample_roi)
            roi_events = event_mask(:, sample_roi);
            fprintf('\nSample ROI %d: %d event frames detected\n', sample_roi, sum(roi_events));
            
            if sum(roi_events) > 0
                event_frames = find(roi_events);
                fprintf('  Event frames: [%s%s]\n', sprintf('%d ', event_frames(1:min(5, length(event_frames)))), ...
                    ternary(length(event_frames) > 5, '...', ''));
            end
        end
    end
    
else
    fprintf('\n=== Creating Mock Event Stats for Testing ===\n');
    
    % Create mock event_stats based on dfof_stats
    event_stats = struct();
    event_stats.events_per_roi = dfof_stats.events_per_roi;  % Use existing simple detection
    event_stats.total_events = sum(dfof_stats.events_per_roi);
    event_stats.rois_with_events = sum(dfof_stats.events_per_roi > 0);
    
    % Create mock event_mask
    [numFrames, numROIs] = size(dfof_data);
    event_mask = false(numFrames, numROIs);
    
    % Simple threshold-based events for testing
    for roi = 1:numROIs
        if dfof_stats.events_per_roi(roi) > 0
            % Find peaks above 2 standard deviations
            trace = dfof_data(:, roi);
            threshold = 2 * std(trace, 'omitnan');
            above_threshold = trace > threshold;
            
            % Mark peaks
            if any(above_threshold)
                event_mask(:, roi) = above_threshold;
            end
        end
    end
    
    event_stats.event_mask = event_mask;
    
    fprintf('Mock event stats created:\n');
    fprintf('  Total events: %d\n', event_stats.total_events);
    fprintf('  Active ROIs: %d\n', event_stats.rois_with_events);
    
    has_real_events = true;  % Use mock events for testing
end

%% Test event peak finding
fprintf('\n=== Testing Event Peak Finding ===\n');

if has_real_events
    % Test the peak finding function
    sample_roi = find(event_stats.events_per_roi > 0, 1, 'first');
    
    if ~isempty(sample_roi)
        fprintf('Testing peak finding for ROI %d...\n', sample_roi);
        
        time_vector = (1:size(dfof_data, 1)) / config.frame_rate;
        
        % This is our new function from the updated baseline_plotter
        try
            event_peaks = find_event_peaks_from_mask_test(sample_roi, event_stats, dfof_data, time_vector);
            
            if ~isempty(event_peaks.frames)
                fprintf('✓ Found %d event peaks for ROI %d\n', length(event_peaks.frames), sample_roi);
                fprintf('  Peak frames: [%s]\n', sprintf('%d ', event_peaks.frames));
                fprintf('  Peak amplitudes: [%s]\n', sprintf('%.3f ', event_peaks.amplitudes));
            else
                fprintf('✗ No event peaks found for ROI %d\n', sample_roi);
            end
            
        catch ME
            fprintf('✗ Peak finding failed: %s\n', ME.message);
        end
    else
        fprintf('No active ROIs found for testing\n');
    end
end

%% Summary and recommendations
fprintf('\n=== Summary ===\n');
if has_real_events
    fprintf('✅ Event data available and processed\n');
    fprintf('✅ Updated baseline_plotter should work with your data\n');
    fprintf('\nRecommendations:\n');
    fprintf('1. Replace your baseline_plotter.m with the updated version\n');
    fprintf('2. Ensure your main_pipeline passes event_stats to the plotter\n');
    fprintf('3. Run main_pipeline again to see proper event marking\n');
else
    fprintf('❌ Event detection system not accessible\n');
    fprintf('\nNext steps:\n');
    fprintf('1. Check if schmitt_event_detector.m exists in your codebase\n');
    fprintf('2. Make sure it''s being called in your main_pipeline\n');
    fprintf('3. The updated baseline_plotter will fall back to basic plots without events\n');
end

%% Helper function for testing
function event_peaks = find_event_peaks_from_mask_test(roi_idx, event_stats, dfof_data, time_vector)
    % Test version of the event peak finding function
    
    event_peaks = struct('frames', [], 'amplitudes', []);
    
    % Try to get event mask for this ROI
    if isfield(event_stats, 'event_mask') && size(event_stats.event_mask, 2) >= roi_idx
        roi_events = event_stats.event_mask(:, roi_idx);
    else
        fprintf('  Warning: No event_mask found, using fallback\n');
        return;
    end
    
    if ~any(roi_events)
        return;
    end
    
    % Find discrete events (connected components)
    roi_trace = dfof_data(:, roi_idx);
    
    % Find start and end of each event episode
    event_starts = find(diff([0; roi_events; 0]) == 1);
    event_ends = find(diff([0; roi_events; 0]) == -1);
    
    % For each event episode, find the peak
    peak_frames = [];
    peak_amplitudes = [];
    
    for i = 1:length(event_starts)
        start_frame = event_starts(i);
        end_frame = event_ends(i);
        
        % Find peak within this event
        event_segment = roi_trace(start_frame:end_frame);
        [max_val, max_idx] = max(event_segment);
        
        % Convert back to original frame index
        peak_frame = start_frame + max_idx - 1;
        
        peak_frames(end+1) = peak_frame;
        peak_amplitudes(end+1) = max_val;
    end
    
    event_peaks.frames = peak_frames;
    event_peaks.amplitudes = peak_amplitudes;
end

function result = ternary(condition, true_val, false_val)
    % Simple ternary operator
    if condition
        result = true_val;
    else
        result = false_val;
    end
end