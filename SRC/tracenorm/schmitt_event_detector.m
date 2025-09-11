function [event_mask, event_stats] = schmitt_event_detector(dfof_data, config)
    % SCHMITT_EVENT_DETECTOR - Detect events using Schmitt trigger thresholding
    % Implements hysteresis-based event detection with sensor kinetics extension
    %
    % Inputs:
    %   dfof_data - [frames x ROIs] normalized dF/F traces
    %   config    - Configuration struct from tracenorm_config()
    %
    % Outputs:
    %   event_mask  - [frames x ROIs] logical mask of detected events
    %   event_stats - Struct with event detection statistics
    
    if nargin < 2
        config = tracenorm_config();
    end
    
    [numFrames, numROIs] = size(dfof_data);
    
    % Get Schmitt trigger parameters from config
    if ~isfield(config, 'event_detection')
        % Default parameters
        upper_threshold_sigma = 3.0;    % Upper threshold (event start)
        lower_threshold_sigma = 1.5;    % Lower threshold (event end)
        decay_extension_frames = 3;     % Extend events by decay time (75ms = 3 frames)
        rise_extension_frames = 1;      % Extend events before start (conservative)
    else
        upper_threshold_sigma = config.event_detection.upper_threshold_sigma;
        lower_threshold_sigma = config.event_detection.lower_threshold_sigma;
        decay_extension_frames = config.event_detection.decay_extension_frames;
        rise_extension_frames = config.event_detection.rise_extension_frames;
    end
    
    if config.verbose
        fprintf('Detecting events using Schmitt trigger (%.1fσ/%.1fσ thresholds)...\n', ...
            upper_threshold_sigma, lower_threshold_sigma);
    end
    
    %% === Calculate Thresholds Per ROI ===
    % Use robust noise estimation
    noise_std = mad(dfof_data, 1, 1) * 1.4826;  % Convert MAD to std estimate
    
    % Calculate thresholds
    upper_threshold = upper_threshold_sigma * noise_std;  % [1 x ROIs]
    lower_threshold = lower_threshold_sigma * noise_std;  % [1 x ROIs]
    
    %% === Apply Schmitt Trigger Detection ===
    event_mask = false(numFrames, numROIs);
    
    % Per-ROI event detection
    for roi = 1:numROIs
        roi_trace = dfof_data(:, roi);
        roi_upper = upper_threshold(roi);
        roi_lower = lower_threshold(roi);
        
        % Skip ROIs with insufficient data or too much noise
        if all(isnan(roi_trace)) || roi_upper <= 0
            continue;
        end
        
        % Apply Schmitt trigger logic
        roi_events = apply_schmitt_trigger(roi_trace, roi_upper, roi_lower);
        
        % Extend events based on sensor kinetics
        roi_events = extend_events_for_kinetics(roi_events, rise_extension_frames, decay_extension_frames);
        
        event_mask(:, roi) = roi_events;
    end
    
    %% === Calculate Event Statistics ===
    event_stats = calculate_event_statistics(event_mask, dfof_data, upper_threshold, lower_threshold, config);
    
    if config.verbose
        fprintf('  Detected %d total events across %d ROIs\n', ...
            event_stats.total_events, event_stats.rois_with_events);
    end
end

function event_mask = apply_schmitt_trigger(trace, upper_thresh, lower_thresh)
    % Apply Schmitt trigger logic to a single trace
    
    numFrames = length(trace);
    event_mask = false(numFrames, 1);
    
    % State machine: 'baseline' or 'in_event'
    state = 'baseline';
    event_start = 0;
    
    for frame = 1:numFrames
        signal = trace(frame);
        
        % Skip NaN values
        if isnan(signal)
            continue;
        end
        
        switch state
            case 'baseline'
                % Look for signal crossing above upper threshold
                if signal > upper_thresh
                    state = 'in_event';
                    event_start = frame;
                end
                
            case 'in_event'
                % Look for signal falling below lower threshold
                if signal < lower_thresh
                    % Event ends - mark all frames from start to here
                    event_mask(event_start:frame) = true;
                    state = 'baseline';
                    event_start = 0;
                end
        end
    end
    
    % Handle case where event extends to end of trace
    if strcmp(state, 'in_event') && event_start > 0
        event_mask(event_start:end) = true;
    end
end

function extended_mask = extend_events_for_kinetics(event_mask, rise_frames, decay_frames)
    % Extend event boundaries based on sensor rise/decay kinetics
    
    numFrames = length(event_mask);
    extended_mask = event_mask;
    
    if ~any(event_mask)
        return;  % No events to extend
    end
    
    % Find event boundaries
    event_starts = find(diff([false; event_mask]) == 1);
    event_ends = find(diff([event_mask; false]) == -1);
    
    % Extend each event
    for i = 1:length(event_starts)
        start_frame = event_starts(i);
        end_frame = event_ends(i);
        
        % Extend backward for rise time (conservative since we can't measure rise)
        extended_start = max(1, start_frame - rise_frames);
        
        % Extend forward for decay time (75ms = 3 frames)
        extended_end = min(numFrames, end_frame + decay_frames);
        
        % Mark extended region
        extended_mask(extended_start:extended_end) = true;
    end
end

function stats = calculate_event_statistics(event_mask, dfof_data, upper_threshold, lower_threshold, config)
    % Calculate comprehensive event detection statistics
    
    [numFrames, numROIs] = size(event_mask);
    
    %% === Basic Event Counts ===
    events_per_roi = sum(event_mask, 1);  % Total event frames per ROI
    rois_with_events = sum(events_per_roi > 0);
    total_event_frames = sum(event_mask, 'all');
    
    %% === Event Characteristics ===
    event_durations = [];
    event_amplitudes = [];
    event_count_per_roi = zeros(1, numROIs);
    
    for roi = 1:numROIs
        roi_events = event_mask(:, roi);
        roi_dfof = dfof_data(:, roi);
        
        if any(roi_events)
            % Find discrete events (groups of consecutive true values)
            event_starts = find(diff([false; roi_events]) == 1);
            event_ends = find(diff([roi_events; false]) == -1);
            
            event_count_per_roi(roi) = length(event_starts);
            
            % Calculate duration and amplitude for each event
            for e = 1:length(event_starts)
                start_frame = event_starts(e);
                end_frame = event_ends(e);
                
                % Event duration (frames)
                duration = end_frame - start_frame + 1;
                event_durations(end+1) = duration;
                
                % Event amplitude (peak dF/F during event)
                event_trace = roi_dfof(start_frame:end_frame);
                amplitude = max(event_trace, [], 'omitnan');
                if ~isnan(amplitude)
                    event_amplitudes(end+1) = amplitude;
                end
            end
        end
    end
    
    %% === Threshold Effectiveness ===
    % Check what fraction of detected events actually exceeded thresholds
    peak_amplitudes_per_roi = zeros(1, numROIs);
    for roi = 1:numROIs
        if events_per_roi(roi) > 0
            roi_events = event_mask(:, roi);
            peak_amplitudes_per_roi(roi) = max(dfof_data(roi_events, roi), [], 'omitnan');
        end
    end
    
    %% === Compile Statistics ===
    stats = struct();
    
    % Basic counts
    stats.total_events = sum(event_count_per_roi);  % Number of discrete events
    stats.total_event_frames = total_event_frames;   % Total frames marked as events
    stats.rois_with_events = rois_with_events;
    stats.events_per_roi = event_count_per_roi;      % Discrete events per ROI
    stats.event_frames_per_roi = events_per_roi;     % Event frames per ROI
    
    % Event characteristics
    if ~isempty(event_durations)
        stats.mean_event_duration = mean(event_durations);
        stats.median_event_duration = median(event_durations);
        stats.event_duration_range = [min(event_durations), max(event_durations)];
    else
        stats.mean_event_duration = 0;
        stats.median_event_duration = 0;
        stats.event_duration_range = [0, 0];
    end
    
    if ~isempty(event_amplitudes)
        stats.mean_event_amplitude = mean(event_amplitudes);
        stats.median_event_amplitude = median(event_amplitudes);
        stats.event_amplitude_range = [min(event_amplitudes), max(event_amplitudes)];
    else
        stats.mean_event_amplitude = 0;
        stats.median_event_amplitude = 0;
        stats.event_amplitude_range = [0, 0];
    end
    
    % Threshold information
    stats.upper_thresholds = upper_threshold;
    stats.lower_thresholds = lower_threshold;
    stats.mean_upper_threshold = mean(upper_threshold, 'omitnan');
    stats.mean_lower_threshold = mean(lower_threshold, 'omitnan');
    
    % Detection efficiency
    stats.fraction_rois_with_events = rois_with_events / numROIs;
    stats.fraction_frames_in_events = total_event_frames / numel(event_mask);
    
    % Summary for integration with existing code
    stats.event_summary = struct();
    stats.event_summary.rois_with_events = rois_with_events;
    stats.event_summary.total_events = stats.total_events;
    stats.event_summary.fraction_active_rois = stats.fraction_rois_with_events;
    if rois_with_events > 0
        stats.event_summary.mean_events_per_active_roi = mean(event_count_per_roi(event_count_per_roi > 0));
    else
        stats.event_summary.mean_events_per_active_roi = 0;
    end
end