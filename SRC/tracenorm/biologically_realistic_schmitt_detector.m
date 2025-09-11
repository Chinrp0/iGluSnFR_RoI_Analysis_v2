function [event_mask, event_stats] = biologically_realistic_schmitt_detector(dfof_data, config, baseline_stats)
    % BIOLOGICALLY_REALISTIC_SCHMITT_DETECTOR 
    % Fixes the issue where events are too short by requiring minimum duration
    % and using more appropriate threshold spacing
    
    if nargin < 2
        config = tracenorm_config();
    end
    
    if nargin < 3
        baseline_stats = [];
    end
    
    [numFrames, numROIs] = size(dfof_data);
    
    % Enhanced parameters for biological realism
    upper_threshold_sigma = config.event_detection.upper_threshold_sigma;
    lower_threshold_sigma = config.event_detection.lower_threshold_sigma;
    
    % NEW: Minimum event duration requirements
    min_event_duration = 3;  % At least 3 frames (75ms at 40Hz) for real events
    max_gap_to_merge = 2;    % Merge events separated by ≤2 frames
    decay_extension_frames = config.event_detection.decay_extension_frames;
    rise_extension_frames = config.event_detection.rise_extension_frames;
    
    if config.verbose
        fprintf('Biologically realistic Schmitt trigger detection...\n');
        fprintf('  Thresholds: %.1fσ/%.1fσ\n', upper_threshold_sigma, lower_threshold_sigma);
        fprintf('  Min event duration: %d frames (%.0f ms)\n', ...
            min_event_duration, min_event_duration * 1000 / config.frame_rate);
    end
    
    %% === Calculate Adaptive Thresholds ===
    [upper_threshold, lower_threshold, noise_metrics] = calculate_adaptive_thresholds(...
        dfof_data, upper_threshold_sigma, lower_threshold_sigma, config);
    
    %% === Apply Enhanced Schmitt Trigger ===
    event_mask = false(numFrames, numROIs);
    
    for roi = 1:numROIs
        roi_trace = dfof_data(:, roi);
        roi_upper = upper_threshold(roi);
        roi_lower = lower_threshold(roi);
        
        if all(isnan(roi_trace)) || roi_upper <= 0 || isnan(roi_upper)
            continue;
        end
        
        % Apply enhanced Schmitt trigger with duration requirements
        roi_events = enhanced_schmitt_trigger_with_duration(roi_trace, roi_upper, roi_lower, ...
            min_event_duration, max_gap_to_merge);
        
        % Apply kinetics extension (now just fine-tuning, not main duration)
        roi_events = extend_events_for_kinetics(roi_events, rise_extension_frames, decay_extension_frames);
        
        event_mask(:, roi) = roi_events;
    end
    
    %% === Calculate Event Statistics ===
    event_stats = calculate_enhanced_event_statistics(event_mask, dfof_data, upper_threshold, ...
        lower_threshold, noise_metrics, config);
    
    event_stats.event_mask = event_mask;
    event_stats.noise_metrics = noise_metrics;
    
    if config.verbose
        fprintf('  Detected %d total events across %d ROIs\n', ...
            event_stats.total_events, event_stats.rois_with_events);
        fprintf('  Mean event duration: %.1f frames (%.0f ms)\n', ...
            event_stats.mean_event_duration, event_stats.mean_event_duration * 1000 / config.frame_rate);
    end
end

function [upper_threshold, lower_threshold, noise_metrics] = calculate_adaptive_thresholds(...
    dfof_data, upper_sigma, lower_sigma, config)
    % Calculate thresholds with better noise estimation and adaptive spacing
    
    [~, numROIs] = size(dfof_data);
    
    upper_threshold = zeros(1, numROIs);
    lower_threshold = zeros(1, numROIs);
    noise_std = zeros(1, numROIs);
    
    for roi = 1:numROIs
        roi_trace = dfof_data(:, roi);
        
        if all(isnan(roi_trace))
            noise_std(roi) = 0.01;
            upper_threshold(roi) = upper_sigma * 0.01;
            lower_threshold(roi) = lower_sigma * 0.01;
            continue;
        end
        
        % Use percentile-based noise estimation (excludes events)
        baseline_percentile = 70;
        threshold_70 = prctile(roi_trace, baseline_percentile);
        baseline_periods = roi_trace <= threshold_70;
        
        if sum(baseline_periods) > 20
            baseline_data = roi_trace(baseline_periods);
            noise_estimate = mad(baseline_data, 1) * 1.4826;
        else
            noise_estimate = mad(roi_trace, 1) * 1.4826;
        end
        
        noise_std(roi) = noise_estimate;
        
        % ENHANCED: Adaptive threshold spacing based on signal characteristics
        signal_range = range(roi_trace, 'omitnan');
        signal_to_noise = signal_range / noise_estimate;
        
        % For very noisy ROIs, increase threshold separation
        if signal_to_noise < 5  % Low SNR ROI
            adjusted_lower_sigma = lower_sigma * 0.8;  % Lower the lower threshold
        elseif signal_to_noise > 20  % High SNR ROI
            adjusted_lower_sigma = lower_sigma * 1.2;  % Raise the lower threshold
        else
            adjusted_lower_sigma = lower_sigma;
        end
        
        upper_threshold(roi) = upper_sigma * noise_estimate;
        lower_threshold(roi) = adjusted_lower_sigma * noise_estimate;
        
        % Ensure minimum threshold separation for hysteresis
        min_separation = 2 * noise_estimate;  % At least 2σ separation
        if (upper_threshold(roi) - lower_threshold(roi)) < min_separation
            lower_threshold(roi) = upper_threshold(roi) - min_separation;
            lower_threshold(roi) = max(lower_threshold(roi), 0.5 * noise_estimate);  % Don't go too low
        end
    end
    
    % Compile noise metrics
    noise_metrics = struct();
    noise_metrics.noise_std = noise_std;
    noise_metrics.upper_thresholds = upper_threshold;
    noise_metrics.lower_thresholds = lower_threshold;
    noise_metrics.threshold_separation = upper_threshold - lower_threshold;
end

function event_mask = enhanced_schmitt_trigger_with_duration(trace, upper_thresh, lower_thresh, ...
    min_duration, max_gap_to_merge)
    % Enhanced Schmitt trigger that enforces minimum duration and merges nearby events
    
    numFrames = length(trace);
    event_mask = false(numFrames, 1);
    
    %% === Phase 1: Basic Schmitt Trigger Detection ===
    state = 'baseline';
    event_start = 0;
    potential_events = [];  % Store [start, end] pairs
    
    for frame = 1:numFrames
        signal = trace(frame);
        
        if isnan(signal)
            continue;
        end
        
        switch state
            case 'baseline'
                if signal > upper_thresh
                    state = 'in_event';
                    event_start = frame;
                end
                
            case 'in_event'
                if signal < lower_thresh
                    % Event candidate found
                    potential_events(end+1, :) = [event_start, frame];
                    state = 'baseline';
                    event_start = 0;
                end
        end
    end
    
    % Handle event extending to end of trace
    if strcmp(state, 'in_event') && event_start > 0
        potential_events(end+1, :) = [event_start, numFrames];
    end
    
    if isempty(potential_events)
        return;  % No events detected
    end
    
    %% === Phase 2: Merge Nearby Events ===
    % Events separated by ≤max_gap_to_merge frames are likely the same biological event
    merged_events = [];
    current_start = potential_events(1, 1);
    current_end = potential_events(1, 2);
    
    for i = 2:size(potential_events, 1)
        next_start = potential_events(i, 1);
        next_end = potential_events(i, 2);
        
        % Check if events should be merged
        gap = next_start - current_end;
        if gap <= max_gap_to_merge
            % Merge events
            current_end = next_end;
        else
            % Save current event and start new one
            merged_events(end+1, :) = [current_start, current_end];
            current_start = next_start;
            current_end = next_end;
        end
    end
    
    % Don't forget the last event
    merged_events(end+1, :) = [current_start, current_end];
    
    %% === Phase 3: Apply Minimum Duration Filter ===
    for i = 1:size(merged_events, 1)
        start_frame = merged_events(i, 1);
        end_frame = merged_events(i, 2);
        duration = end_frame - start_frame + 1;
        
        % Only accept events that meet minimum duration
        if duration >= min_duration
            event_mask(start_frame:end_frame) = true;
        end
    end
end

function extended_mask = extend_events_for_kinetics(event_mask, rise_frames, decay_frames)
    % Extend events for sensor kinetics (now just fine-tuning)
    
    numFrames = length(event_mask);
    extended_mask = event_mask;
    
    if ~any(event_mask)
        return;
    end
    
    % Find event boundaries
    event_starts = find(diff([false; event_mask]) == 1);
    event_ends = find(diff([event_mask; false]) == -1);
    
    % Extend each event
    for i = 1:length(event_starts)
        start_frame = event_starts(i);
        end_frame = event_ends(i);
        
        % Conservative extension
        extended_start = max(1, start_frame - rise_frames);
        extended_end = min(numFrames, end_frame + decay_frames);
        
        extended_mask(extended_start:extended_end) = true;
    end
end

function stats = calculate_enhanced_event_statistics(event_mask, dfof_data, upper_threshold, ...
    lower_threshold, noise_metrics, config)
    % Calculate comprehensive event statistics
    
    [numFrames, numROIs] = size(event_mask);
    
    % Basic counts
    event_frames_per_roi = sum(event_mask, 1);
    rois_with_events = sum(event_frames_per_roi > 0);
    total_event_frames = sum(event_mask, 'all');
    
    % Enhanced event characterization
    event_durations = [];
    event_amplitudes = [];
    events_per_roi = zeros(1, numROIs);
    
    for roi = 1:numROIs
        roi_events = event_mask(:, roi);
        roi_dfof = dfof_data(:, roi);
        
        if any(roi_events)
            % Find discrete events
            event_starts = find(diff([false; roi_events]) == 1);
            event_ends = find(diff([roi_events; false]) == -1);
            events_per_roi(roi) = length(event_starts);
            
            % Analyze each event
            for e = 1:length(event_starts)
                start_frame = event_starts(e);
                end_frame = event_ends(e);
                
                % Duration
                duration = end_frame - start_frame + 1;
                event_durations(end+1) = duration;
                
                % Peak amplitude
                event_segment = roi_dfof(start_frame:end_frame);
                max_amplitude = max(event_segment, [], 'omitnan');
                if ~isnan(max_amplitude)
                    event_amplitudes(end+1) = max_amplitude;
                end
            end
        end
    end
    
    % Compile statistics
    stats = struct();
    
    % Basic counts
    stats.total_events = sum(events_per_roi);
    stats.total_event_frames = total_event_frames;
    stats.rois_with_events = rois_with_events;
    stats.events_per_roi = events_per_roi;
    stats.event_frames_per_roi = event_frames_per_roi;
    
    % Enhanced characteristics
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
    
    % Backward compatibility
    stats.event_summary = struct();
    stats.event_summary.rois_with_events = rois_with_events;
    stats.event_summary.total_events = stats.total_events;
    stats.event_summary.fraction_active_rois = stats.fraction_rois_with_events;
    if rois_with_events > 0
        stats.event_summary.mean_events_per_active_roi = mean(events_per_roi(events_per_roi > 0));
    else
        stats.event_summary.mean_events_per_active_roi = 0;
    end
end