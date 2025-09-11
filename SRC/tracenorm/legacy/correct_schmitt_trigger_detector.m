function [event_mask, event_stats] = correct_schmitt_trigger_detector(dfof_data, config, baseline_stats)
    % CORRECT_SCHMITT_TRIGGER_DETECTOR - Proper biological event detection
    % 
    % CORRECT LOGIC:
    % 1. Event starts when signal crosses above upper threshold
    % 2. Signal must stay above lower threshold for ≥3 frames (validation period)
    % 3. If validation fails → reject as noise
    % 4. If validation passes → continue until signal drops below lower threshold
    % 5. Event ends when signal finally drops below lower threshold
    %
    % This ensures only sustained biological events are detected, not noise spikes.
    
    if nargin < 2
        config = tracenorm_config();
    end
    
    if nargin < 3
        baseline_stats = [];
    end
    
    [numFrames, numROIs] = size(dfof_data);
    
    % Get parameters
    upper_threshold_sigma = config.event_detection.upper_threshold_sigma;  % 3.0
    lower_threshold_sigma = config.event_detection.lower_threshold_sigma;  % 1.5
    min_sustained_frames = 3;  % Must stay above lower threshold for ≥3 frames
    max_gap_to_merge = 2;      % Merge events separated by ≤2 frames
    
    if config.verbose
        fprintf('CORRECT Schmitt trigger with validation period:\n');
        fprintf('  Thresholds: %.1fσ upper / %.1fσ lower\n', upper_threshold_sigma, lower_threshold_sigma);
        fprintf('  Validation: ≥%d frames above lower threshold required\n', min_sustained_frames);
        fprintf('  Merge gap: ≤%d frames\n', max_gap_to_merge);
    end
    
    %% === Calculate Thresholds ===
    [upper_threshold, lower_threshold, noise_metrics] = calculate_research_thresholds(...
        dfof_data, upper_threshold_sigma, lower_threshold_sigma, config);
    
    %% === Apply Correct Schmitt Trigger to Each ROI ===
    event_mask = false(numFrames, numROIs);
    debug_roi = 613;  % For detailed debugging
    
    for roi = 1:numROIs
        roi_trace = dfof_data(:, roi);
        roi_upper = upper_threshold(roi);
        roi_lower = lower_threshold(roi);
        
        % Skip invalid ROIs
        if all(isnan(roi_trace)) || roi_upper <= 0 || isnan(roi_upper)
            continue;
        end
        
        % Apply correct Schmitt trigger with validation
        roi_events = correct_schmitt_with_validation(roi_trace, roi_upper, roi_lower, ...
            min_sustained_frames, max_gap_to_merge, roi == debug_roi && config.verbose);
        
        event_mask(:, roi) = roi_events;
    end
    
    %% === Calculate Event Statistics ===
    event_stats = calculate_event_statistics(event_mask, dfof_data, upper_threshold, ...
        lower_threshold, noise_metrics, config);
    
    % Add metadata
    event_stats.event_mask = event_mask;
    event_stats.noise_metrics = noise_metrics;
    event_stats.method = 'correct_schmitt_with_validation';
    event_stats.min_sustained_frames = min_sustained_frames;
    
    if config.verbose
        fprintf('  RESULTS: %d total events across %d ROIs\n', ...
            event_stats.total_events, event_stats.rois_with_events);
        fprintf('  Mean duration: %.1f frames (%.0f ms)\n', ...
            event_stats.mean_event_duration, event_stats.mean_event_duration * 1000 / config.frame_rate);
    end
end

function [upper_threshold, lower_threshold, noise_metrics] = calculate_research_thresholds(...
    dfof_data, upper_sigma, lower_sigma, config)
    % Calculate thresholds using research paper method (baseline-only noise estimation)
    
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
        
        % Research paper approach: exclude likely events from noise calculation
        valid_data = roi_trace(~isnan(roi_trace));
        
        if length(valid_data) > 20
            % Use bottom 70% of data (excludes events from noise estimation)
            baseline_70th = prctile(valid_data, 70);
            baseline_data = valid_data(valid_data <= baseline_70th);
            
            if length(baseline_data) > 10
                noise_estimate = mad(baseline_data, 1) * 1.4826;
            else
                noise_estimate = mad(valid_data, 1) * 1.4826;
            end
        else
            noise_estimate = mad(valid_data, 1) * 1.4826;
        end
        
        noise_std(roi) = noise_estimate;
        upper_threshold(roi) = upper_sigma * noise_estimate;
        lower_threshold(roi) = lower_sigma * noise_estimate;
    end
    
    % Compile noise metrics
    noise_metrics = struct();
    noise_metrics.noise_std = noise_std;
    noise_metrics.upper_thresholds = upper_threshold;
    noise_metrics.lower_thresholds = lower_threshold;
    noise_metrics.mean_noise_std = mean(noise_std, 'omitnan');
end

function event_mask = correct_schmitt_with_validation(trace, upper_thresh, lower_thresh, ...
    min_sustained_frames, max_gap_to_merge, debug_mode)
    % CORRECT Schmitt trigger with validation period
    
    numFrames = length(trace);
    event_mask = false(numFrames, 1);
    
    if debug_mode
        fprintf('\n--- DEBUG ROI 613 Correct Schmitt Trigger ---\n');
        fprintf('Thresholds: upper=%.6f, lower=%.6f\n', upper_thresh, lower_thresh);
    end
    
    %% === Phase 1: Schmitt Trigger with Validation ===
    state = 'baseline';
    event_start = 0;
    validated_events = [];  % Store [start, end] pairs for validated events
    
    for frame = 1:numFrames
        signal = trace(frame);
        
        if isnan(signal)
            continue;
        end
        
        if strcmp(state, 'baseline')
            if signal > upper_thresh
                % Potential event start
                state = 'validating';
                event_start = frame;
                
                if debug_mode
                    fprintf('Potential event START at frame %d (%.3fs): signal=%.6f > %.6f\n', ...
                        frame, frame/40, signal, upper_thresh);
                end
            end
            
        elseif strcmp(state, 'validating')
            % Check if we've sustained above lower threshold long enough
            frames_since_start = frame - event_start + 1;
            
            if signal < lower_thresh
                % Failed validation - signal dropped below lower threshold too quickly
                if debug_mode
                    fprintf('  VALIDATION FAILED at frame %d (%.3fs): signal=%.6f < %.6f after %d frames\n', ...
                        frame, frame/40, signal, lower_thresh, frames_since_start);
                end
                state = 'baseline';
                event_start = 0;
                
            elseif frames_since_start >= min_sustained_frames
                % Validation passed - this is a real event!
                state = 'validated_event';
                
                if debug_mode
                    fprintf('  VALIDATION PASSED at frame %d: sustained ≥%d frames, continuing event...\n', ...
                        frame, min_sustained_frames);
                end
            end
            % If still validating and above lower threshold, continue
            
        elseif strcmp(state, 'validated_event')
            if signal < lower_thresh
                % Event ends - we've validated it and now it's dropping below lower threshold
                event_end = frame - 1;  % Last frame above lower threshold
                validated_events(end+1, :) = [event_start, event_end];
                
                duration = event_end - event_start + 1;
                if debug_mode
                    fprintf('  EVENT END at frame %d (%.3fs): signal=%.6f < %.6f, duration=%d frames\n', ...
                        event_end, event_end/40, signal, lower_thresh, duration);
                end
                
                state = 'baseline';
                event_start = 0;
            end
            % If still above lower threshold, continue the validated event
        end
    end
    
    % Handle event extending to end of trace
    if strcmp(state, 'validated_event') && event_start > 0
        validated_events(end+1, :) = [event_start, numFrames];
        if debug_mode
            fprintf('  EVENT extends to end: duration=%d frames\n', numFrames - event_start + 1);
        end
    elseif strcmp(state, 'validating') && debug_mode
        fprintf('  Potential event at end rejected: insufficient validation\n');
    end
    
    if debug_mode
        fprintf('Found %d validated events before merging\n', size(validated_events, 1));
    end
    
    if isempty(validated_events)
        return;  % No validated events
    end
    
    %% === Phase 2: Merge Nearby Validated Events ===
    merged_events = [];
    current_start = validated_events(1, 1);
    current_end = validated_events(1, 2);
    
    for i = 2:size(validated_events, 1)
        next_start = validated_events(i, 1);
        next_end = validated_events(i, 2);
        
        gap = next_start - current_end - 1;
        if gap <= max_gap_to_merge
            % Merge events
            current_end = next_end;
            if debug_mode
                fprintf('  Merging events: gap=%d frames ≤ %d\n', gap, max_gap_to_merge);
            end
        else
            % Save current event and start new one
            merged_events(end+1, :) = [current_start, current_end];
            current_start = next_start;
            current_end = next_end;
        end
    end
    merged_events(end+1, :) = [current_start, current_end];
    
    if debug_mode
        fprintf('After merging: %d final events\n', size(merged_events, 1));
    end
    
    %% === Phase 3: Create Event Mask ===
    for i = 1:size(merged_events, 1)
        start_frame = merged_events(i, 1);
        end_frame = merged_events(i, 2);
        event_mask(start_frame:end_frame) = true;
        
        if debug_mode
            duration = end_frame - start_frame + 1;
            fprintf('  Final Event %d: [%d-%d] duration=%d frames (%.0f ms)\n', ...
                i, start_frame, end_frame, duration, duration * 1000 / 40);
        end
    end
end

function stats = calculate_event_statistics(event_mask, dfof_data, upper_threshold, ...
    lower_threshold, noise_metrics, config)
    % Calculate comprehensive event statistics
    
    [numFrames, numROIs] = size(event_mask);
    
    % Basic counts
    event_frames_per_roi = sum(event_mask, 1);
    rois_with_events = sum(event_frames_per_roi > 0);
    total_event_frames = sum(event_mask, 'all');
    
    % Detailed event analysis
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
                if any(~isnan(event_segment))
                    max_amplitude = max(event_segment, [], 'omitnan');
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
