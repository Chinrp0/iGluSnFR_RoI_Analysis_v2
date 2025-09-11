function [event_mask, event_stats] = fixed_pure_schmitt_trigger_detector(dfof_data, config, baseline_stats)
    % FIXED_PURE_SCHMITT_TRIGGER_DETECTOR - Corrected Schmitt trigger implementation
    % Based on the research paper description:
    % "we identified putative fluorescence signals when they exceeded an upper threshold of 3.5σ; 
    % the signals terminated when they decayed below a lower threshold of 1.5σ"
    %
    % FIXES:
    % 1. Proper threshold calculation using baseline periods
    % 2. Correct state machine logic 
    % 3. Events end exactly when signal drops below lower threshold
    % 4. No artificial extensions - pure biological events
    
    if nargin < 2
        config = tracenorm_config();
    end
    
    if nargin < 3
        baseline_stats = [];
    end
    
    [numFrames, numROIs] = size(dfof_data);
    
    % Get parameters from config
    upper_threshold_sigma = config.event_detection.upper_threshold_sigma;  % Default: 3.0
    lower_threshold_sigma = config.event_detection.lower_threshold_sigma;  % Default: 1.5
    
    % Minimum duration filter (still filter obvious noise)
    min_event_duration = 3;  % Minimum 3 frames (75 ms at 40 Hz)
    max_gap_to_merge = 2;    % Merge events separated by ≤2 frames
    
    if config.verbose
        fprintf('FIXED Pure Schmitt trigger detection:\n');
        fprintf('  Thresholds: %.1fσ upper / %.1fσ lower\n', upper_threshold_sigma, lower_threshold_sigma);
        fprintf('  Min duration: %d frames (%.0f ms)\n', min_event_duration, min_event_duration * 1000 / config.frame_rate);
        fprintf('  No kinetic extensions applied\n');
    end
    
    %% === Calculate Thresholds Using Research Paper Method ===
    [upper_threshold, lower_threshold, noise_metrics] = calculate_research_paper_thresholds(...
        dfof_data, upper_threshold_sigma, lower_threshold_sigma, config);
    
    %% === Apply Pure Schmitt Trigger to Each ROI ===
    event_mask = false(numFrames, numROIs);
    
    for roi = 1:numROIs
        roi_trace = dfof_data(:, roi);
        roi_upper = upper_threshold(roi);
        roi_lower = lower_threshold(roi);
        
        % Skip invalid ROIs
        if all(isnan(roi_trace)) || roi_upper <= 0 || isnan(roi_upper)
            continue;
        end
        
        % Apply corrected Schmitt trigger
        roi_events = corrected_schmitt_trigger(roi_trace, roi_upper, roi_lower, ...
            min_event_duration, max_gap_to_merge);
        
        event_mask(:, roi) = roi_events;
    end
    
    %% === Calculate Event Statistics ===
    event_stats = calculate_enhanced_event_statistics(event_mask, dfof_data, upper_threshold, ...
        lower_threshold, noise_metrics, config);
    
    % Add metadata
    event_stats.event_mask = event_mask;
    event_stats.noise_metrics = noise_metrics;
    event_stats.method = 'fixed_pure_schmitt';
    event_stats.kinetic_extension_applied = false;
    event_stats.min_duration_filter = min_event_duration;
    
    if config.verbose
        fprintf('  RESULTS: %d total events across %d ROIs\n', ...
            event_stats.total_events, event_stats.rois_with_events);
        fprintf('  Mean duration: %.1f frames (%.0f ms)\n', ...
            event_stats.mean_event_duration, event_stats.mean_event_duration * 1000 / config.frame_rate);
    end
end

function [upper_threshold, lower_threshold, noise_metrics] = calculate_research_paper_thresholds(...
    dfof_data, upper_sigma, lower_sigma, config)
    % Calculate thresholds using the research paper approach
    % Key insight: Use baseline periods to estimate noise, excluding events
    
    [~, numROIs] = size(dfof_data);
    
    upper_threshold = zeros(1, numROIs);
    lower_threshold = zeros(1, numROIs);
    noise_std = zeros(1, numROIs);
    
    if config.verbose && numROIs > 1000
        fprintf('    Calculating thresholds for %d ROIs...\n', numROIs);
    end
    
    for roi = 1:numROIs
        roi_trace = dfof_data(:, roi);
        
        % Handle completely invalid ROIs
        if all(isnan(roi_trace))
            noise_std(roi) = 0.01;  % Default noise level
            upper_threshold(roi) = upper_sigma * 0.01;
            lower_threshold(roi) = lower_sigma * 0.01;
            continue;
        end
        
        % Remove NaN values
        valid_data = roi_trace(~isnan(roi_trace));
        
        if length(valid_data) < 20
            % Too little data - use simple MAD
            noise_estimate = mad(valid_data, 1) * 1.4826;
        else
            % Research paper approach: Use lower percentiles to exclude events
            % This estimates noise from baseline periods, not from event periods
            baseline_percentile = 70;  % Use bottom 70% of data for noise estimation
            threshold_70 = prctile(valid_data, baseline_percentile);
            
            % Get baseline periods (likely non-event data)
            baseline_data = valid_data(valid_data <= threshold_70);
            
            if length(baseline_data) > 10
                % Use MAD (median absolute deviation) for robust noise estimation
                noise_estimate = mad(baseline_data, 1) * 1.4826;  % MAD → std conversion
            else
                % Fallback to full data MAD
                noise_estimate = mad(valid_data, 1) * 1.4826;
            end
        end
        
        % Store results
        noise_std(roi) = noise_estimate;
        upper_threshold(roi) = upper_sigma * noise_estimate;
        lower_threshold(roi) = lower_sigma * noise_estimate;
    end
    
    % Compile noise metrics
    noise_metrics = struct();
    noise_metrics.noise_std = noise_std;
    noise_metrics.upper_thresholds = upper_threshold;
    noise_metrics.lower_thresholds = lower_threshold;
    noise_metrics.threshold_separation = upper_threshold - lower_threshold;
    noise_metrics.mean_noise_std = mean(noise_std, 'omitnan');
    noise_metrics.median_noise_std = median(noise_std, 'omitnan');
end

function event_mask = corrected_schmitt_trigger(trace, upper_thresh, lower_thresh, ...
    min_duration, max_gap_to_merge)
    % CORRECTED Schmitt trigger implementation
    % Events start when crossing ABOVE upper threshold
    % Events end when dropping BELOW lower threshold
    
    numFrames = length(trace);
    event_mask = false(numFrames, 1);
    
    %% === Phase 1: Pure Schmitt Trigger State Machine ===
    state = 'baseline';
    event_start = 0;
    raw_events = [];  % Store [start, end] pairs
    
    for frame = 1:numFrames
        signal = trace(frame);
        
        % Skip NaN values but maintain state
        if isnan(signal)
            continue;
        end
        
        switch state
            case 'baseline'
                if signal > upper_thresh
                    % Event starts when crossing above upper threshold
                    state = 'in_event';
                    event_start = frame;
                end
                
            case 'in_event'
                if signal < lower_thresh
                    % Event ends when dropping below lower threshold
                    % Include the frame where it was still above lower threshold
                    event_end = frame - 1;  % Last frame above lower threshold
                    
                    % Store this event
                    if event_end >= event_start  % Valid event
                        raw_events(end+1, :) = [event_start, event_end];
                    end
                    
                    % Return to baseline state
                    state = 'baseline';
                    event_start = 0;
                end
        end
    end
    
    % Handle event that extends to end of trace
    if strcmp(state, 'in_event') && event_start > 0
        raw_events(end+1, :) = [event_start, numFrames];
    end
    
    % Early return if no events detected
    if isempty(raw_events)
        return;
    end
    
    %% === Phase 2: Merge Nearby Events (Biological Realism) ===
    % Events separated by ≤ max_gap_to_merge frames are likely the same biological event
    merged_events = [];
    current_start = raw_events(1, 1);
    current_end = raw_events(1, 2);
    
    for i = 2:size(raw_events, 1)
        next_start = raw_events(i, 1);
        next_end = raw_events(i, 2);
        
        % Calculate gap between current event end and next event start
        gap = next_start - current_end - 1;
        
        if gap <= max_gap_to_merge
            % Merge events by extending current event
            current_end = next_end;
        else
            % Save current event and start tracking new one
            merged_events(end+1, :) = [current_start, current_end];
            current_start = next_start;
            current_end = next_end;
        end
    end
    
    % Don't forget the last event
    merged_events(end+1, :) = [current_start, current_end];
    
    %% === Phase 3: Apply Minimum Duration Filter ===
    % Filter out very short events (likely noise spikes)
    for i = 1:size(merged_events, 1)
        start_frame = merged_events(i, 1);
        end_frame = merged_events(i, 2);
        duration = end_frame - start_frame + 1;
        
        % Only accept events that meet minimum duration requirement
        if duration >= min_duration
            event_mask(start_frame:end_frame) = true;
        end
    end
end

function stats = calculate_enhanced_event_statistics(event_mask, dfof_data, upper_threshold, ...
    lower_threshold, noise_metrics, config)
    % Calculate comprehensive event statistics for the corrected implementation
    
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
            % Find discrete events in this ROI
            event_starts = find(diff([false; roi_events]) == 1);
            event_ends = find(diff([roi_events; false]) == -1);
            events_per_roi(roi) = length(event_starts);
            
            % Analyze each discrete event
            for e = 1:length(event_starts)
                start_frame = event_starts(e);
                end_frame = event_ends(e);
                
                % Event duration
                duration = end_frame - start_frame + 1;
                event_durations(end+1) = duration;
                
                % Peak amplitude during event
                event_segment = roi_dfof(start_frame:end_frame);
                if any(~isnan(event_segment))
                    max_amplitude = max(event_segment, [], 'omitnan');
                    event_amplitudes(end+1) = max_amplitude;
                end
            end
        end
    end
    
    % Compile comprehensive statistics
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
        stats.std_event_duration = std(event_durations);
        stats.event_duration_range = [min(event_durations), max(event_durations)];
        stats.event_durations = event_durations;  % For detailed analysis
    else
        stats.mean_event_duration = 0;
        stats.median_event_duration = 0;
        stats.std_event_duration = 0;
        stats.event_duration_range = [0, 0];
        stats.event_durations = [];
    end
    
    if ~isempty(event_amplitudes)
        stats.mean_event_amplitude = mean(event_amplitudes);
        stats.median_event_amplitude = median(event_amplitudes);
        stats.std_event_amplitude = std(event_amplitudes);
        stats.event_amplitude_range = [min(event_amplitudes), max(event_amplitudes)];
        stats.event_amplitudes = event_amplitudes;  % For detailed analysis
    else
        stats.mean_event_amplitude = 0;
        stats.median_event_amplitude = 0;
        stats.std_event_amplitude = 0;
        stats.event_amplitude_range = [0, 0];
        stats.event_amplitudes = [];
    end
    
    % Threshold information
    stats.upper_thresholds = upper_threshold;
    stats.lower_thresholds = lower_threshold;
    stats.mean_upper_threshold = mean(upper_threshold, 'omitnan');
    stats.mean_lower_threshold = mean(lower_threshold, 'omitnan');
    
    % Detection efficiency metrics
    stats.fraction_rois_with_events = rois_with_events / numROIs;
    stats.fraction_frames_in_events = total_event_frames / numel(event_mask);
    
    % Activity summary
    if rois_with_events > 0
        stats.mean_events_per_active_roi = mean(events_per_roi(events_per_roi > 0));
        stats.median_events_per_active_roi = median(events_per_roi(events_per_roi > 0));
    else
        stats.mean_events_per_active_roi = 0;
        stats.median_events_per_active_roi = 0;
    end
    
    % Backward compatibility with existing pipeline
    stats.event_summary = struct();
    stats.event_summary.rois_with_events = rois_with_events;
    stats.event_summary.total_events = stats.total_events;
    stats.event_summary.fraction_active_rois = stats.fraction_rois_with_events;
    stats.event_summary.mean_events_per_active_roi = stats.mean_events_per_active_roi;
end
