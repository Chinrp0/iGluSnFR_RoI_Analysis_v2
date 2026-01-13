function [event_mask, event_stats] = corrected_schmitt_detector(dfof_data, config, baseline_stats)
    % CORRECTED_SCHMITT_DETECTOR - Proper biological event detection with corrected noise estimation
    % 
    % KEY CORRECTION: Include outliers in noise calculation (they ARE noise, not events)
    % Only exclude sustained biological events (≥7 frames) from noise calculation
    %
    % Algorithm:
    % 1. Use baseline detector's outlier mask - these ARE noise, include in calculation
    % 2. Identify sustained events (≥7 frames above 85th percentile) - exclude from noise
    % 3. Calculate noise from: baseline periods + outliers - sustained events  
    % 4. Apply Schmitt trigger with 3-frame validation period
    % 5. Events start at upper threshold crossing, end when dropping below lower threshold
    
    if nargin < 2
        config = tracenorm_config();
    end
    
    if nargin < 3
        baseline_stats = [];
        warning('No baseline_stats provided - creating fallback outlier mask');
    end
    
    [numFrames, numROIs] = size(dfof_data);
    
    % Get corrected parameters
    upper_threshold_sigma = config.corrected_schmitt.upper_threshold_sigma;  % 3.5
    lower_threshold_sigma = config.corrected_schmitt.lower_threshold_sigma;  % 1.5
    noise_exclusion_window = config.corrected_schmitt.noise_exclusion_window; % 7 frames
    min_event_duration = config.corrected_schmitt.min_event_duration;       % 3 frames
    sustained_percentile = config.corrected_schmitt.sustained_percentile;   % 85th percentile
    
    if config.verbose
        fprintf('CORRECTED Schmitt trigger with proper noise estimation:\n');
        fprintf('  Thresholds: %.1fσ upper / %.1fσ lower\n', upper_threshold_sigma, lower_threshold_sigma);
        fprintf('  Noise calculation: INCLUDES outliers, EXCLUDES sustained events (≥%d frames)\n', noise_exclusion_window);
        fprintf('  Validation: ≥%d frames above lower threshold required\n', min_event_duration);
    end
    
    %% === Extract Outlier Mask from Baseline Detector ===
    if isfield(baseline_stats, 'outlier_mask') && ~isempty(baseline_stats.outlier_mask)
        outlier_mask = baseline_stats.outlier_mask;  % [frames x ROIs]
        if config.verbose
            fprintf('  Using baseline detector outlier mask: %.1f%% outliers\n', ...
                100 * sum(outlier_mask, 'all') / numel(outlier_mask));
        end
    else
        % Fallback: create simple outlier mask using MAD
        outlier_mask = create_fallback_outlier_mask(dfof_data, config);
        warning('Created fallback outlier mask - recommend providing baseline_stats');
    end
    
    %% === Calculate CORRECTED Thresholds ===
    [upper_threshold, lower_threshold, noise_metrics] = calculate_corrected_thresholds(...
        dfof_data, outlier_mask, upper_threshold_sigma, lower_threshold_sigma, ...
        sustained_percentile, noise_exclusion_window, config);
    
    %% === Apply Corrected Schmitt Trigger to Each ROI ===
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
        
        % Apply corrected Schmitt trigger with validation
        roi_events = corrected_schmitt_with_validation(roi_trace, roi_upper, roi_lower, ...
            min_event_duration, roi == debug_roi && config.verbose);
        
        event_mask(:, roi) = roi_events;
    end
    
    %% === Calculate Event Statistics ===
    event_stats = calculate_corrected_event_statistics(event_mask, dfof_data, upper_threshold, ...
        lower_threshold, noise_metrics, config);
    
    % Add metadata
    event_stats.event_mask = event_mask;
    event_stats.noise_metrics = noise_metrics;
    event_stats.method = 'corrected_schmitt_detector';
    event_stats.includes_outliers_in_noise = true;
    event_stats.excludes_sustained_events = true;
    event_stats.min_sustained_frames = noise_exclusion_window;
    event_stats.min_validation_frames = min_event_duration;
    
    if config.verbose
        fprintf('  RESULTS: %d total events across %d ROIs\n', ...
            event_stats.total_events, event_stats.rois_with_events);
        fprintf('  Mean duration: %.1f frames (%.0f ms)\n', ...
            event_stats.mean_event_duration, event_stats.mean_event_duration * 1000 / config.frame_rate);
        
        % NEW: Display frequency and peak amplitude information
        if event_stats.rois_with_events > 0
            fprintf('  Event frequency: %.3f Hz avg (range: %.3f-%.3f Hz)\n', ...
                event_stats.mean_frequency_across_active_rois, ...
                event_stats.frequency_range(1), event_stats.frequency_range(2));
            fprintf('  Peak amplitude: %.4f dF/F avg (range: %.4f-%.4f dF/F)\n', ...
                event_stats.mean_peak_amplitude_across_active_rois, ...
                event_stats.peak_amplitude_range_across_rois(1), event_stats.peak_amplitude_range_across_rois(2));
        else
            fprintf('  Event frequency: 0 Hz (no active ROIs)\n');
            fprintf('  Peak amplitude: 0 dF/F (no events detected)\n');
        end
        
        fprintf('  Recording time: %.1f s (%.1f min)\n', ...
            event_stats.recording_time_seconds, event_stats.recording_time_seconds / 60);
        fprintf('  Noise method: CORRECTED (outliers included, sustained excluded)\n');
    end
end

function [upper_threshold, lower_threshold, noise_metrics] = calculate_corrected_thresholds(...
    dfof_data, outlier_mask, upper_sigma, lower_sigma, sustained_percentile, ...
    exclusion_window, config)
    % CORRECTED threshold calculation: Include outliers, exclude sustained events
    
    [numFrames, numROIs] = size(dfof_data);
    
    upper_threshold = zeros(1, numROIs);
    lower_threshold = zeros(1, numROIs);
    noise_std = zeros(1, numROIs);
    noise_method = cell(1, numROIs);
    
    if config.verbose && numROIs > 1000
        fprintf('    Calculating CORRECTED thresholds for %d ROIs...\n', numROIs);
    end
    
    for roi = 1:numROIs
        roi_trace = dfof_data(:, roi);
        roi_outliers = outlier_mask(:, roi);
        
        % Handle completely invalid ROIs
        if all(isnan(roi_trace))
            noise_std(roi) = 0.01;
            upper_threshold(roi) = upper_sigma * 0.01;
            lower_threshold(roi) = lower_sigma * 0.01;
            noise_method{roi} = 'fallback_default';
            continue;
        end
        
        % STEP 1: Get valid frames (non-NaN)
        valid_frames = ~isnan(roi_trace);
        
        % STEP 2: Get outlier frames (these ARE noise - include them)
        outlier_frames = valid_frames & roi_outliers;
        
        % STEP 3: Get baseline frames (valid but not outliers)
        baseline_frames = valid_frames & ~roi_outliers;
        
        if sum(baseline_frames) < 20
            % Too few baseline frames - use all valid data
            noise_data = roi_trace(valid_frames);
            noise_method{roi} = 'all_valid_data';
        else
            % STEP 4: Identify sustained events to EXCLUDE from noise calculation
            baseline_data = roi_trace(baseline_frames);
            sustained_threshold = prctile(baseline_data, sustained_percentile);
            
            % Find sustained periods (≥exclusion_window frames above threshold)
            sustained_mask = identify_sustained_periods(roi_trace, sustained_threshold, exclusion_window);
            
            % STEP 5: CORRECTED noise calculation
            % Include: baseline frames + outlier frames
            % Exclude: sustained event frames
            noise_frames = (baseline_frames | outlier_frames) & ~sustained_mask;
            
            if sum(noise_frames) < 10
                % Fallback if too few noise frames
                noise_data = roi_trace(baseline_frames | outlier_frames);
                noise_method{roi} = 'baseline_plus_outliers';
            else
                noise_data = roi_trace(noise_frames);
                noise_method{roi} = 'corrected_with_exclusion';
            end
        end
        
        % Calculate robust noise estimate using MAD
        if length(noise_data) > 5
            noise_estimate = mad(noise_data, 1) * 1.4826;  % MAD → std conversion
        else
            noise_estimate = 0.01;  % Minimal fallback
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
    noise_metrics.noise_method = noise_method;
    noise_metrics.mean_noise_std = mean(noise_std, 'omitnan');
    noise_metrics.median_noise_std = median(noise_std, 'omitnan');
    noise_metrics.outliers_included = true;
    noise_metrics.sustained_events_excluded = true;
    
    % Calculate signal quality metrics
    signal_peaks = max(dfof_data, [], 1, 'omitnan');
    signal_quality = signal_peaks ./ noise_std;
    noise_metrics.signal_quality = signal_quality;
    noise_metrics.mean_signal_quality = mean(signal_quality, 'omitnan');
end

function sustained_mask = identify_sustained_periods(trace, threshold, min_frames)
    % Identify sustained periods above threshold (≥min_frames duration)
    
    sustained_mask = false(size(trace));
    
    % Find periods above threshold
    above_threshold = trace > threshold & ~isnan(trace);
    
    if ~any(above_threshold)
        return;
    end
    
    % Find start and end of periods
    period_starts = find(diff([false; above_threshold]) == 1);
    period_ends = find(diff([above_threshold; false]) == -1);
    
    % Mark sustained periods (≥min_frames duration)
    for i = 1:length(period_starts)
        start_frame = period_starts(i);
        end_frame = period_ends(i);
        duration = end_frame - start_frame + 1;
        
        if duration >= min_frames
            sustained_mask(start_frame:end_frame) = true;
        end
    end
end

function outlier_mask = create_fallback_outlier_mask(dfof_data, config)
    % Create simple outlier mask as fallback when baseline_stats unavailable
    
    [numFrames, numROIs] = size(dfof_data);
    outlier_mask = false(numFrames, numROIs);
    
    for roi = 1:numROIs
        roi_trace = dfof_data(:, roi);
        valid_data = roi_trace(~isnan(roi_trace));
        
        if length(valid_data) > 10
            median_val = median(valid_data);
            mad_val = mad(valid_data, 1);
            threshold = median_val + 2.0 * mad_val * 1.4826;  % 2σ threshold
            
            outlier_mask(:, roi) = roi_trace > threshold & ~isnan(roi_trace);
        end
    end
end

function event_mask = corrected_schmitt_with_validation(trace, upper_thresh, lower_thresh, ...
    min_sustained_frames, debug_mode)
    % Corrected Schmitt trigger with proper validation period
    
    numFrames = length(trace);
    event_mask = false(numFrames, 1);
    
    if debug_mode
        fprintf('\n--- DEBUG ROI 613 CORRECTED Schmitt Trigger ---\n');
        fprintf('Thresholds: upper=%.6f, lower=%.6f\n', upper_thresh, lower_thresh);
    end
    
    %% === Phase 1: Schmitt Trigger with Validation ===
    state = 'baseline';
    event_start = 0;
    validated_events = [];
    
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
                    fprintf('Potential event START at frame %d: signal=%.6f > %.6f\n', ...
                        frame, signal, upper_thresh);
                end
            end
            
        elseif strcmp(state, 'validating')
            % Check validation period
            frames_since_start = frame - event_start + 1;
            
            if signal < lower_thresh
                % Failed validation - signal dropped too quickly
                if debug_mode
                    fprintf('  VALIDATION FAILED at frame %d: signal=%.6f < %.6f after %d frames\n', ...
                        frame, signal, lower_thresh, frames_since_start);
                end
                state = 'baseline';
                event_start = 0;
                
            elseif frames_since_start >= min_sustained_frames
                % Validation passed - this is a real event!
                state = 'validated_event';
                
                if debug_mode
                    fprintf('  VALIDATION PASSED at frame %d: sustained ≥%d frames\n', ...
                        frame, min_sustained_frames);
                end
            end
            
        elseif strcmp(state, 'validated_event')
            if signal < lower_thresh
                % Event ends
                event_end = frame - 1;
                validated_events(end+1, :) = [event_start, event_end];
                
                duration = event_end - event_start + 1;
                if debug_mode
                    fprintf('  EVENT END at frame %d: duration=%d frames\n', event_end, duration);
                end
                
                state = 'baseline';
                event_start = 0;
            end
        end
    end
    
    % Handle event extending to end
    if strcmp(state, 'validated_event') && event_start > 0
        validated_events(end+1, :) = [event_start, numFrames];
    end
    
    if debug_mode
        fprintf('Total validated events: %d\n', size(validated_events, 1));
    end
    
    %% === Phase 2: Merge Nearby Events (≤2 frame gap) ===
    if isempty(validated_events)
        return;
    end
    
    merged_events = [];
    current_start = validated_events(1, 1);
    current_end = validated_events(1, 2);
    
    for i = 2:size(validated_events, 1)
        next_start = validated_events(i, 1);
        next_end = validated_events(i, 2);
        
        gap = next_start - current_end - 1;
        if gap <= 2  % Merge nearby events
            current_end = next_end;
        else
            merged_events(end+1, :) = [current_start, current_end];
            current_start = next_start;
            current_end = next_end;
        end
    end
    merged_events(end+1, :) = [current_start, current_end];
    
    %% === Phase 3: Create Final Event Mask ===
    for i = 1:size(merged_events, 1)
        start_frame = merged_events(i, 1);
        end_frame = merged_events(i, 2);
        event_mask(start_frame:end_frame) = true;
        
        if debug_mode
            duration = end_frame - start_frame + 1;
            fprintf('  Final Event %d: [%d-%d] duration=%d frames\n', ...
                i, start_frame, end_frame, duration);
        end
    end
end

function stats = calculate_corrected_event_statistics(event_mask, dfof_data, upper_threshold, ...
    lower_threshold, noise_metrics, config)
    % Calculate comprehensive event statistics with ENHANCED per-ROI peak amplitudes and frequency
    
    [numFrames, numROIs] = size(event_mask);
    
    % Calculate total recording time
    total_time_seconds = numFrames / config.frame_rate;  % For frequency calculation
    
    % Basic counts
    event_frames_per_roi = sum(event_mask, 1);
    rois_with_events = sum(event_frames_per_roi > 0);
    total_event_frames = sum(event_mask, 'all');
    
    % ENHANCED: Per-ROI detailed event analysis
    event_durations = [];  % Flat array for backward compatibility
    event_amplitudes = []; % Flat array for backward compatibility
    events_per_roi = zeros(1, numROIs);
    
    % NEW: Per-ROI data structures
    peak_amplitudes_per_roi = cell(1, numROIs);     % Cell array: {roi} = [peak1, peak2, ...]
    mean_peak_amplitude_per_roi = zeros(1, numROIs); % Mean peak amplitude per ROI
    frequency_hz = zeros(1, numROIs);                % Event frequency per ROI
    event_durations_per_roi = cell(1, numROIs);     % Duration of each event per ROI
    
    for roi = 1:numROIs
        roi_events = event_mask(:, roi);
        roi_dfof = dfof_data(:, roi);
        
        % Initialize per-ROI storage
        roi_peak_amplitudes = [];
        roi_event_durations = [];
        
        if any(roi_events)
            % Find discrete events
            event_starts = find(diff([false; roi_events]) == 1);
            event_ends = find(diff([roi_events; false]) == -1);
            events_per_roi(roi) = length(event_starts);
            
            % Analyze each event for this ROI
            for e = 1:length(event_starts)
                start_frame = event_starts(e);
                end_frame = event_ends(e);
                
                % Duration
                duration = end_frame - start_frame + 1;
                event_durations(end+1) = duration;  % Flat array (backward compatibility)
                roi_event_durations(end+1) = duration;  % Per-ROI array
                
                % Peak amplitude within event
                event_segment = roi_dfof(start_frame:end_frame);
                if any(~isnan(event_segment))
                    max_amplitude = max(event_segment, [], 'omitnan');
                    event_amplitudes(end+1) = max_amplitude;  % Flat array (backward compatibility)
                    roi_peak_amplitudes(end+1) = max_amplitude;  % Per-ROI array
                end
            end
        end
        
        % Store per-ROI data
        peak_amplitudes_per_roi{roi} = roi_peak_amplitudes;
        event_durations_per_roi{roi} = roi_event_durations;
        
        % Calculate mean peak amplitude for this ROI
        if ~isempty(roi_peak_amplitudes)
            mean_peak_amplitude_per_roi(roi) = mean(roi_peak_amplitudes);
        else
            mean_peak_amplitude_per_roi(roi) = 0;  % No events
        end
        
        % Calculate frequency in Hz for this ROI
        frequency_hz(roi) = events_per_roi(roi) / total_time_seconds;
    end
    
    % Compile statistics
    stats = struct();
    
    % Basic counts (unchanged for backward compatibility)
    stats.total_events = sum(events_per_roi);
    stats.total_event_frames = total_event_frames;
    stats.rois_with_events = rois_with_events;
    stats.events_per_roi = events_per_roi;
    stats.event_frames_per_roi = event_frames_per_roi;
    
    % Event characteristics (backward compatibility)
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
    
    % NEW: Enhanced per-ROI event metrics
    stats.peak_amplitudes_per_roi = peak_amplitudes_per_roi;        % Cell array: {roi} = [peak1, peak2, ...]
    stats.mean_peak_amplitude_per_roi = mean_peak_amplitude_per_roi; % [1 x numROIs] mean peak per ROI
    stats.frequency_hz = frequency_hz;                              % [1 x numROIs] event frequency in Hz
    stats.event_durations_per_roi = event_durations_per_roi;        % Cell array: {roi} = [dur1, dur2, ...]
    
    % NEW: Summary metrics for comparison analysis
    stats.recording_time_seconds = total_time_seconds;
    stats.active_roi_frequencies = frequency_hz(frequency_hz > 0);  % Only ROIs with events
    stats.active_roi_mean_peaks = mean_peak_amplitude_per_roi(mean_peak_amplitude_per_roi > 0);  % Only ROIs with events
    
    if ~isempty(stats.active_roi_frequencies)
        stats.mean_frequency_across_active_rois = mean(stats.active_roi_frequencies);
        stats.median_frequency_across_active_rois = median(stats.active_roi_frequencies);
        stats.frequency_range = [min(stats.active_roi_frequencies), max(stats.active_roi_frequencies)];
    else
        stats.mean_frequency_across_active_rois = 0;
        stats.median_frequency_across_active_rois = 0;
        stats.frequency_range = [0, 0];
    end
    
    if ~isempty(stats.active_roi_mean_peaks)
        stats.mean_peak_amplitude_across_active_rois = mean(stats.active_roi_mean_peaks);
        stats.median_peak_amplitude_across_active_rois = median(stats.active_roi_mean_peaks);
        stats.peak_amplitude_range_across_rois = [min(stats.active_roi_mean_peaks), max(stats.active_roi_mean_peaks)];
    else
        stats.mean_peak_amplitude_across_active_rois = 0;
        stats.median_peak_amplitude_across_active_rois = 0;
        stats.peak_amplitude_range_across_rois = [0, 0];
    end
    
    % Threshold information (unchanged)
    stats.upper_thresholds = upper_threshold;
    stats.lower_thresholds = lower_threshold;
    stats.mean_upper_threshold = mean(upper_threshold, 'omitnan');
    stats.mean_lower_threshold = mean(lower_threshold, 'omitnan');
    
    % Detection efficiency (unchanged)
    stats.fraction_rois_with_events = rois_with_events / numROIs;
    stats.fraction_frames_in_events = total_event_frames / numel(event_mask);
    
    % Backward compatibility (unchanged)
    stats.event_summary = struct();
    stats.event_summary.rois_with_events = rois_with_events;
    stats.event_summary.total_events = stats.total_events;
    stats.event_summary.fraction_active_rois = stats.fraction_rois_with_events;
    if rois_with_events > 0
        stats.event_summary.mean_events_per_active_roi = mean(events_per_roi(events_per_roi > 0));
    else
        stats.event_summary.mean_events_per_active_roi = 0;
    end
    
    % NEW: Clearer field names for peak amplitude access
    stats.all_event_peaks = event_amplitudes;  % CLEAR: All individual peaks flattened across ROIs
    stats.num_total_events = length(event_amplitudes);  % Total number of individual events detected
    
    % Summary statistics for all individual events
    if ~isempty(event_amplitudes)
        stats.mean_all_peaks = mean(event_amplitudes);      % Mean of ALL individual peaks
        stats.median_all_peaks = median(event_amplitudes);  % Median of ALL individual peaks
        stats.std_all_peaks = std(event_amplitudes);        % Standard deviation of ALL peaks
    else
        stats.mean_all_peaks = 0;
        stats.median_all_peaks = 0;
        stats.std_all_peaks = 0;
    end
end