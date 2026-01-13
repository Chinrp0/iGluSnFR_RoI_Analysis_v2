function [event_mask, event_stats] = optimal_schmitt_event_detector(dfof_data, config, baseline_stats)
    % OPTIMAL_SCHMITT_EVENT_DETECTOR - Best implementation using all available toolboxes
    % Uses Parallel Computing, Signal Processing, and Statistics toolboxes for optimal performance
    %
    % Features:
    % - Parallel processing across ROIs
    % - Enhanced noise estimation with iterative refinement
    % - Signal processing for robust filtering
    % - Advanced peak detection with findpeaks()
    % - Comprehensive noise characterization
    
    if nargin < 2
        config = tracenorm_config();
    end
    
    if nargin < 3
        baseline_stats = [];
    end
    
    [numFrames, numROIs] = size(dfof_data);
    
    % Get parameters
    upper_threshold_sigma = config.event_detection.upper_threshold_sigma;
    lower_threshold_sigma = config.event_detection.lower_threshold_sigma;
    decay_extension_frames = config.event_detection.decay_extension_frames;
    rise_extension_frames = config.event_detection.rise_extension_frames;
    
    if config.verbose
        fprintf('Optimal Schmitt trigger detection using all toolboxes...\n');
        fprintf('  Parallel processing: %d workers\n', gcp('nocreate').NumWorkers);
        fprintf('  Thresholds: %.1fσ/%.1fσ\n', upper_threshold_sigma, lower_threshold_sigma);
    end
    
    %% === OPTIMAL NOISE ESTIMATION (Parallel + Statistics + Signal Processing) ===
    tic;
    [upper_threshold, lower_threshold, noise_metrics] = optimal_noise_estimation(...
        dfof_data, upper_threshold_sigma, lower_threshold_sigma, config);
    noise_time = toc;
    
    if config.verbose
        fprintf('  Noise estimation: %.3f s (%.1f μs/ROI)\n', noise_time, noise_time*1e6/numROIs);
        fprintf('  Mean noise: %.4f ± %.4f dF/F\n', ...
            mean(noise_metrics.noise_std, 'omitnan'), std(noise_metrics.noise_std, 'omitnan'));
    end
    
    %% === PARALLEL EVENT DETECTION ===
    tic;
    event_mask = false(numFrames, numROIs);
    
    % Parallel processing across ROIs
    parfor roi = 1:numROIs
        roi_trace = dfof_data(:, roi);
        roi_upper = upper_threshold(roi);
        roi_lower = lower_threshold(roi);
        
        if all(isnan(roi_trace)) || roi_upper <= 0 || isnan(roi_upper)
            continue;
        end
        
        % Enhanced Schmitt trigger with signal processing
        roi_events = enhanced_schmitt_trigger(roi_trace, roi_upper, roi_lower);
        
        % Extend events based on sensor kinetics
        roi_events = extend_events_for_kinetics(roi_events, rise_extension_frames, decay_extension_frames);
        
        event_mask(:, roi) = roi_events;
    end
    detection_time = toc;
    
    if config.verbose
        fprintf('  Event detection: %.3f s (%.1f μs/ROI)\n', detection_time, detection_time*1e6/numROIs);
    end
    
    %% === ENHANCED EVENT CHARACTERIZATION ===
    tic;
    event_stats = enhanced_event_statistics(event_mask, dfof_data, upper_threshold, lower_threshold, ...
        noise_metrics, config);
    stats_time = toc;
    
    % Include event_mask for consistency
    event_stats.event_mask = event_mask;
    event_stats.noise_metrics = noise_metrics;
    
    if config.verbose
        fprintf('  Event analysis: %.3f s\n', stats_time);
        fprintf('  Detected %d total events across %d ROIs\n', ...
            event_stats.total_events, event_stats.rois_with_events);
        fprintf('  Total processing: %.3f s\n', noise_time + detection_time + stats_time);
    end
end

function [upper_threshold, lower_threshold, noise_metrics] = optimal_noise_estimation(...
    dfof_data, upper_sigma, lower_sigma, config)
    % Optimal noise estimation using all three toolboxes
    
    [numFrames, numROIs] = size(dfof_data);
    
    % Initialize outputs
    upper_threshold = zeros(1, numROIs);
    lower_threshold = zeros(1, numROIs);
    
    % Noise characterization arrays
    noise_std = zeros(1, numROIs);
    baseline_points = zeros(1, numROIs);
    convergence_iterations = zeros(1, numROIs);
    noise_method = cell(1, numROIs);
    signal_quality = zeros(1, numROIs);
    
    % Pre-filter data using Signal Processing Toolbox for better noise estimation
    if config.verbose
        fprintf('    Pre-filtering traces for noise estimation...\n');
    end
    
    % Apply gentle median filtering to reduce high-frequency noise in noise estimation
    filtered_data = medfilt1(dfof_data, 3, [], 1, 'truncate');  % 3-point median filter
    
    % Parallel noise estimation across ROIs
    parfor roi = 1:numROIs
        roi_trace = filtered_data(:, roi);
        
        if all(isnan(roi_trace))
            noise_std(roi) = 0.01;
            baseline_points(roi) = 0;
            convergence_iterations(roi) = 0;
            noise_method{roi} = 'invalid';
            signal_quality(roi) = 0;
            continue;
        end
        
        % Method 1: Iterative exclusion with convergence
        [noise_iter, points_iter, iters_iter] = iterative_noise_estimation(roi_trace, upper_sigma);
        
        % Method 2: Robust percentile method (fallback)
        [noise_pct, points_pct] = percentile_noise_estimation(roi_trace, 70);
        
        % Method 3: Hampel filter approach (Signal Processing Toolbox)
        [noise_hampel, points_hampel] = hampel_noise_estimation(roi_trace);
        
        % Choose best method based on convergence and data quality
        if iters_iter < 5 && points_iter > 50  % Good iterative convergence
            noise_std(roi) = noise_iter;
            baseline_points(roi) = points_iter;
            convergence_iterations(roi) = iters_iter;
            noise_method{roi} = 'iterative';
        elseif points_pct > 30  % Sufficient percentile data
            noise_std(roi) = noise_pct;
            baseline_points(roi) = points_pct;
            convergence_iterations(roi) = 0;
            noise_method{roi} = 'percentile';
        else  % Fallback to Hampel
            noise_std(roi) = noise_hampel;
            baseline_points(roi) = points_hampel;
            convergence_iterations(roi) = 0;
            noise_method{roi} = 'hampel';
        end
        
        % Calculate signal quality metric
        signal_range = range(roi_trace, 'omitnan');
        signal_quality(roi) = signal_range / noise_std(roi);  % Signal-to-noise ratio
    end
    
    % Handle invalid noise estimates
    invalid_noise = noise_std <= 0 | isnan(noise_std);
    if any(invalid_noise)
        median_noise = median(noise_std(~invalid_noise));
        if isnan(median_noise) || median_noise <= 0
            median_noise = 0.01;
        end
        noise_std(invalid_noise) = median_noise;
    end
    
    % Calculate final thresholds
    upper_threshold = upper_sigma * noise_std;
    lower_threshold = lower_sigma * noise_std;
    
    % Compile noise metrics
    noise_metrics = struct();
    noise_metrics.noise_std = noise_std;
    noise_metrics.baseline_points = baseline_points;
    noise_metrics.convergence_iterations = convergence_iterations;
    noise_metrics.noise_method = noise_method;
    noise_metrics.signal_quality = signal_quality;
    noise_metrics.upper_thresholds = upper_threshold;
    noise_metrics.lower_thresholds = lower_threshold;
    
    % Summary statistics
    noise_metrics.mean_noise = mean(noise_std, 'omitnan');
    noise_metrics.std_noise = std(noise_std, 'omitnan');
    noise_metrics.mean_snr = mean(signal_quality, 'omitnan');
    noise_metrics.methods_used = unique(noise_method(~strcmp(noise_method, 'invalid')));
end

function [noise_std, baseline_points, iterations] = iterative_noise_estimation(trace, initial_sigma)
    % Iterative noise estimation with convergence checking
    
    max_iterations = 10;
    convergence_threshold = 0.02;  % 2% change threshold
    
    current_noise = mad(trace, 1) * 1.4826;  % Initial estimate
    
    for iter = 1:max_iterations
        % Identify potential events
        threshold = initial_sigma * current_noise;
        potential_events = trace > threshold;
        baseline_periods = ~potential_events;
        
        if sum(baseline_periods) < 20
            break;  % Insufficient baseline data
        end
        
        % Calculate new noise estimate
        baseline_data = trace(baseline_periods);
        new_noise = mad(baseline_data, 1) * 1.4826;
        
        % Check convergence
        if abs(new_noise - current_noise) / current_noise < convergence_threshold
            current_noise = new_noise;
            break;
        end
        
        current_noise = new_noise;
    end
    
    noise_std = current_noise;
    baseline_points = sum(baseline_periods);
    iterations = iter;
end

function [noise_std, baseline_points] = percentile_noise_estimation(trace, percentile)
    % Robust percentile-based noise estimation
    
    threshold = prctile(trace, percentile);
    baseline_periods = trace <= threshold;
    
    if sum(baseline_periods) > 20
        baseline_data = trace(baseline_periods);
        noise_std = mad(baseline_data, 1) * 1.4826;
        baseline_points = length(baseline_data);
    else
        noise_std = mad(trace, 1) * 1.4826;
        baseline_points = length(trace);
    end
end

function [noise_std, baseline_points] = hampel_noise_estimation(trace)
    % Use Hampel identifier to find outliers, estimate noise from non-outliers
    
    try
        % Hampel filter to identify outliers (Signal Processing Toolbox)
        window_size = min(50, floor(length(trace)/4));
        [~, outlier_idx] = hampel(trace, window_size, 3);
        
        baseline_periods = ~outlier_idx;
        
        if sum(baseline_periods) > 20
            baseline_data = trace(baseline_periods);
            noise_std = mad(baseline_data, 1) * 1.4826;
            baseline_points = length(baseline_data);
        else
            noise_std = mad(trace, 1) * 1.4826;
            baseline_points = length(trace);
        end
    catch
        % Fallback if hampel fails
        noise_std = mad(trace, 1) * 1.4826;
        baseline_points = length(trace);
    end
end

function event_mask = enhanced_schmitt_trigger(trace, upper_thresh, lower_thresh)
    % Enhanced Schmitt trigger with signal processing improvements
    
    numFrames = length(trace);
    event_mask = false(numFrames, 1);
    
    % Optional: Apply gentle smoothing for more stable detection
    % smooth_trace = smoothdata(trace, 'movmean', 3);  % Signal Processing Toolbox
    smooth_trace = trace;  % Use original for now
    
    % Standard Schmitt trigger state machine
    state = 'baseline';
    event_start = 0;
    
    for frame = 1:numFrames
        signal = smooth_trace(frame);
        
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
                    event_mask(event_start:frame) = true;
                    state = 'baseline';
                    event_start = 0;
                end
        end
    end
    
    % Handle event extending to end
    if strcmp(state, 'in_event') && event_start > 0
        event_mask(event_start:end) = true;
    end
end

function extended_mask = extend_events_for_kinetics(event_mask, rise_frames, decay_frames)
    % Extend event boundaries (unchanged from original)
    
    numFrames = length(event_mask);
    extended_mask = event_mask;
    
    if ~any(event_mask)
        return;
    end
    
    event_starts = find(diff([false; event_mask]) == 1);
    event_ends = find(diff([event_mask; false]) == -1);
    
    for i = 1:length(event_starts)
        start_frame = event_starts(i);
        end_frame = event_ends(i);
        
        extended_start = max(1, start_frame - rise_frames);
        extended_end = min(numFrames, end_frame + decay_frames);
        
        extended_mask(extended_start:extended_end) = true;
    end
end

function stats = enhanced_event_statistics(event_mask, dfof_data, upper_threshold, lower_threshold, ...
    noise_metrics, config)
    % Enhanced event statistics using Signal Processing Toolbox
    
    [numFrames, numROIs] = size(event_mask);
    
    % Basic counts
    event_frames_per_roi = sum(event_mask, 1);
    rois_with_events = sum(event_frames_per_roi > 0);
    total_event_frames = sum(event_mask, 'all');
    
    % Enhanced event characterization
    event_durations = [];
    event_amplitudes = [];
    event_rise_times = [];
    event_decay_times = [];
    events_per_roi = zeros(1, numROIs);
    
    % Use Signal Processing Toolbox for enhanced peak detection
    for roi = 1:numROIs
        roi_events = event_mask(:, roi);
        roi_dfof = dfof_data(:, roi);
        
        if any(roi_events)
            % Find discrete events
            event_starts = find(diff([false; roi_events]) == 1);
            event_ends = find(diff([roi_events; false]) == -1);
            events_per_roi(roi) = length(event_starts);
            
            % Enhanced characterization for each event
            for e = 1:length(event_starts)
                start_frame = event_starts(e);
                end_frame = event_ends(e);
                event_trace = roi_dfof(start_frame:end_frame);
                
                % Duration
                duration = end_frame - start_frame + 1;
                event_durations(end+1) = duration;
                
                % Peak amplitude using findpeaks for precision
                try
                    [peaks, peak_locs] = findpeaks(event_trace, 'NPeaks', 1, 'SortStr', 'descend');
                    if ~isempty(peaks)
                        event_amplitudes(end+1) = peaks(1);
                        
                        % Rise and decay times
                        peak_idx = peak_locs(1);
                        half_max = peaks(1) / 2;
                        
                        % Rise time (10% to 90% of peak)
                        rise_trace = event_trace(1:peak_idx);
                        rise_90 = find(rise_trace >= 0.9 * peaks(1), 1, 'first');
                        rise_10 = find(rise_trace >= 0.1 * peaks(1), 1, 'first');
                        if ~isempty(rise_90) && ~isempty(rise_10)
                            event_rise_times(end+1) = rise_90 - rise_10;
                        end
                        
                        % Decay time (90% to 10% of peak after peak)
                        if peak_idx < length(event_trace)
                            decay_trace = event_trace(peak_idx:end);
                            decay_90 = find(decay_trace <= 0.9 * peaks(1), 1, 'first');
                            decay_10 = find(decay_trace <= 0.1 * peaks(1), 1, 'first');
                            if ~isempty(decay_90) && ~isempty(decay_10)
                                event_decay_times(end+1) = decay_10 - decay_90;
                            end
                        end
                    else
                        event_amplitudes(end+1) = max(event_trace, [], 'omitnan');
                    end
                catch
                    % Fallback to simple max
                    event_amplitudes(end+1) = max(event_trace, [], 'omitnan');
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
    
    % Enhanced event characteristics
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
    
    % Kinetics (new with Signal Processing)
    if ~isempty(event_rise_times)
        stats.mean_rise_time = mean(event_rise_times);
        stats.median_rise_time = median(event_rise_times);
    else
        stats.mean_rise_time = 0;
        stats.median_rise_time = 0;
    end
    
    if ~isempty(event_decay_times)
        stats.mean_decay_time = mean(event_decay_times);
        stats.median_decay_time = median(event_decay_times);
    else
        stats.mean_decay_time = 0;
        stats.median_decay_time = 0;
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