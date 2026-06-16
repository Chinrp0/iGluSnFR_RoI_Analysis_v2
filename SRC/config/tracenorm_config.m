function config = tracenorm_config()
    % TRACENORM_CONFIG - Configuration for fluorescence trace normalization pipeline
    % CORRECTED: 3000 frames, 10ms exposure (100 Hz), 30s total recording
    %
    % Returns:
    %   config - struct with all pipeline parameters
    
    config = struct();

    %% === Analysis Mode ===
    % Choose the type of experiment to analyze:
    %   'Spont' - spontaneous activity (rolling-median baseline, whole-trace
    %             Schmitt event detection, WT vs R213W comparison)
    %   '1AP'   - single action potential evoked responses (pre-stimulus
    %             baseline, evoked-peak + sync/async + AUC + tau analysis)
    % This is the single switch that selects the whole analysis path.
    config.dataType = '1AP';                 % '1AP' or 'Spont'

    %% === CORRECTED: Imaging Parameters ===
    % Frame rate is the SINGLE acquisition knob: all 1AP time windows below are
    % specified in milliseconds and converted to frames at run time, so they
    % follow frame_rate automatically. Only the two experiment-specific frame
    % counts (oneAP.expected_frames, oneAP.stim_frame) are set per dataset.
    %   100 Hz -> 10 ms/frame      200 Hz -> 5 ms/frame
    config.frame_rate = 200;                   % Hz (1000/frame_rate = ms per frame)

    %% === 1AP (evoked single action potential) Parameters ===
    % All windows are in MILLISECONDS RELATIVE TO STIMULUS ONSET (frame-rate
    % independent). Only stim_frame and expected_frames are in frames because
    % they are properties of the specific recording, not the frame rate.
    config.oneAP = struct();
    config.oneAP.expected_frames   = 600;      % Total frames per 1AP trace (this dataset)
    config.oneAP.stim_frame        = 266;      % Frame of stimulus onset (PROGRAMMABLE, this dataset)
    config.oneAP.sync_window_ms    = [0, 15];  % Synchronous release window (timing label only)
    config.oneAP.async_window_ms   = [15, 250];% Asynchronous release window (timing label only)
                                               % NOTE: the evoked iGlu peak lands ~15-35 ms post-stim
                                               % (median 25 ms, data-measured), so the responder call
                                               % and amplitude use the LARGEST peak, not the sync window.
    config.oneAP.bin_size_ms       = 5;       % Temporal histogram bin
    config.oneAP.auc_window_ms     = [0, 250]; % Window for area-under-curve (baseline-subtracted)
    % Decay kinetics have NO fixed fit-window parameter:
    %   - averaged-responder fit (DOUBLE exponential) runs from the peak to the
    %     end of the trace, so the slow component is fully captured
    %   - per-ROI fit (mono) runs from the peak until the trace returns below the
    %     lower threshold, i.e. it spans that ROI's own decay
    config.oneAP.tau_peak_window_ms = 50;      % Decay-onset peak searched in [0, this] post-stim
    config.oneAP.tau_min_r2        = 0.5;      % Min fit R^2 to keep a decay tau
    config.oneAP.baseline_guard_ms = 0;        % Time before stim to drop from baseline (settling)
    % Multi-peak detection thresholds for 1AP (sigma of pre-stim noise)
    config.oneAP.upper_threshold_sigma = 3.0;  % Event / responder threshold (peak height)
    config.oneAP.lower_threshold_sigma = 1.5;  % Lower threshold (reported for reference)
    config.oneAP.min_peak_distance_ms  = 10;   % Min spacing between detected peaks (findpeaks)
    % Option B responder cleaning (data-driven, this dataset):
    config.oneAP.min_frames_above_upper = 2;   % Peak must stay >= upper threshold for >= this many
                                               % frames. Rejects single-frame noise spikes; the real
                                               % discriminator (sustained ABOVE the upper threshold,
                                               % not the lower one). Set to 1 for old height-only behavior.
    config.oneAP.max_noise_sd          = 0.02; % Drop ROIs whose pre-stim dF/F robust SD exceeds this
                                               % (too noisy to score). median sigma ~0.006, p90 ~0.009.

    %% === Mode-dependent frame counts ===
    if strcmpi(config.dataType, '1AP')
        config.expected_frames = config.oneAP.expected_frames;  % 450
    else
        config.expected_frames = 3000;         % Total frames in recording
    end
    config.recording_duration_s = config.expected_frames / config.frame_rate;
    config.min_valid_frames = 400;             % Minimum frames required for baseline calc
    
    %% === Baseline Calculation Parameters ===
    config.baseline_method = 'iterative_rolling_median';
    config.rolling_window_sec = 0.50;          % Rolling window in seconds (500ms)
    config.rolling_window_frames = round(config.rolling_window_sec * config.frame_rate); % 50 frames at 100 Hz
    config.outlier_threshold_sigma = 2.0;      % Standard deviations for outlier detection
    config.max_iterations = 3;                 % Refinement iterations
    
    % NEW: Missing fields required by baseline_detector.m
    config.verbose = false;                    % Suppress detailed console output
    config.smooth_baseline = false;            % Don't smooth final baseline (already using rolling median)
    config.smooth_method = 'movmedian';        % Method if smoothing enabled
    config.smooth_window = 5;                  % Window size for smoothing (frames)
    config.transport_slope_threshold = 0.01;   % Threshold for detecting transport ROIs
    config.min_baseline_frames = 400;          % Minimum frames for valid baseline
    
    %% === dF/F Calculation ===
    config.dfof_method = 'divide';           % (F - F0) / F0
    config.min_baseline_value = 0.01;          % Avoid division by zero
    
    %% === CORRECTED Schmitt Trigger Event Detection ===
    config.corrected_schmitt = struct();
    config.corrected_schmitt.upper_threshold_sigma = 3.5;  % Upper threshold (event start)
    config.corrected_schmitt.lower_threshold_sigma = 1.5;  % Lower threshold (event end)
    config.corrected_schmitt.min_event_duration = 3;       % Minimum 3 frames (30ms at 100 Hz)
    config.corrected_schmitt.merge_gap_frames = 2;         % Merge events ≤2 frames apart (20ms)
    config.corrected_schmitt.noise_exclusion_window = 7;   % Exclude sustained events ≥7 frames (70ms)
    config.corrected_schmitt.sustained_percentile = 85;    % ADD THIS LINE - 85th percentile threshold
    config.corrected_schmitt.includes_outliers_in_noise = true;
    config.corrected_schmitt.excludes_sustained_events = true;
        
    %% === Legacy Schmitt Trigger (for comparison) ===
    config.event_detection = struct();
    config.event_detection.method = 'schmitt_trigger';
    config.event_detection.upper_threshold_sigma = 3.0;    % Legacy: lower threshold
    config.event_detection.lower_threshold_sigma = 1.5;
    config.event_detection.min_event_duration = 3;         % frames
    config.event_detection.merge_gap_frames = 2;           % frames
    
    %% === Quality Metrics ===
    config.quality = struct();
    config.quality.min_baseline_stability = 0.8;     % Minimum stability score
    config.quality.max_noise_level = 0.2;            % Maximum acceptable noise (dF/F)
    config.quality.min_snr = 2.0;                    % Minimum signal-to-noise ratio
    
    %% === Temporal Parameters (CORRECTED for 100 Hz) ===
    config.temporal = struct();
    config.temporal.frame_rate = config.frame_rate;         % 100 Hz
    config.temporal.frame_duration_ms = 1000 / config.frame_rate; % 10 ms per frame
    config.temporal.recording_duration_s = config.recording_duration_s; % 30 seconds
    config.temporal.max_iei_s = 10;                         % Maximum IEI to analyze (10s)
    config.temporal.bin_size_ms = 50;                       % For IEI histograms (50ms bins)
    
    %% === Visualization ===
    config.visualization = struct();
    config.visualization.plot_baseline = true;
    config.visualization.plot_events = true;
    config.visualization.max_traces_to_plot = 100;    % Limit for raster plots
    config.visualization.colors = struct(...
        'wt', [0.2, 0.6, 1.0], ...               % Blue for WT
        'mut', [1.0, 0.4, 0.2]);                 % Red/Orange for mutant
    
    %% === Output Options ===
    config.output = struct();
    config.output.save_figures = false;
    config.output.figure_format = 'png';
    config.output.figure_dpi = 300;
    config.output.save_data = false;
    
    %% === Summary Display ===
    fprintf('\n=== Configuration Loaded ===\n');
    fprintf('Analysis Mode: %s\n', config.dataType);
    fprintf('Frame Rate: %d Hz (%.1f ms exposure)\n', config.frame_rate, config.temporal.frame_duration_ms);
    fprintf('Recording: %d frames = %.1f seconds\n', config.expected_frames, config.recording_duration_s);
    if strcmpi(config.dataType, '1AP')
        fprintf('1AP Stim Frame: %d (%.0f ms), Sync: %d-%d ms, Async: %d-%d ms, Bin: %d ms\n', ...
            config.oneAP.stim_frame, (config.oneAP.stim_frame - 1) * config.temporal.frame_duration_ms, ...
            config.oneAP.sync_window_ms(1), config.oneAP.sync_window_ms(2), ...
            config.oneAP.async_window_ms(1), config.oneAP.async_window_ms(2), ...
            config.oneAP.bin_size_ms);
    end
    fprintf('Baseline Window: %.0f ms (%d frames)\n', ...
        config.rolling_window_sec * 1000, config.rolling_window_frames);
    fprintf('Event Detection: Corrected Schmitt %.1f/%.1fσ\n', ...
        config.corrected_schmitt.upper_threshold_sigma, ...
        config.corrected_schmitt.lower_threshold_sigma);
    fprintf('Min Event Duration: %d frames (%.0f ms)\n', ...
        config.corrected_schmitt.min_event_duration, ...
        config.corrected_schmitt.min_event_duration * config.temporal.frame_duration_ms);
    fprintf('============================\n\n');
end