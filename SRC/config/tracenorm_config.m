function config = tracenorm_config()
    % TRACENORM_CONFIG - Configuration for fluorescence trace normalization pipeline
    % CORRECTED: 3000 frames, 10ms exposure (100 Hz), 30s total recording
    %
    % Returns:
    %   config - struct with all pipeline parameters
    
    config = struct();
    
    %% === CORRECTED: Imaging Parameters ===
    % 10ms exposure = 100 Hz frame rate
    % 3000 frames ÷ 100 Hz = 30 seconds total recording
    config.frame_rate = 40;                   % Hz - 10ms exposure = 100 fps
    config.expected_frames = 1200;             % Total frames in recording
    config.recording_duration_s = config.expected_frames / config.frame_rate; % 30 seconds
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
    fprintf('Frame Rate: %d Hz (%.1f ms exposure)\n', config.frame_rate, config.temporal.frame_duration_ms);
    fprintf('Recording: %d frames = %.1f seconds\n', config.expected_frames, config.recording_duration_s);
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