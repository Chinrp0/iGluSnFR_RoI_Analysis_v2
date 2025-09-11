function config = tracenorm_config()
    % TRACENORM_CONFIG - Updated configuration with Schmitt trigger event detection
    % Implements modular event detection and quality assessment
    
    config = struct();
    
    %% === Baseline Calculation Parameters ===
    config.baseline_method = 'iterative_rolling_median';
    config.frame_rate = 40;                    % Hz - imaging frequency
    config.rolling_window_sec = 0.75;          % Rolling window in seconds
    config.rolling_window_frames = round(config.rolling_window_sec * config.frame_rate); % 30 frames
    config.outlier_threshold_sigma = 2.0;      % Standard deviations for outlier detection
    config.max_iterations = 3;                 % Refinement iterations (paper uses 3 stages)
    
    %% === Data Validation ===
    config.expected_frames = 1200;             % Expected number of frames
    config.min_valid_frames = 1000;            % Minimum frames required for baseline calc
    
    %% === dF/F Calculation ===
    config.dfof_method = 'divide';             % (F - F0) / F0
    config.smooth_baseline = false;            % Apply smoothing to final baseline
    config.smooth_method = 'movmean';          % 'movmean', 'movmedian', 'gaussian'
    config.smooth_window = 5;                  % Smoothing window size
    
    %% === EVENT DETECTION - Schmitt Trigger Parameters ===
    config.event_detection = struct();
    config.event_detection.method = 'schmitt_trigger';         % Detection method
    config.event_detection.upper_threshold_sigma = 3.0;        % Upper threshold (event start)
    config.event_detection.lower_threshold_sigma = 1.5;        % Lower threshold (event end)
    config.event_detection.decay_extension_frames = 3;         % Extend events by 3 frames (75ms decay)
    config.event_detection.rise_extension_frames = 1;          % Extend events before start (conservative)
    
    % Legacy event detection (simple threshold) - for comparison
    config.simple_event_threshold_sigma = 2.0;  % Simple threshold for comparison
    
    %% === QUALITY CRITERIA ===
    config.quality_criteria = struct();
    config.quality_criteria.min_snr = 2.0;                     % Minimum signal-to-noise ratio
    config.quality_criteria.max_baseline_cv = 0.2;             % Maximum baseline CV (20%)
    config.quality_criteria.min_valid_fraction = 0.8;          % Minimum valid baseline data (80%)
    config.quality_criteria.min_dynamic_range = 0.01;          % Minimum dF/F range (1%)
    config.quality_criteria.max_transport_slope = 0.1;         % Transport detection threshold
    config.quality_criteria.min_dfof_valid = 0.8;              % Minimum valid dF/F data (80%)
    config.quality_criteria.min_events_for_active = 2;         % Minimum events for "active" ROI
    config.quality_criteria.max_event_fraction = 0.5;          % Max fraction of frames as events (50%)
    
    %% === Transport ROI Detection ===
    config.detect_transport_rois = true;       % Flag ROIs with gradual increase
    config.transport_slope_threshold = config.quality_criteria.max_transport_slope;
    config.min_baseline_frames = 100;          % Minimum non-outlier frames for valid baseline
    
    %% === Visualization ===
    config.plot_sample_traces = true;          % Show sample traces in baseline plots
    config.num_sample_traces = 6;              % Number of traces to plot
    config.plot_validation = true;             % Show baseline validation plots
    config.plot_events = true;                 % Show event detection plots
    
    %% === Performance ===
    config.use_parallel = false;               % Rolling median is already optimized
    config.memory_efficient = true;            % Use memory-efficient processing
    
    %% === Debug Options ===
    config.verbose = false;                    % Print detailed progress
    config.save_intermediate = false;          % Save intermediate baseline estimates
    config.compare_event_methods = false;      % Compare Schmitt vs simple threshold
    
    fprintf('Config loaded: Schmitt trigger (%.1f/%.1fσ), %d-frame decay extension\n', ...
        config.event_detection.upper_threshold_sigma, ...
        config.event_detection.lower_threshold_sigma, ...
        config.event_detection.decay_extension_frames);
end