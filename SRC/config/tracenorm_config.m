function config = tracenorm_config()
    % TRACENORM_CONFIG - Configuration for baseline calculation and dF/F
    % Implements iterative rolling median algorithm from iGluSnFR3 paper
    
    config = struct();
    
    % === Baseline Calculation Parameters ===
    config.baseline_method = 'iterative_rolling_median';
    config.frame_rate = 40;                    % Hz - imaging frequency
    config.rolling_window_sec = 0.75;          % Rolling window in seconds
    config.rolling_window_frames = round(config.rolling_window_sec * config.frame_rate); % 30 frames
    config.outlier_threshold_sigma = 2.0;      % Standard deviations for outlier detection
    config.max_iterations = 3;                 % Refinement iterations (paper uses 3 stages)
    
    % === Data Validation ===
    config.expected_frames = 1200;             % Expected number of frames
    config.min_valid_frames = 1000;            % Minimum frames required for baseline calc
    
    % === dF/F Calculation ===
    config.dfof_method = 'divide';             % (F - F0) / F0
    config.smooth_baseline = false;            % Apply smoothing to final baseline
    config.smooth_method = 'movmean';          % 'movmean', 'movmedian', 'gaussian'
    config.smooth_window = 5;                  % Smoothing window size
    
    % === Quality Control ===
    config.detect_transport_rois = true;       % Flag ROIs with gradual increase
    config.transport_slope_threshold = 0.1;    % Slope threshold for transport detection
    config.min_baseline_frames = 100;          % Minimum non-outlier frames for valid baseline
    
    % === Visualization ===
    config.plot_sample_traces = true;          % Show sample traces in baseline plots
    config.num_sample_traces = 6;              % Number of traces to plot
    config.plot_validation = true;             % Show baseline validation plots
    
    % === Performance ===
    config.use_parallel = false;               % Rolling median is already optimized
    config.memory_efficient = true;            % Use memory-efficient processing
    
    % === Debug Options ===
    config.verbose = false;                    % Print detailed progress
    config.save_intermediate = false;          % Save intermediate baseline estimates
    
    fprintf('Baseline config: %.1fs window (%d frames), %.1fσ threshold\n', ...
        config.rolling_window_sec, config.rolling_window_frames, config.outlier_threshold_sigma);
end