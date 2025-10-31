function config = tracenorm_config()
    % TRACENORM_CONFIG - Updated configuration for corrected Schmitt trigger pipeline
    % Includes new corrected_schmitt parameters for proper noise estimation
    
    config = struct();
    
    %% === Baseline Calculation Parameters ===
    config.baseline_method = 'iterative_rolling_median';
    config.frame_rate = 10;                    % Hz - imaging frequency
    config.rolling_window_sec = 0.50;          % Rolling window in seconds
    config.rolling_window_frames = round(config.rolling_window_sec * config.frame_rate); % 30 frames
    config.outlier_threshold_sigma = 2.0;      % Standard deviations for outlier detection
    config.max_iterations = 3;                 % Refinement iterations
    
    %% === Data Validation ===
    config.expected_frames = 3000;             % Expected number of frames
    config.min_valid_frames = 400;            % Minimum frames required for baseline calc
    
    %% === dF/F Calculation ===
    config.dfof_method = 'divide';             % (F - F0) / F0
    config.smooth_baseline = false;            % Apply smoothing to final baseline
    config.smooth_method = 'movmean';          % 'movmean', 'movmedian', 'gaussian'
    config.smooth_window = 5;                  % Smoothing window size
    
    %% === CORRECTED SCHMITT TRIGGER - NEW IMPLEMENTATION ===
    config.corrected_schmitt = struct();
    config.corrected_schmitt.upper_threshold_sigma = 3.5;      % Upper threshold (start events)
    config.corrected_schmitt.lower_threshold_sigma = 1.5;      % Lower threshold (end events)
    config.corrected_schmitt.noise_exclusion_window = 7;       % Frames to exclude from noise calc
    config.corrected_schmitt.min_event_duration = 3;           % Minimum event duration
    config.corrected_schmitt.sustained_percentile = 85;        % Percentile for sustained event detection
    config.corrected_schmitt.use_outlier_mask = true;          % Use baseline detector's outlier mask
    config.corrected_schmitt.include_outliers_in_noise = true; % CRITICAL: Include outliers in noise calc
    
    %% === LEGACY EVENT DETECTION - Schmitt Trigger Parameters ===
    config.event_detection = struct();
    config.event_detection.method = 'schmitt_trigger';         % Detection method
    config.event_detection.upper_threshold_sigma = 3.0;        % Upper threshold (event start)
    config.event_detection.lower_threshold_sigma = 1.5;        % Lower threshold (event end)
    config.event_detection.decay_extension_frames = 4;         % Extend events by 3 frames (75ms decay)
    config.event_detection.rise_extension_frames = 1;          % Extend events before start (conservative)
    
    % Legacy event detection (simple threshold) - for comparison
    config.simple_event_threshold_sigma = 2.0;  % Simple threshold for comparison
    
    %% === QUALITY CRITERIA ===
    config.quality_criteria = struct();
    config.quality_criteria.min_snr = 3.0;                     % Minimum signal-to-noise ratio
    config.quality_criteria.max_baseline_cv = 0.2;             % Maximum baseline CV (20%)
    config.quality_criteria.min_valid_fraction = 0.8;          % Minimum valid baseline data (80%)
    config.quality_criteria.min_dynamic_range = 0.01;          % Minimum dF/F range (1%)
    config.quality_criteria.max_transport_slope = 0.1;         % Transport detection threshold
    config.quality_criteria.min_dfof_valid = 0.8;              % Minimum valid dF/F data (80%)
    config.quality_criteria.min_events_for_active = 1;         % Minimum events for "active" ROI
    config.quality_criteria.max_event_fraction = 0.5;          % Max fraction of frames as events (50%)
    
    %% === Transport ROI Detection ===
    config.detect_transport_rois = true;       % Flag ROIs with gradual increase
    config.transport_slope_threshold = config.quality_criteria.max_transport_slope;
    config.min_baseline_frames = 100;          % Minimum non-outlier frames for valid baseline
    
    %% === VISUALIZATION - Fixed for 2x4 layout ===
    config.plot_sample_traces = true;          % Show sample traces in baseline plots
    config.num_sample_traces = 8;              % FIXED: Always 8 traces for 2x4 layout
    config.plot_validation = true;             % Show baseline validation plots
    config.plot_events = true;                 % Show event detection plots
    config.plot_layout = struct();             % NEW: Standardized plot layout
    config.plot_layout.rows = 4;               % FIXED: 4 rows
    config.plot_layout.cols = 2;               % FIXED: 2 columns
    config.plot_layout.subplots = 8;           % FIXED: 8 total subplots
    
    %% === Peak Marker Settings ===
    config.peak_markers = struct();            % NEW: Peak marker configuration
    config.peak_markers.show_above_events = true;      % Show markers above event peaks
    config.peak_markers.height_offset = 0.02;          % Offset above peak (fraction of range)
    config.peak_markers.size = 60;                     % Marker size
    config.peak_markers.color = [0 0.8 0];            % Green color for peak markers
    config.peak_markers.symbol = '^';                  % Triangle symbol
    
    %% === Performance ===
    config.use_parallel = false;               % Rolling median is already optimized
    config.memory_efficient = true;            % Use memory-efficient processing
    
    %% === Debug Options ===
    config.verbose = false;                    % Print detailed progress
    config.save_intermediate = false;          % Save intermediate baseline estimates
    config.compare_event_methods = false;      % Compare Schmitt vs simple threshold
    
    %% === Module Integration Settings ===
    config.ensure_data_consistency = true;     % NEW: Ensure consistent data passing
    config.validate_event_mask = true;         % NEW: Validate event_mask in event_stats
    config.backward_compatibility = true;      % NEW: Maintain old field names for compatibility
    
    %% === Error Handling ===
    config.continue_on_module_error = true;    % Continue pipeline if individual modules fail
    config.create_fallback_data = true;        % Create empty data structures on failure
    
    if config.verbose
        fprintf('Config loaded: Corrected Schmitt trigger pipeline\n');
        fprintf('  Corrected Schmitt: %.1f/%.1fσ thresholds (includes outliers in noise)\n', ...
            config.corrected_schmitt.upper_threshold_sigma, ...
            config.corrected_schmitt.lower_threshold_sigma);
        fprintf('  Legacy Schmitt: %.1f/%.1fσ thresholds\n', ...
            config.event_detection.upper_threshold_sigma, ...
            config.event_detection.lower_threshold_sigma);
        fprintf('  Visualization: %dx%d layout (%d subplots)\n', ...
            config.plot_layout.rows, config.plot_layout.cols, config.plot_layout.subplots);
    end
end