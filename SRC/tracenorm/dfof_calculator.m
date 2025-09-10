function [dfof_data, dfof_stats] = dfof_calculator(raw_data, baseline, config)
    % DFOF_CALCULATOR - Calculate dF/F traces from raw data and baseline
    % Implements: dF/F = (F - F0) / F0 where F0 is the baseline
    %
    % Inputs:
    %   raw_data - [frames x ROIs] raw fluorescence intensity
    %   baseline - [frames x ROIs] baseline fluorescence from baseline_detector
    %   config   - Configuration struct from tracenorm_config()
    %
    % Outputs:
    %   dfof_data  - [frames x ROIs] normalized dF/F traces
    %   dfof_stats - Struct with dF/F calculation statistics
    
    if nargin < 3
        config = tracenorm_config();
    end
    
    [numFrames, numROIs] = size(raw_data);
    
    if config.verbose
        fprintf('Calculating dF/F for %d ROIs...\n', numROIs);
    end
    
    % Validate inputs
    if ~isequal(size(raw_data), size(baseline))
        error('Raw data and baseline must have same dimensions');
    end
    
    % === Calculate dF/F ===
    switch config.dfof_method
        case 'divide'
            % Standard: dF/F = (F - F0) / F0
            dfof_data = (raw_data - baseline) ./ baseline;
            
        case 'subtract'
            % Alternative: dF/F = F - F0 (absolute change)
            dfof_data = raw_data - baseline;
            
        case 'percent'
            % Percentage change: dF/F = 100 * (F - F0) / F0
            dfof_data = 100 * (raw_data - baseline) ./ baseline;
            
        otherwise
            error('Unknown dF/F method: %s', config.dfof_method);
    end
    
    % Handle division by zero or very small baselines
    invalid_baseline = baseline <= 0 | baseline < (0.01 * mean(baseline, 'all', 'omitnan'));
    dfof_data(invalid_baseline) = NaN;
    
    % === Calculate dF/F Statistics ===
    dfof_stats = calculate_dfof_stats(raw_data, baseline, dfof_data, config);
    
    if config.verbose
        fprintf('  dF/F range: [%.3f, %.3f], %.1f%% valid values\n', ...
            min(dfof_data, [], 'all', 'omitnan'), ...
            max(dfof_data, [], 'all', 'omitnan'), ...
            100 * dfof_stats.fraction_valid);
    end
end

function stats = calculate_dfof_stats(raw_data, baseline, dfof_data, config)
    % Calculate comprehensive dF/F statistics
    
    [numFrames, numROIs] = size(dfof_data);
    
    % === Basic Statistics ===
    % Per-ROI statistics
    dfof_mean = mean(dfof_data, 1, 'omitnan');        % [1 x ROIs]
    dfof_std = std(dfof_data, 0, 1, 'omitnan');       % [1 x ROIs] 
    dfof_max = max(dfof_data, [], 1, 'omitnan');      % [1 x ROIs]
    dfof_min = min(dfof_data, [], 1, 'omitnan');      % [1 x ROIs]
    
    % Baseline statistics
    baseline_mean = mean(baseline, 1, 'omitnan');
    baseline_std = std(baseline, 0, 1, 'omitnan');
    baseline_stability = baseline_std ./ baseline_mean;  % Coefficient of variation
    
    % === Event Detection Metrics ===
    % Simple peak detection for validation
    event_threshold = 2 * dfof_std;  % 2 sigma above noise
    potential_events = dfof_data > event_threshold;
    events_per_roi = sum(potential_events, 1);
    
    % === Signal Quality Metrics ===
    % Signal-to-noise ratio
    signal_power = dfof_std;  % Use std as proxy for signal power
    noise_floor = baseline_stability .* baseline_mean;  % Baseline noise
    snr = signal_power ./ noise_floor;
    
    % Dynamic range
    dynamic_range = dfof_max - dfof_min;
    
    % === Data Validity ===
    valid_data = ~isnan(dfof_data) & ~isinf(dfof_data);
    fraction_valid_per_roi = sum(valid_data, 1) / numFrames;
    
    % === Compile Statistics ===
    stats = struct();
    
    % Per-ROI metrics
    stats.dfof_mean = dfof_mean;
    stats.dfof_std = dfof_std;
    stats.dfof_max = dfof_max;
    stats.dfof_min = dfof_min;
    stats.dynamic_range = dynamic_range;
    stats.snr = snr;
    stats.events_per_roi = events_per_roi;
    stats.fraction_valid_per_roi = fraction_valid_per_roi;
    
    % Baseline metrics
    stats.baseline_mean = baseline_mean;
    stats.baseline_stability = baseline_stability;
    
    % Overall metrics
    stats.overall_mean = mean(dfof_mean, 'omitnan');
    stats.overall_std = mean(dfof_std, 'omitnan');
    stats.overall_max = max(dfof_max, [], 'omitnan');
    stats.overall_min = min(dfof_min, [], 'omitnan');
    stats.fraction_valid = sum(valid_data, 'all') / numel(dfof_data);
    stats.mean_events_per_roi = mean(events_per_roi, 'omitnan');
    stats.mean_snr = mean(snr, 'omitnan');
    
    % === Quality Flags ===
    stats.quality_flags = struct();
    stats.quality_flags.low_snr = snr < 2;                           % SNR < 2
    stats.quality_flags.unstable_baseline = baseline_stability > 0.2; % CV > 20%
    stats.quality_flags.low_signal = dynamic_range < 0.1;            % Small dynamic range
    stats.quality_flags.invalid_data = fraction_valid_per_roi < 0.8;  % <80% valid data
    
    % Summary quality metrics
    stats.num_low_quality = sum(stats.quality_flags.low_snr | ...
                               stats.quality_flags.unstable_baseline | ...
                               stats.quality_flags.low_signal | ...
                               stats.quality_flags.invalid_data);
    stats.fraction_good_rois = 1 - (stats.num_low_quality / numROIs);
    
    % === Event Detection Summary ===
    stats.event_summary = struct();
    stats.event_summary.rois_with_events = sum(events_per_roi > 0);
    stats.event_summary.total_events = sum(events_per_roi);
    stats.event_summary.mean_events_per_active_roi = mean(events_per_roi(events_per_roi > 0));
    stats.event_summary.fraction_active_rois = stats.event_summary.rois_with_events / numROIs;
end