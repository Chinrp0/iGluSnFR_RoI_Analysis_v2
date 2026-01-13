function [dfof_data, basic_stats] = dfof_calculator(raw_data, baseline, config)
    % DFOF_CALCULATOR - Calculate dF/F traces from raw data and baseline
    % REFACTORED VERSION: Only handles dF/F calculation, no event detection
    %
    % Inputs:
    %   raw_data - [frames x ROIs] raw fluorescence intensity
    %   baseline - [frames x ROIs] baseline fluorescence from baseline_detector
    %   config   - Configuration struct from tracenorm_config()
    %
    % Outputs:
    %   dfof_data    - [frames x ROIs] normalized dF/F traces
    %   basic_stats  - Basic dF/F statistics (no events, no quality assessment)
    
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
    
    %% === Calculate dF/F ===
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
    
    %% === Handle Invalid Baselines ===
    % Set very small or negative baselines to NaN
    invalid_baseline = baseline <= 0 | baseline < (0.01 * mean(baseline, 'all', 'omitnan'));
    dfof_data(invalid_baseline) = NaN;
    
    %% === Calculate Basic Statistics Only ===
    basic_stats = calculate_basic_dfof_stats(raw_data, baseline, dfof_data);
    
    if config.verbose
        fprintf('  dF/F range: [%.3f, %.3f], %.1f%% valid values\n', ...
            min(dfof_data, [], 'all', 'omitnan'), ...
            max(dfof_data, [], 'all', 'omitnan'), ...
            100 * basic_stats.fraction_valid);
    end
end

function stats = calculate_basic_dfof_stats(raw_data, baseline, dfof_data)
    % Calculate basic dF/F statistics only (no events, no quality flags)
    
    [numFrames, numROIs] = size(dfof_data);
    
    %% === Basic Statistics ===
    % Per-ROI statistics
    dfof_mean = mean(dfof_data, 1, 'omitnan');        % [1 x ROIs]
    dfof_std = std(dfof_data, 0, 1, 'omitnan');       % [1 x ROIs] 
    dfof_max = max(dfof_data, [], 1, 'omitnan');      % [1 x ROIs]
    dfof_min = min(dfof_data, [], 1, 'omitnan');      % [1 x ROIs]
    
    % Baseline statistics
    baseline_mean = mean(baseline, 1, 'omitnan');
    baseline_std = std(baseline, 0, 1, 'omitnan');
    baseline_stability = baseline_std ./ baseline_mean;  % Coefficient of variation
    
    %% === Signal Quality Metrics ===
    % Signal-to-noise ratio (corrected calculation)
    peak_response = max(dfof_data, [], 1);
    median_dfof = median(dfof_data, 1);
    signal_amplitude = peak_response - median_dfof;
    noise_estimate = mad(dfof_data, 1, 1) * 1.4826;  % Robust noise estimate
    snr = signal_amplitude ./ noise_estimate;
    
    % Dynamic range
    dynamic_range = dfof_max - dfof_min;
    
    %% === Data Validity ===
    valid_data = ~isnan(dfof_data) & ~isinf(dfof_data);
    fraction_valid_per_roi = sum(valid_data, 1) / numFrames;
    
    %% === Compile Basic Statistics ===
    stats = struct();
    
    % Per-ROI metrics
    stats.dfof_mean = dfof_mean;
    stats.dfof_std = dfof_std;
    stats.dfof_max = dfof_max;
    stats.dfof_min = dfof_min;
    stats.dynamic_range = dynamic_range;
    stats.snr = snr;
    stats.signal_amplitude = signal_amplitude;
    stats.noise_estimate = noise_estimate;
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
    stats.mean_snr = mean(snr, 'omitnan');
    
    % NOTE: No event detection or quality flags - these are handled by separate modules
end