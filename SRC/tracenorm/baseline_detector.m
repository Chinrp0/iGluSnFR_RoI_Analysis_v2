function [baseline, outlier_mask, stats] = baseline_detector(data, config)
    % BASELINE_DETECTOR - Iterative rolling median baseline calculation
    % Implements the algorithm from iGluSnFR3 paper (peakFinder.R equivalent)
    %
    % Inputs:
    %   data   - [frames x ROIs] fluorescence intensity matrix
    %   config - Configuration struct from tracenorm_config()
    %
    % Outputs:
    %   baseline     - [frames x ROIs] baseline fluorescence estimate
    %   outlier_mask - [frames x ROIs] logical mask of detected outliers
    %   stats        - Struct with baseline calculation statistics
    
    if nargin < 2
        config = tracenorm_config();
    end
    
    [numFrames, numROIs] = size(data);
    
    % Validate input
    if numFrames ~= config.expected_frames
        warning('Expected %d frames, got %d', config.expected_frames, numFrames);
    end
    
    if config.verbose
        fprintf('Computing iterative rolling median baseline for %d ROIs...\n', numROIs);
    end
    
    % === Stage 1: Initial outlier detection ===
    if config.verbose, fprintf('  Stage 1: Initial rolling median and outlier detection\n'); end
    
    % Calculate initial rolling median across all ROIs (vectorized)
    initial_median = movmedian(data, config.rolling_window_frames, 1, 'omitnan');
    
    % Calculate residuals and standard deviation for outlier detection
    residuals = data - initial_median;
    
    % Calculate rolling standard deviation (approximate using MAD * 1.4826)
    rolling_mad = movmedian(abs(residuals), config.rolling_window_frames, 1, 'omitnan');
    rolling_std = rolling_mad * 1.4826;  % Convert MAD to std approximation
    
    % Flag outliers: points above threshold*sigma from rolling median
    outlier_mask = residuals > (config.outlier_threshold_sigma * rolling_std);
    
    % === Stage 2: Refined baseline excluding outliers ===
    if config.verbose, fprintf('  Stage 2: Refined baseline excluding outliers\n'); end
    
    % Create data with outliers set to NaN
    data_no_outliers = data;
    data_no_outliers(outlier_mask) = NaN;
    
    % Recalculate rolling median excluding outliers
    refined_baseline = movmedian(data_no_outliers, config.rolling_window_frames, 1, 'omitnan');
    
    % === Stage 3: Fill outlier gaps with "last observation carried forward" ===
    if config.verbose, fprintf('  Stage 3: Filling outlier gaps\n'); end
    
    % Apply LOCF to fill NaN values in the refined baseline
    baseline = fillmissing(refined_baseline, 'previous', 1);
    
    % Handle any remaining NaN at the beginning (forward fill)
    baseline = fillmissing(baseline, 'next', 1);
    
    % === Optional: Additional refinement iterations ===
    if config.max_iterations > 3
        for iter = 4:config.max_iterations
            if config.verbose, fprintf('  Iteration %d: Further refinement\n', iter); end
            
            % Recalculate outliers with new baseline
            new_residuals = data - baseline;
            new_rolling_mad = movmedian(abs(new_residuals), config.rolling_window_frames, 1, 'omitnan');
            new_rolling_std = new_rolling_mad * 1.4826;
            outlier_mask = new_residuals > (config.outlier_threshold_sigma * new_rolling_std);
            
            % Update baseline
            data_no_outliers = data;
            data_no_outliers(outlier_mask) = NaN;
            refined_baseline = movmedian(data_no_outliers, config.rolling_window_frames, 1, 'omitnan');
            baseline = fillmissing(refined_baseline, 'previous', 1);
            baseline = fillmissing(baseline, 'next', 1);
        end
    end
    
    % === Optional: Smooth final baseline ===
    if config.smooth_baseline
        switch config.smooth_method
            case 'movmean'
                baseline = movmean(baseline, config.smooth_window, 1, 'omitnan');
            case 'movmedian'
                baseline = movmedian(baseline, config.smooth_window, 1, 'omitnan');
            case 'gaussian'
                % Apply Gaussian smoothing using imgaussfilt1
                for roi = 1:numROIs
                    baseline(:, roi) = imgaussfilt1(baseline(:, roi), config.smooth_window/4);
                end
        end
    end
    
    % === Calculate statistics ===
    stats = calculate_baseline_stats(data, baseline, outlier_mask, config);
    
    if config.verbose
        fprintf('  Completed: %.1f%% outliers detected, %.1f%% valid baseline points\n', ...
            100 * sum(outlier_mask(:)) / numel(outlier_mask), ...
            100 * stats.mean_valid_fraction);
    end
end

function stats = calculate_baseline_stats(data, baseline, outlier_mask, config)
    % Calculate baseline quality statistics
    
    [numFrames, numROIs] = size(data);
    
    % Per-ROI statistics
    outlier_fraction = sum(outlier_mask, 1) / numFrames;  % [1 x ROIs]
    valid_fraction = 1 - outlier_fraction;
    
    % Baseline stability (coefficient of variation)
    baseline_mean = mean(baseline, 1, 'omitnan');
    baseline_std = std(baseline, 0, 1, 'omitnan');
    baseline_cv = baseline_std ./ baseline_mean;
    
    % Signal-to-baseline ratio (max signal / baseline)
    max_signal = max(data, [], 1);
    signal_to_baseline = max_signal ./ baseline_mean;
    
    % Transport detection (linear trend in baseline)
    time_vector = (1:numFrames)';
    transport_slopes = zeros(1, numROIs);
    for roi = 1:numROIs
        p = polyfit(time_vector, baseline(:, roi), 1);
        transport_slopes(roi) = p(1);  % Slope of linear fit
    end
    
    % Identify potential transport ROIs
    transport_rois = transport_slopes > config.transport_slope_threshold;
    
    % Compile statistics
    stats = struct();
    stats.outlier_fraction = outlier_fraction;
    stats.valid_fraction = valid_fraction;
    stats.mean_outlier_fraction = mean(outlier_fraction);
    stats.mean_valid_fraction = mean(valid_fraction);
    
    stats.baseline_cv = baseline_cv;
    stats.signal_to_baseline = signal_to_baseline;
    stats.transport_slopes = transport_slopes;
    stats.transport_rois = transport_rois;
    stats.num_transport_rois = sum(transport_rois);
    
    stats.quality_flags = struct();
    stats.quality_flags.low_valid_data = valid_fraction < (config.min_baseline_frames / size(data, 1));
    stats.quality_flags.high_noise = baseline_cv > 0.2;  % >20% CV indicates high noise
    stats.quality_flags.potential_transport = transport_rois;
    
    % Overall quality summary
    stats.num_low_quality = sum(stats.quality_flags.low_valid_data | ...
                               stats.quality_flags.high_noise);
    stats.fraction_good_rois = 1 - (stats.num_low_quality / numROIs);
end