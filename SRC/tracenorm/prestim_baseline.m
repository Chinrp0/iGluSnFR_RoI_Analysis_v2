function [baseline, baseline_stats] = prestim_baseline(raw_data, config)
    % PRESTIM_BASELINE - Per-ROI baseline from the pre-stimulus period (1AP mode)
    %
    % For evoked (1AP) recordings the baseline is the resting fluorescence
    % BEFORE the stimulus. We use the flat pre-stim segment (frames
    % 1 .. stim_frame-1-guard) to estimate, per ROI:
    %   - F0       (mean of pre-stim raw fluorescence, excluding upward
    %               outliers by the same rolling-median / MAD rule as
    %               baseline_detector)                  -> used for dF/F
    %   - noise SD (std of pre-stim raw fluorescence)
    %
    % The returned baseline is a [frames x ROIs] matrix with each column held
    % constant at that ROI's F0, so it plugs directly into dfof_calculator.
    %
    % Inputs:
    %   raw_data - [frames x ROIs] raw fluorescence intensity
    %   config   - tracenorm_config() struct (uses config.oneAP)
    %
    % Outputs:
    %   baseline       - [frames x ROIs] constant pre-stim F0 per column
    %   baseline_stats - struct with per-ROI pre-stim statistics

    if nargin < 2
        config = tracenorm_config();
    end

    [numFrames, numROIs] = size(raw_data);

    % The pre-stim window ends just before the FIRST stimulus. Pick the stim
    % frame / guard from the active paradigm so 1AP, 2AP (PPF) and HFS share this
    % code (2AP and HFS baseline off the first pulse, stim1_frame).
    switch upper(config.dataType)
        case '2AP'
            stim_frame = config.twoAP.stim1_frame;
            guard_ms   = config.twoAP.baseline_guard_ms;
        case 'HFS'
            stim_frame = config.hfs.stim1_frame;
            guard_ms   = config.hfs.baseline_guard_ms;
        otherwise   % '1AP'
            stim_frame = config.oneAP.stim_frame;
            guard_ms   = config.oneAP.baseline_guard_ms;
    end
    frame_dur = config.frame_dur_ms;                             % ms per frame
    guard     = round(guard_ms / frame_dur);

    % Pre-stimulus frame window: 1 .. (stim_frame - 1 - guard)
    prestim_end = stim_frame - 1 - guard;
    if prestim_end < 5
        error('prestim_baseline:tooFewFrames', ...
            'Pre-stim window has only %d frames (stim_frame=%d, guard=%d). Need >=5.', ...
            prestim_end, stim_frame, guard);
    end
    if prestim_end > numFrames
        prestim_end = numFrames;
    end
    prestim_idx = 1:prestim_end;

    prestim = raw_data(prestim_idx, :);                 % [prestimFrames x ROIs]

    % --- Outlier exclusion (same rule as baseline_detector, spont path) ---
    % Rolling median over rolling_window_sec; local sigma = rolling MAD * 1.4826;
    % frames more than outlier_threshold_sigma ABOVE the rolling median (one-
    % sided: stray spontaneous events) are excluded from F0.
    win = config.rolling_window_frames;
    if ~isfinite(win) || win < 1
        win = round(config.rolling_window_sec * 1000 / frame_dur);
    end
    win = min(win, prestim_end);
    roll_med  = movmedian(prestim, win, 1, 'omitnan');
    residuals = prestim - roll_med;
    roll_sd   = movmedian(abs(residuals), win, 1, 'omitnan') * 1.4826;
    outlier_mask = residuals > (config.outlier_threshold_sigma * roll_sd);

    prestim_clean = prestim;
    prestim_clean(outlier_mask) = NaN;

    % Per-ROI pre-stim statistics
    f0_mean_all = mean(prestim, 1, 'omitnan');          % [1 x ROIs] incl. outliers
    f0_mean   = mean(prestim_clean, 1, 'omitnan');      % [1 x ROIs] outliers excluded
    f0_median = median(prestim, 1, 'omitnan');          % [1 x ROIs]
    raw_sd    = std(prestim, 0, 1, 'omitnan');          % [1 x ROIs] classic SD
    % Robust noise estimate (MAD->SD) - resistant to a stray pre-stim event
    robust_sd = mad(prestim, 1, 1) * 1.4826;            % [1 x ROIs]

    % F0 for dF/F = mean of the non-outlier pre-stim frames
    f0 = f0_mean;

    % Build [frames x ROIs] constant baseline for dfof_calculator
    baseline = repmat(f0, numFrames, 1);

    % Guard against zero / negative baselines breaking dF/F downstream
    bad_f0 = ~isfinite(f0) | f0 <= 0;
    if any(bad_f0)
        % Leave them as-is; dfof_calculator will NaN out invalid baselines.
    end

    baseline_stats = struct();
    baseline_stats.method            = 'prestim_baseline';
    baseline_stats.stim_frame        = stim_frame;
    baseline_stats.prestim_frames    = prestim_end;
    baseline_stats.prestim_idx       = prestim_idx;
    baseline_stats.f0_mean           = f0_mean;
    baseline_stats.f0_mean_all       = f0_mean_all;
    baseline_stats.f0_median         = f0_median;
    baseline_stats.outlier_mask      = outlier_mask;                  % [prestimFrames x ROIs]
    baseline_stats.n_outlier_frames  = sum(outlier_mask, 1);          % [1 x ROIs]
    baseline_stats.outlier_fraction  = mean(outlier_mask, 1);         % [1 x ROIs]
    baseline_stats.prestim_raw_sd    = raw_sd;
    baseline_stats.prestim_robust_sd = robust_sd;
    baseline_stats.fraction_good_rois = sum(~bad_f0) / max(numROIs, 1);

    if isfield(config, 'verbose') && config.verbose
        fprintf('      Pre-stim baseline: %d frames (1..%d), %d/%d ROIs valid, %.1f%% outlier frames excluded\n', ...
            prestim_end, prestim_end, sum(~bad_f0), numROIs, 100 * mean(outlier_mask(:)));
    end
end
