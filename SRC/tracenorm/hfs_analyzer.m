function res = hfs_analyzer(dfof_data, baseline_stats, config)
    % HFS_ANALYZER - Per-ROI high-frequency train (HFS) analysis
    %
    % For each ROI: pre-stim noise QC, a responder call, and the per-pulse and
    % train-level metrics of hfs_train_metrics (phasic amplitude and tonic level
    % per pulse, train charge, plateau, steady-state ratio, cumulative-amplitude
    % pool / replenishment). Single-ROI single-trial traces are noisy, so these
    % are per-ROI diagnostics; the primary values are measured on the coverslip
    % grand-average trace in hfs_batch_analysis (no responder gating).
    %
    % Responder: the trace smoothed over responder_smooth_ms must, somewhere in
    % the train, clear the baseline by upper_threshold_sigma x the SD of the
    % equally smoothed pre-stim trace. A train response is sustained, so the
    % smoothed comparison is fairer than a single-frame peak.
    %
    % Inputs:
    %   dfof_data      - [frames x ROIs] dF/F epoch (stim1 at config.hfs.stim1_frame)
    %   baseline_stats - struct from prestim_baseline (prestim_idx, outlier_mask)
    %   config         - tracenorm_config() with apply_train_timing applied
    %
    % Output:
    %   res - struct of [1 x ROIs] metrics, [nPulses x ROIs] per-pulse matrices

    p  = config.hfs;
    fd = config.frame_dur_ms;
    [numFrames, numROIs] = size(dfof_data);
    s1   = p.stim1_frame;
    nP   = numel(p.pulse_frames);
    tend = min(p.train_end_frame, numFrames);

    if isstruct(baseline_stats) && isfield(baseline_stats, 'prestim_idx') && ~isempty(baseline_stats.prestim_idx)
        prestim_idx = baseline_stats.prestim_idx;
    else
        prestim_idx = 1:max(s1 - 1 - round(p.baseline_guard_ms / fd), 1);
    end
    pre_outliers = prestim_outlier_mask(baseline_stats, numel(prestim_idx), numROIs);
    sm = max(1, round(p.responder_smooth_ms / fd));

    prestim_sd = nan(1, numROIs); prestim_base = nan(1, numROIs);
    noise_excluded = false(1, numROIs); analyzed = false(1, numROIs);
    is_responder = false(1, numROIs);
    names = {'q_train', 'plateau', 'peak_train', 'a1', 'ss_amp', 'ss_tonic', 'ss_ratio', ...
        'pool_est', 'repl_per_pulse', 'median_lat_ms'};
    for k = 1:numel(names), out.(names{k}) = nan(1, numROIs); end
    amp = nan(nP, numROIs); tonic = nan(nP, numROIs);

    for roi = 1:numROIs
        trace = double(dfof_data(:, roi));
        pre_raw = trace(prestim_idx);
        pre = pre_raw(isfinite(pre_raw));
        if numel(pre) < 5, continue; end
        base = prestim_base_mean(pre_raw, pre_outliers(:, roi));
        sd = mad(pre, 1) * 1.4826;
        if ~isfinite(sd) || sd <= 0, sd = std(pre); end
        if ~isfinite(sd) || sd <= 0, continue; end
        prestim_sd(roi) = sd;
        if sd > p.max_noise_sd, noise_excluded(roi) = true; continue; end
        prestim_base(roi) = base;
        analyzed(roi) = true;

        % Responder on the smoothed trace, against the smoothed pre-stim noise.
        smt = movmean(trace, sm, 'omitnan');
        sd_s = std(smt(prestim_idx), 'omitnan');
        resp = max(smt(s1:tend - 1)) - base;
        is_responder(roi) = isfinite(resp) && sd_s > 0 && resp >= p.upper_threshold_sigma * sd_s;

        m = hfs_train_metrics(trace, base, config, false);
        for k = 1:numel(names), out.(names{k})(roi) = m.(names{k}); end
        amp(:, roi) = m.amp(:); tonic(:, roi) = m.tonic(:);
    end

    res = out;
    res.dataType       = 'HFS';
    res.numROIs        = numROIs;
    res.numFrames      = numFrames;
    res.frame_dur_ms   = fd;
    res.stim1_frame    = s1;
    res.pulse_frames   = p.pulse_frames;
    res.train_end_frame = p.train_end_frame;
    res.resolvable     = hfs_pulses_resolvable(config);
    res.prestim_idx    = prestim_idx;
    res.prestim_sd     = prestim_sd;
    res.prestim_base   = prestim_base;
    res.noise_excluded = noise_excluded;
    res.analyzed       = analyzed;            % passed noise QC (used for the averages)
    res.is_responder   = is_responder;
    res.amp            = amp;                 % [nPulses x ROIs] phasic amplitude
    res.tonic          = tonic;               % [nPulses x ROIs] pre-pulse level above baseline
    res.n_analyzable   = sum(analyzed);
    res.n_responders   = sum(is_responder);
    res.responder_fraction = res.n_responders / max(res.n_analyzable, 1);
end
