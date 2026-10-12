function res = twoap_analyzer(dfof_data, baseline_stats, config)
    % TWOAP_ANALYZER - Per-ROI paired-pulse (2AP / PPF) analysis
    %
    % For each ROI computes the response to pulse 1 and pulse 2 and the
    % paired-pulse ratio (PPR = A2 / A1) by THREE methods, so peak-selection can
    % be compared (they diverge most at short ISI, where pulse 2 rides on the
    % decay of pulse 1). In all three, A1 uses the pre-stim1 baseline and A2 is
    % measured above a LOCAL residual baseline taken just before stim2:
    %
    %   1) peak-max      : A1 = max(dF/F) in [0, peak_search] after stim1 (capped
    %                      at the ISI so it can't reach pulse 2). A2 = max(dF/F) in
    %                      [0, peak_search] after stim2. Peak-to-peak.
    %   2) sync time-lock: A1 / A2 read at a FIXED latency (sync_latency_ms, ~10 ms
    %                      = the synchronous-release phase) after each stim. Both
    %                      pulses sampled at the same early time -> apples-to-apples.
    %   3) p1-latency lock: A1 = pulse-1 peak-max (as in method 1); A2 read at
    %                      pulse 1's OWN peak latency after stim2 (each ROI locks
    %                      pulse 2 to where its pulse 1 peaked). Adaptive per ROI.
    %
    % PPR is only defined for pulse-1 responders (A1 above threshold and > 0);
    % otherwise NaN. Thresholds come from each ROI's own pre-stim noise. The
    % stimulus frames come from config.twoAP (stim1_frame, stim2_frame set per
    % file by twoap_pipeline from the filename ISI).
    %
    % Inputs:
    %   dfof_data      - [frames x ROIs] dF/F traces
    %   baseline_stats - struct from prestim_baseline (prestim_idx); optional
    %   config         - tracenorm_config() struct with config.twoAP.stim2_frame /
    %                    isi_ms populated for this file
    %
    % Output:
    %   res - struct of [1 x ROIs] metric arrays plus per-ROI detail

    if nargin < 3 || isempty(config)
        config = tracenorm_config();
    end

    % findpeaks warns when a window has nothing above threshold; silence just that
    % warning (onCleanup restores prior state on return/error), same as 1AP.
    wstate = warning('off', 'signal:findpeaks:largeMinPeakHeight');
    restoreWarn = onCleanup(@() warning(wstate));

    [numFrames, numROIs] = size(dfof_data);
    p          = config.twoAP;
    frame_dur  = config.frame_dur_ms;
    stim1      = p.stim1_frame;
    if ~isfield(p, 'stim2_frame') || isempty(p.stim2_frame)
        error('twoap_analyzer:noStim2', ...
            'config.twoAP.stim2_frame is not set (twoap_pipeline sets it per file).');
    end
    stim2  = p.stim2_frame;
    isi_ms = getfielddef(p, 'isi_ms', (stim2 - stim1) * frame_dur);

    % ms -> frames
    search_frames  = max(1, round(p.peak_search_ms / frame_dur));
    sync_frames    = max(1, round(p.sync_latency_ms / frame_dur));
    isi_frames     = stim2 - stim1;
    % Pulse-1 search must not reach pulse 2: cap at one frame before stim2.
    search1_end    = min(search_frames, max(isi_frames - 1, 1));
    % AUC = TOTAL CHARGE over ONE window encompassing both pulses: [stim1,
    % stim2 + auc_post_frames], integrated above the single pre-stim1 baseline.
    % Not a per-pulse split and not a ratio (see config.twoAP.auc_post_ms).
    auc_post_frames = max(1, round(getfielddef(p, 'auc_post_ms', 250) / frame_dur));

    % Schmitt validation frames (above the LOWER threshold), same scheme as 1AP.
    min_frames_above = 1;
    if isfield(p, 'min_frames_above_lower'), min_frames_above = max(1, p.min_frames_above_lower); end
    % findpeaks min peak spacing (a pulse peak must be a genuine local maximum).
    min_peak_dist = max(1, round(getfielddef(p, 'min_peak_distance_ms', 10) / frame_dur));
    % Noise-ceiling QC: peak1 must exceed peak_over_baseline_k * this percentile of
    % the ROI's own pre-stim1 excursions above baseline (0 disables).
    base_exc_pct  = getfielddef(p, 'baseline_excursion_pct', 99);
    peak_over_k   = getfielddef(p, 'peak_over_baseline_k', 0);

    % Pre-stim index (before pulse 1)
    if isstruct(baseline_stats) && isfield(baseline_stats, 'prestim_idx') ...
            && ~isempty(baseline_stats.prestim_idx)
        prestim_idx = baseline_stats.prestim_idx;
    else
        guard_frames = round(getfielddef(p, 'baseline_guard_ms', 0) / frame_dur);
        prestim_idx  = 1:max(stim1 - 1 - guard_frames, 1);
    end

    % --- Preallocate ------------------------------------------------------
    prestim_sd      = nan(1, numROIs);
    upper_threshold = nan(1, numROIs);
    lower_threshold = nan(1, numROIs);      % Schmitt decay/hysteresis threshold
    base_mean_roi   = nan(1, numROIs);
    residual_base   = nan(1, numROIs);      % local baseline just before stim2
    base_exc_ceil   = nan(1, numROIs);      % pre-stim1 excursion ceiling (noise-ceiling QC)

    % Method 1: peak-max
    peak1_amp = nan(1, numROIs); peak1_time_ms = nan(1, numROIs);
    peak2_amp = nan(1, numROIs); peak2_time_ms = nan(1, numROIs);
    ppr_max   = nan(1, numROIs);
    % Method 2: sync time-lock (fixed sync_latency_ms)
    a1_sync = nan(1, numROIs); a2_sync = nan(1, numROIs); ppr_sync = nan(1, numROIs);
    % Method 3: p1-latency lock (A1 = pulse-1 peak; A2 read at pulse-1 latency)
    a2_p1lock = nan(1, numROIs); ppr_p1lock = nan(1, numROIs);
    % Method 4: total charge (dF/F area over ONE window spanning both pulses)
    q_pair = nan(1, numROIs);

    is_responder1  = false(1, numROIs);     % responded to pulse 1
    is_responder2  = false(1, numROIs);     % responded to pulse 2
    noise_excluded = false(1, numROIs);
    analyzed       = false(1, numROIs);

    % Local baseline window before stim2: the 2 frames immediately before it
    % (best available estimate of the residual level; at very short ISI this
    % sits on the pulse-1 rise and slightly overestimates the residual).
    b2_lo = max(stim2 - 2, stim1 + 1);
    b2_hi = stim2 - 1;

    for roi = 1:numROIs
        trace = dfof_data(:, roi);
        if all(isnan(trace)), continue; end

        pre = trace(prestim_idx);
        pre = pre(isfinite(pre));
        if numel(pre) < 5, continue; end
        base_mean = mean(pre);
        sd = mad(pre, 1) * 1.4826;
        if ~isfinite(sd) || sd <= 0, sd = std(pre); end
        if ~isfinite(sd) || sd <= 0, continue; end

        if isfield(p, 'max_noise_sd') && sd > p.max_noise_sd
            prestim_sd(roi) = sd; noise_excluded(roi) = true; continue;
        end
        prestim_sd(roi)      = sd;
        base_mean_roi(roi)   = base_mean;
        u = base_mean + p.upper_threshold_sigma * sd;
        l = base_mean + getfielddef(p, 'lower_threshold_sigma', 1.5) * sd;
        upper_threshold(roi) = u;
        lower_threshold(roi) = l;
        analyzed(roi)        = true;

        % Noise-ceiling QC: the largest "typical" pre-stim excursion above baseline.
        % A real pulse-1 peak must beat peak_over_k * this (rejects uniform-noise
        % ROIs whose response is no bigger than their own baseline swings).
        bexc = prctile(pre - base_mean, base_exc_pct);
        base_exc_ceil(roi) = bexc;

        % Local residual baseline before pulse 2
        if b2_hi >= b2_lo
            seg = trace(b2_lo:b2_hi); seg = seg(isfinite(seg));
            if ~isempty(seg), residual_base(roi) = mean(seg); end
        end
        b2 = residual_base(roi);
        if ~isfinite(b2), b2 = base_mean; end

        % --- Pulse 1 (local-max + Schmitt, with the noise-ceiling floor) ---
        w1 = stim1:min(stim1 + search1_end, numFrames);
        [p1v, p1i, sustained1] = window_peak(trace, w1, u, l, min_frames_above, min_peak_dist);
        t1_frames = NaN;
        if isfinite(p1v)
            p1_abs = w1(p1i);
            t1_frames          = p1_abs - stim1;                 % pulse-1 peak latency (frames)
            peak1_amp(roi)     = p1v - base_mean;
            peak1_time_ms(roi) = t1_frames * frame_dur;
            % Responder needs: a sustained local max >= upper AND the peak must
            % exceed the ROI's own baseline-excursion ceiling (noise-ceiling QC).
            floor_ok = (peak_over_k <= 0) || ~isfinite(bexc) || (peak1_amp(roi) > peak_over_k * max(bexc, 0));
            is_responder1(roi) = sustained1 && (p1v >= u) && floor_ok;
        end

        % --- Pulse 2 (local-max + Schmitt), amplitude above the local residual ---
        % Negative incremental A2 means pulse 2 did not respond (its peak sits below
        % the pulse-1-elevated residual): clip A2 to 0 and flag non-responder, so the
        % per-ROI PPR is floored at 0 rather than going nonsensically negative.
        w2 = stim2:min(stim2 + search_frames, numFrames);
        [p2v, p2i, sustained2] = window_peak(trace, w2, u, l, min_frames_above, min_peak_dist);
        if isfinite(p2v)
            a2 = p2v - b2;
            peak2_time_ms(roi) = (w2(p2i) - stim2) * frame_dur;
            if a2 < 0
                peak2_amp(roi)     = 0;
                is_responder2(roi) = false;
            else
                peak2_amp(roi)     = a2;
                is_responder2(roi) = sustained2 && (p2v >= u);
            end
        end

        % --- Method 1: peak-max PPR (pulse-1 responders only) ---
        if is_responder1(roi) && peak1_amp(roi) > 0 && isfinite(peak2_amp(roi))
            ppr_max(roi) = peak2_amp(roi) / peak1_amp(roi);
        end

        % --- Method 2: sync time-lock (both pulses at sync_latency_ms) ---
        f1s = stim1 + sync_frames;
        f2s = stim2 + sync_frames;
        if f1s <= numFrames && isfinite(trace(f1s)), a1_sync(roi) = trace(f1s) - base_mean; end
        if f2s <= numFrames && isfinite(trace(f2s)), a2_sync(roi) = trace(f2s) - b2; end
        if is_responder1(roi) && isfinite(a1_sync(roi)) && a1_sync(roi) > 0 && isfinite(a2_sync(roi))
            ppr_sync(roi) = a2_sync(roi) / a1_sync(roi);
        end

        % --- Method 3: p1-latency lock (A1 = pulse-1 peak; A2 at pulse-1 latency) ---
        if isfinite(t1_frames)
            f2p = stim2 + t1_frames;
            if f2p <= numFrames && isfinite(trace(f2p)), a2_p1lock(roi) = trace(f2p) - b2; end
            if is_responder1(roi) && peak1_amp(roi) > 0 && isfinite(a2_p1lock(roi))
                ppr_p1lock(roi) = a2_p1lock(roi) / peak1_amp(roi);
            end
        end

        % --- Method 4: total charge over the encompassing window ---
        % Integrate dF/F above the SINGLE pre-stim1 baseline from stim1 to
        % auc_post_frames after stim2 (one window spanning both pulses). This is
        % the total glutamate charge of the paired response, not a ratio, and it
        % avoids the short-ISI pulse-splitting / residual problems entirely.
        % Computed for every analyzable ROI (no responder gating; the primary
        % figures use the ratio-of-means value on the coverslip average trace).
        wq = stim1:min(stim2 + auc_post_frames, numFrames);
        q_pair(roi) = trapz_ms(trace, wq, base_mean, frame_dur);
    end

    % --- Assemble result --------------------------------------------------
    res = struct();
    res.dataType        = '2AP';
    res.numROIs         = numROIs;
    res.numFrames       = numFrames;
    res.frame_rate      = config.frame_rate;
    res.frame_dur_ms    = frame_dur;
    res.isi_ms          = isi_ms;
    res.stim1_frame     = stim1;
    res.stim2_frame     = stim2;
    res.sync_latency_ms = p.sync_latency_ms;
    res.prestim_idx     = prestim_idx;

    res.prestim_sd      = prestim_sd;
    res.upper_threshold = upper_threshold;
    res.lower_threshold = lower_threshold;
    res.prestim_base    = base_mean_roi;
    res.residual_base   = residual_base;
    res.base_exc_ceil   = base_exc_ceil;   % pre-stim1 excursion ceiling (noise-ceiling QC)

    res.peak1_amp       = peak1_amp;      % A1 (dF/F above pre-stim1 baseline)
    res.peak1_time_ms   = peak1_time_ms;
    res.peak2_amp       = peak2_amp;      % A2 (dF/F above local pre-stim2 residual)
    res.peak2_time_ms   = peak2_time_ms;
    res.ppr_max         = ppr_max;        % method 1: peak-max

    res.a1_sync         = a1_sync;
    res.a2_sync         = a2_sync;
    res.ppr_sync        = ppr_sync;       % method 2: sync time-lock (10 ms)

    res.a2_p1lock       = a2_p1lock;
    res.ppr_p1lock      = ppr_p1lock;     % method 3: p1-latency lock

    res.q_pair          = q_pair;         % method 4: total charge (dF/F*ms) over
                                          % [stim1, stim2 + auc_post_ms], above pre-stim1 baseline
    res.auc_post_ms     = auc_post_frames * frame_dur;

    res.is_responder1   = is_responder1;
    res.is_responder2   = is_responder2;
    res.noise_excluded  = noise_excluded;

    res.n_responders1   = sum(is_responder1);
    res.n_ppr_valid_max    = sum(isfinite(ppr_max));
    res.n_ppr_valid_sync   = sum(isfinite(ppr_sync));
    res.n_ppr_valid_p1lock = sum(isfinite(ppr_p1lock));
    res.n_analyzable    = sum(analyzed);
    res.responder1_fraction = res.n_responders1 / max(numROIs, 1);
end

function [pv, pidx, sustained] = window_peak(trace, win, u, l, min_frames_above, min_peak_dist)
    % Genuine local-maximum detection with Schmitt dual-threshold validation, the
    % SAME scheme as oneap_analyzer (findpeaks). A pulse peak must be a LOCAL
    % MAXIMUM that clears the UPPER threshold u (findpeaks, MinPeakHeight=u,
    % MinPeakDistance=min_peak_dist), after which the run of consecutive frames
    % >= the LOWER threshold l spanning it must be >= min_frames_above. Requiring a
    % true local max (not just the window max) rejects monotonic noise ramps that
    % happen to end high. Keeps fast transients whose decay lingers above the lower.
    % Returns the peak value, its index WITHIN win, and the sustained flag; when no
    % qualifying local max exists it returns the window max with sustained=false, so
    % the amplitude is still measurable but the ROI fails the responder gate (caller
    % additionally requires pv >= u).
    pv = NaN; pidx = NaN; sustained = false;
    if isempty(win), return; end
    if nargin < 6 || isempty(min_peak_dist), min_peak_dist = 1; end
    seg = trace(win);
    fin = isfinite(seg);
    if ~any(fin), return; end
    seg(~fin) = min(seg(fin));      % pin NaNs below everything (they can't be peaks)

    pv_all = []; loc_all = [];
    if any(seg >= u) && numel(seg) >= 3
        mpd = max(1, min(min_peak_dist, numel(seg) - 1));
        [pv_all, loc_all] = findpeaks(seg, 'MinPeakHeight', u, 'MinPeakDistance', mpd);
    end
    if isempty(pv_all)
        % No interior local max above threshold -> report the window max for the
        % amplitude estimate, but leave sustained=false so it fails the gate.
        [pv, pidx] = max(seg);
        return;
    end
    [pv, k] = max(pv_all); pidx = loc_all(k);   % largest above-threshold local max
    % Count the consecutive >= l run around the peak (in absolute frames)
    li = win(pidx); ri = win(pidx);
    while li > 1 && isfinite(trace(li-1)) && trace(li-1) >= l, li = li - 1; end
    while ri < numel(trace) && isfinite(trace(ri+1)) && trace(ri+1) >= l, ri = ri + 1; end
    sustained = (ri - li + 1) >= min_frames_above;
end

function a = trapz_ms(trace, win, base, frame_dur)
    % Trapezoidal integral of (trace - base) over frame window win, in dF/F*ms.
    % NaN frames are dropped (rare); needs >= 2 finite samples.
    a = NaN;
    if numel(win) < 2, return; end
    seg = trace(win) - base;
    seg = seg(isfinite(seg));
    if numel(seg) < 2, return; end
    a = trapz(seg) * frame_dur;
end

function v = getfielddef(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
