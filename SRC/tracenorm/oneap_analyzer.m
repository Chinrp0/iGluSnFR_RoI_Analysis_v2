function res = oneap_analyzer(dfof_data, baseline_stats, config)
    % ONEAP_ANALYZER - Per-ROI analysis of single action potential (1AP) traces
    %
    % For each ROI (column of dfof_data) this computes, relative to the stimulus
    % frame (config.oneAP.stim_frame):
    %   - ALL peaks above threshold (multi-peak local-maxima detection), each
    %     classified as synchronous or asynchronous by its post-stim peak time
    %   - synchronous vs asynchronous peak counts + a post-stim time histogram
    %   - evoked peak amplitude / time-to-peak (the largest peak above threshold,
    %     wherever it falls post-stim; sync/async is only a timing label of it)
    %   - area under the curve (AUC, baseline-subtracted dF/F)
    %   - decay time constant tau (mono-exponential fit from the peak) + fit R^2
    %
    % Detection uses the Signal Processing Toolbox findpeaks with Schmitt
    % dual-threshold validation: a local maximum that clears the UPPER threshold
    % (base_mean + upper_threshold_sigma * pre-stim SD) is counted only if the
    % response then stays above the LOWER threshold (lower_threshold_sigma) for at
    % least min_frames_above_lower consecutive frames - i.e. it shows a real decay.
    % This keeps fast real transients that clear the upper for only one frame while
    % rejecting single-frame noise spikes. Subject to a minimum peak-to-peak
    % spacing; there is no event merging, so sustained peaks that are far enough
    % apart are counted individually.
    %
    % Thresholds are derived from the PRE-STIMULUS portion of each dF/F trace,
    % so each ROI is judged against its own resting noise.
    %
    % Inputs:
    %   dfof_data      - [frames x ROIs] dF/F traces (aligned: same stim frame)
    %   baseline_stats - struct from prestim_baseline (for prestim_idx); optional
    %   config         - tracenorm_config() struct (uses config.oneAP)
    %
    % Output:
    %   res - struct of [1 x ROIs] metric arrays plus per-ROI detail

    if nargin < 3 || isempty(config)
        config = tracenorm_config();
    end

    % findpeaks warns when an ROI has nothing above threshold (a subthreshold
    % ROI we handle as zero peaks). Silence just that warning for this function;
    % onCleanup restores the prior state on return or error.
    wstate = warning('off', 'signal:findpeaks:largeMinPeakHeight');
    restoreWarn = onCleanup(@() warning(wstate));

    [numFrames, numROIs] = size(dfof_data);

    p          = config.oneAP;
    stim_frame = p.stim_frame;
    frame_dur  = 1000 / config.frame_rate;          % ms per frame

    % Derive frame counts from the ms-based config (frame_rate is the only knob)
    guard_frames        = round(p.baseline_guard_ms / frame_dur);
    min_peak_dist_frames = max(1, round(p.min_peak_distance_ms / frame_dur));

    % --- Frame windows (relative to stim) ---------------------------------
    if isstruct(baseline_stats) && isfield(baseline_stats, 'prestim_idx') ...
            && ~isempty(baseline_stats.prestim_idx)
        prestim_idx = baseline_stats.prestim_idx;
    else
        prestim_end = stim_frame - 1 - guard_frames;
        prestim_idx = 1:max(prestim_end, 1);
    end

    ms2frame = @(t) stim_frame + round(t / frame_dur);

    auc_f    = [max(ms2frame(p.auc_window_ms(1)), 1),  min(ms2frame(p.auc_window_ms(2)), numFrames)];
    % AUC split by release timing: synchronous vs asynchronous windows (the same
    % windows used for the event-count split). With the default config these tile
    % the total AUC window, so auc_sync + auc_async ~= auc.
    sync_auc_f  = [max(ms2frame(p.sync_window_ms(1)), 1),  min(ms2frame(p.sync_window_ms(2)), numFrames)];
    async_auc_f = [max(ms2frame(p.async_window_ms(1)), 1), min(ms2frame(p.async_window_ms(2)), numFrames)];
    post_end = min(ms2frame(p.async_window_ms(2)), numFrames);   % peak search end
    post_idx = stim_frame:post_end;

    % Temporal histogram (ms, post-stim) -----------------------------------
    % Centre the bins ON the frame grid (0, bin_size, 2*bin_size, ...) so each
    % bin sits on a real possible event time. Edges fall at the half-bin points
    % (-bin/2, +bin/2, ...), so events - which land exactly on frame times -
    % never sit on an edge and are counted at the matching centre.
    half_bin       = p.bin_size_ms / 2;
    bin_centers_ms = 0 : p.bin_size_ms : p.async_window_ms(2);
    bin_edges_ms   = [bin_centers_ms - half_bin, bin_centers_ms(end) + half_bin];
    nBins = numel(bin_centers_ms);

    % --- Preallocate outputs ----------------------------------------------
    peak_amplitude  = nan(1, numROIs);
    peak_frame      = nan(1, numROIs);
    time_to_peak_ms = nan(1, numROIs);
    is_responder    = false(1, numROIs);
    prestim_sd      = nan(1, numROIs);
    base_mean_roi   = nan(1, numROIs);   % per-ROI pre-stim mean (for averaging)
    upper_threshold = nan(1, numROIs);
    lower_threshold = nan(1, numROIs);
    auc             = nan(1, numROIs);
    auc_sync        = nan(1, numROIs);   % AUC over synchronous window
    auc_async       = nan(1, numROIs);   % AUC over asynchronous window (slow tail)
    tau_ms          = nan(1, numROIs);
    tau_r2          = nan(1, numROIs);
    tau_A           = nan(1, numROIs);   % fitted amplitude  (for plotting)
    tau_offset      = nan(1, numROIs);   % fitted offset + baseline (dF/F units)
    tau_peak_frame  = nan(1, numROIs);   % decay-onset peak frame (within tau window)
    sync_count      = zeros(1, numROIs);
    async_count     = zeros(1, numROIs);
    total_events    = zeros(1, numROIs);
    temporal_hist   = zeros(nBins, numROIs);
    event_peak_times_ms = cell(1, numROIs);

    % QC bookkeeping: did the trace ever clear the upper (3.5 sigma) gate, and
    % was it analyzable at all (enough finite pre-stim, valid noise SD)?
    crossed_upper   = false(1, numROIs);
    analyzed        = false(1, numROIs);
    noise_excluded  = false(1, numROIs);   % dropped: pre-stim noise > max_noise_sd

    sync_lo = p.sync_window_ms(1);  sync_hi = p.sync_window_ms(2);
    async_lo = p.async_window_ms(1); async_hi = p.async_window_ms(2);

    for roi = 1:numROIs
        trace = dfof_data(:, roi);

        if all(isnan(trace))
            continue;
        end

        % --- Pre-stim noise -> thresholds (dF/F units) ---
        pre = trace(prestim_idx);
        pre = pre(isfinite(pre));
        if numel(pre) < 5
            continue;
        end
        base_mean = mean(pre);
        sd = mad(pre, 1) * 1.4826;          % robust SD
        if ~isfinite(sd) || sd <= 0
            sd = std(pre);
        end
        if ~isfinite(sd) || sd <= 0
            continue;
        end
        % ROI-selection noise ceiling: an ROI whose resting noise exceeds
        % max_noise_sd is too noisy to score reliably (its 3.5 sigma gate sits
        % so high that real responses can't clear it and noise excursions can).
        % Record it separately so it is counted, not silently dropped.
        if isfield(p, 'max_noise_sd') && sd > p.max_noise_sd
            prestim_sd(roi)    = sd;
            noise_excluded(roi) = true;
            continue;
        end
        prestim_sd(roi)      = sd;
        base_mean_roi(roi)   = base_mean;
        u = base_mean + p.upper_threshold_sigma * sd;
        l = base_mean + p.lower_threshold_sigma * sd;
        upper_threshold(roi) = u;
        lower_threshold(roi) = l;
        analyzed(roi)        = true;   % cleared NaN / pre-stim / SD checks

        % --- Multi-peak detection over the post-stim window ---
        % Local maxima >= upper threshold are detected, then validated below by a
        % Schmitt dual-threshold check: after crossing the UPPER (3σ) threshold the
        % response must stay above the LOWER (1.5σ) threshold for
        % min_frames_above_lower frames. NaNs are pinned to the baseline so they
        % cannot form spurious peaks.
        post_trace = trace(post_idx);
        post_trace(~isfinite(post_trace)) = base_mean;
        crossed_upper(roi) = any(post_trace > u);   % cleared threshold at all?

        % Signal Processing Toolbox multi-peak detector. Guard on crossed_upper
        % so findpeaks is skipped (and stays quiet) when nothing clears threshold.
        if crossed_upper(roi)
            [pk_vals, pk_locs] = findpeaks(post_trace, ...
                'MinPeakHeight', u, 'MinPeakDistance', min_peak_dist_frames);
        else
            pk_vals = []; pk_locs = [];
        end

        min_frames_above = 1;
        if isfield(p, 'min_frames_above_lower')
            min_frames_above = max(1, p.min_frames_above_lower);
        end

        ev_peak_times = [];                  % only validated peaks are recorded
        best_peak_val = -Inf;
        for e = 1:numel(pk_locs)
            % Schmitt dual-threshold validation: the peak has already crossed the
            % UPPER threshold (findpeaks MinPeakHeight = u). A real evoked transient
            % then decays THROUGH the lower threshold over several frames, while a
            % single-frame noise spike drops straight back to baseline. Count the
            % run of consecutive frames >= the LOWER threshold spanning this peak
            % and reject peaks that don't stay elevated - so fast real transients
            % that clear the upper for only one frame are KEPT (their decay lingers
            % above the lower threshold), while single-frame fliers are rejected.
            li = pk_locs(e);
            while li > 1 && post_trace(li-1) >= l, li = li - 1; end
            ri = pk_locs(e);
            while ri < numel(post_trace) && post_trace(ri+1) >= l, ri = ri + 1; end
            if (ri - li + 1) < min_frames_above
                continue;                    % not sustained above lower -> not counted
            end

            pk_abs = post_idx(1) + pk_locs(e) - 1;       % absolute frame
            t_ms   = (pk_abs - stim_frame) * frame_dur;  % time from stim
            ev_peak_times(end+1) = t_ms; %#ok<AGROW>

            % Classify each peak by timing (synchronous vs asynchronous release).
            % This is a label only - it does NOT decide responder status.
            if t_ms >= sync_lo && t_ms <= sync_hi
                sync_count(roi) = sync_count(roi) + 1;
            elseif t_ms > async_lo && t_ms <= async_hi
                async_count(roi) = async_count(roi) + 1;
            end

            % The evoked response = the LARGEST sustained peak above threshold,
            % wherever it falls in the post-stim window. Amplitude, time-to-peak
            % and the responder call all key off this peak (the indicator-
            % convolved evoked peak lands ~15-35 ms post-stim, outside a narrow
            % synchronous window, so restricting to the sync window would measure
            % the rising edge instead of the response).
            if pk_vals(e) > best_peak_val
                best_peak_val        = pk_vals(e);
                peak_amplitude(roi)  = pk_vals(e) - base_mean;
                peak_frame(roi)      = pk_abs;
                time_to_peak_ms(roi) = t_ms;
            end
        end
        event_peak_times_ms{roi} = ev_peak_times;
        total_events(roi) = numel(ev_peak_times);

        % Responder = at least one sustained evoked peak above threshold (the
        % largest of which is the measured response). Synchronous vs asynchronous
        % is a timing sub-label of that peak, tracked separately in sync_count /
        % async_count - not a gate on responder status.
        is_responder(roi) = total_events(roi) >= 1;
        tau_onset = peak_frame(roi);       % decay starts at the largest evoked peak
        tau_peak_frame(roi) = peak_frame(roi);

        % Temporal histogram of post-stim peak times
        if ~isempty(ev_peak_times)
            temporal_hist(:, roi) = histcounts(ev_peak_times, bin_edges_ms)';
        end

        % --- AUC (baseline-subtracted dF/F), total + sync / async split ---
        % (dF/F)*ms over each window. The sync/async split separates the fast
        % evoked transient from the slow asynchronous-release tail.
        auc(roi)       = auc_window(trace, auc_f,       base_mean, stim_frame, frame_dur);
        auc_sync(roi)  = auc_window(trace, sync_auc_f,  base_mean, stim_frame, frame_dur);
        auc_async(roi) = auc_window(trace, async_auc_f, base_mean, stim_frame, frame_dur);

        % --- Decay tau (mono-exponential from the peak) -----------------------
        % No fixed fit-window parameter: the window spans this ROI's own decay,
        % from the peak until the trace returns below the lower threshold (with a
        % short margin). A single exponential fit over the full trace is unstable
        % (the long noisy tail collapses the R^2), so we fit only the decay phase.
        if isfinite(tau_onset)
            tau_end_frame = tau_onset;
            while tau_end_frame < numFrames && trace(tau_end_frame + 1) >= l
                tau_end_frame = tau_end_frame + 1;
            end
            tau_end_frame = min(tau_end_frame + 2, numFrames);   % small tail margin
            if tau_end_frame > tau_onset + 2
                dec_idx = tau_onset:tau_end_frame;
                y = trace(dec_idx) - base_mean;             % decay toward 0
                t = (0:numel(dec_idx)-1)' * frame_dur;      % ms from peak
                win_ms = t(end);                            % data-driven fit span
                [tau_ms(roi), tau_r2(roi), Afit, Cfit] = fit_decay_tau(t, y);

                % Discard unreliable fits so noise-dominated decays (near-
                % threshold responders) don't pile up near 0 ms in the tau
                % plots. Keep tau only when it is resolvable (>= 1 frame, <= the
                % fit span) and the fit is acceptable (R^2 >= tau_min_r2).
                % tau_r2 is left intact so the QC R^2 histogram still shows
                % every fit's quality.
                if tau_ms(roi) >= frame_dur && tau_ms(roi) <= win_ms ...
                        && tau_r2(roi) >= p.tau_min_r2
                    % Store the curve in absolute dF/F units for plotting:
                    %   yhat(t) = tau_offset + tau_A * exp(-t / tau_ms)
                    tau_A(roi)      = Afit;
                    tau_offset(roi) = base_mean + Cfit;
                else
                    tau_ms(roi) = NaN;
                end
            end
        end
    end

    % --- Response-type QC classification ----------------------------------
    % Codes (index qc_category_names by code+1). Codes 1 and 2 are BOTH
    % responders (is_responder = any sustained evoked peak); the split is by the
    % timing of the LARGEST peak (synchronous vs asynchronous release):
    %   0 excluded        - not analyzable (all-NaN / <5 pre-stim / bad SD)
    %   1 responder_sync  - responder whose largest peak is in the sync window
    %   2 responder_async - responder whose largest peak is in the async window
    %   3 subthreshold    - analyzable but no peak cleared the threshold
    %   4 noise_excluded  - dropped by the pre-stim noise ceiling (max_noise_sd)
    peak_is_sync = is_responder & (time_to_peak_ms >= sync_lo) & (time_to_peak_ms <= sync_hi);
    qc_category = zeros(1, numROIs);
    qc_category(analyzed)          = 3;   % subthreshold (default)
    qc_category(total_events >= 1) = 2;   % responder, largest peak asynchronous
    qc_category(peak_is_sync)      = 1;   % responder, largest peak synchronous
    qc_category(noise_excluded)    = 4;   % too noisy to score (analyzed is false here)
    qc_category_names = {'excluded', 'responder_sync', 'responder_async', ...
        'subthreshold', 'noise_excluded'};

    % NOTE: the group-averaged double-exponential decay (tau_fast / tau_slow /
    % slow fraction) is computed once per group by fit_group_double in
    % oneap_batch_analysis, on the pooled responder average. It is intentionally
    % NOT computed per file here.

    % --- Assemble result ---------------------------------------------------
    res = struct();
    res.dataType        = '1AP';
    res.numROIs         = numROIs;
    res.numFrames       = numFrames;
    res.stim_frame      = stim_frame;
    res.frame_rate      = config.frame_rate;
    res.frame_dur_ms    = frame_dur;
    res.prestim_idx     = prestim_idx;

    res.peak_amplitude  = peak_amplitude;
    res.peak_frame      = peak_frame;
    res.time_to_peak_ms = time_to_peak_ms;
    res.is_responder    = is_responder;
    res.prestim_sd      = prestim_sd;
    res.upper_threshold = upper_threshold;
    res.lower_threshold = lower_threshold;
    res.auc             = auc;
    res.auc_sync        = auc_sync;   % AUC over synchronous window
    res.auc_async       = auc_async;  % AUC over asynchronous window (slow tail)
    res.tau_ms          = tau_ms;
    res.tau_r2          = tau_r2;
    res.tau_A           = tau_A;
    res.tau_offset      = tau_offset;
    res.tau_peak_frame  = tau_peak_frame;
    res.prestim_base    = base_mean_roi;     % per-ROI pre-stim mean (dF/F units)

    res.sync_count      = sync_count;
    res.async_count     = async_count;
    res.total_events    = total_events;

    res.temporal_hist   = temporal_hist;        % [nBins x ROIs]
    res.bin_centers_ms  = bin_centers_ms;
    res.bin_edges_ms    = bin_edges_ms;
    res.event_peak_times_ms = event_peak_times_ms;

    res.n_responders    = sum(is_responder);          % any sustained evoked peak
    res.responder_fraction = res.n_responders / max(numROIs, 1);
    % Fraction over only the ROIs we could actually score (excludes all-NaN,
    % too-few-pre-stim, bad-SD and noise-ceiling ROIs) - the cleaner denominator.
    n_analyzable = sum(qc_category >= 1);   % codes 1/2/3 = responder/sub-threshold
    res.responder_fraction_analyzable = res.n_responders / max(n_analyzable, 1);

    % Release-timing split of the responders (by their LARGEST peak) plus the
    % raw per-window peak counts, for the synchronous/asynchronous sub-analysis.
    res.n_sync_responders   = sum(qc_category == 1);  % largest peak in sync window
    res.n_async_responders  = sum(qc_category == 2);  % largest peak in async window
    res.total_sync_events   = sum(sync_count);        % all peaks in sync window
    res.total_async_events  = sum(async_count);       % all peaks in async window

    % Gate-by-gate QC diagnostic
    res.crossed_upper      = crossed_upper;
    res.noise_excluded     = noise_excluded;
    res.qc_category        = qc_category;
    res.qc_category_names  = qc_category_names;
    res.qc_counts = struct( ...
        'excluded',         sum(qc_category == 0), ...
        'responder_sync',   sum(qc_category == 1), ...
        'responder_async',  sum(qc_category == 2), ...
        'subthreshold',     sum(qc_category == 3), ...
        'noise_excluded',   sum(qc_category == 4));
end

function a = auc_window(trace, f_win, base_mean, stim_frame, frame_dur)
    % AUC_WINDOW - Trapezoidal area under the baseline-subtracted dF/F over a
    % frame window f_win = [f1 f2]. Returns NaN if the window is too short.
    % Units: (dF/F)*ms. NaNs in the segment are pinned to 0 (no contribution).
    a  = NaN;
    f1 = max(f_win(1), 1);
    f2 = min(f_win(2), numel(trace));
    idx = f1:f2;
    if numel(idx) >= 2
        seg = trace(idx) - base_mean;
        seg(~isfinite(seg)) = 0;
        t = (idx - stim_frame) * frame_dur;     % ms from stim
        a = trapz(t, seg);
    end
end

function [tau, r2, A, C] = fit_decay_tau(t, y)
    % FIT_DECAY_TAU - Mono-exponential decay fit y ~ A*exp(-t/tau) + C
    % Uses the Curve Fitting Toolbox (nonlinear least squares) with A >= 0 and
    % tau > 0. Returns tau (ms), fit R^2, and the amplitude A and offset C so
    % the fitted curve can be reconstructed for plotting.

    tau = NaN; r2 = NaN; A = NaN; C = NaN;

    good = isfinite(t) & isfinite(y);
    t = t(good); y = y(good);
    if numel(y) < 4
        return;
    end
    t = double(t(:)); y = double(y(:));   % CFT fit() expects double

    % Initial guesses from the data
    A0   = max(y) - min(y);
    if ~(A0 > 0), A0 = max(abs(y)); end
    C0   = min(y);
    tau0 = max((t(end) - t(1)) / 2, eps);   % ~half the decay span

    ft = fittype('A*exp(-x/tau) + C', 'independent', 'x', ...
        'coefficients', {'A', 'tau', 'C'});
    opts = fitoptions(ft);
    opts.StartPoint = [A0, tau0, C0];
    opts.Lower      = [0,  1e-3, -Inf];     % A >= 0, tau > 0
    opts.Upper      = [Inf, Inf,  Inf];

    try
        [fobj, gof] = fit(t, y, ft, opts);
    catch
        return;                              % fit failed to converge
    end

    if isfinite(fobj.tau) && fobj.tau > 0
        tau = fobj.tau;
        r2  = gof.rsquare;
        A   = fobj.A;
        C   = fobj.C;
    end
end
