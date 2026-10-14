function m = hfs_train_metrics(trace, base, config, do_decay)
    % HFS_TRAIN_METRICS - Per-pulse and train-level metrics of ONE dF/F trace
    %
    % Shared by hfs_analyzer (each ROI) and hfs_batch_analysis (coverslip
    % grand-average traces), so both use the same windows and definitions. The
    % train timing comes from config.hfs (set per file / per set by hfs_pipeline
    % and hfs_batch_analysis):
    %   stim1_frame     - first pulse (epoch index)
    %   pulse_frames    - [1 x nPulses] frame of every pulse
    %   isi_ms          - inter-pulse interval (1000 / Train_Hz)
    %   train_end_frame - one inter-pulse interval after the last pulse
    %
    % Per pulse n (only when the pulses are resolvable, i.e. at least
    % config.hfs.min_frames_per_pulse frames apart):
    %   tonic(n)    - mean of the last residual_window_ms before pulse n, above base
    %   amp(n)      - phasic amplitude: peak in [pulse n, pulse n + peak_search_ms]
    %                 (capped one frame before pulse n+1) minus that pre-pulse level
    %   lat_ms(n)   - latency of that peak after pulse n
    % Train-level (any frame rate):
    %   q_train     - integral of (trace - base) over [stim1, train end + auc_post_ms],
    %                 dF/F*ms (total glutamate charge of the train)
    %   plateau     - mean (trace - base) over the last plateau_window_ms of the train
    %   peak_train  - max (trace - base) during the train
    % Derived from the per-pulse values:
    %   a1, norm_amp (amp / a1), ss_amp (mean amp over the last steady_state_pulses),
    %   ss_ratio (ss_amp / a1), ss_tonic, cum_amp (cumsum of amp), and the
    %   back-extrapolated line through the last cumulative_fit_pulses of cum_amp:
    %   pool_est (intercept at pulse 0) and repl_per_pulse (slope).
    % Post-train decay (do_decay = true; meant for averaged traces):
    %   decay_tau_ms (single exponential + offset, kept when R^2 >= decay_min_r2),
    %   decay_r2, t_half_ms (model-free: time after train end until the trace first
    %   falls to half its end-of-train level).
    %
    % Inputs:
    %   trace    - [frames x 1] dF/F
    %   base     - pre-train baseline level (dF/F); [] = mean of frames before stim1
    %   config   - tracenorm_config() struct with config.hfs timing set
    %   do_decay - fit the post-train decay (default false)

    if nargin < 4 || isempty(do_decay), do_decay = false; end
    p    = config.hfs;
    fd   = config.frame_dur_ms;
    trace = double(trace(:));
    nF   = numel(trace);
    s1   = p.stim1_frame;
    pf   = p.pulse_frames(:)';
    nP   = numel(pf);
    tend = min(p.train_end_frame, nF);

    if isempty(base)
        pre = trace(1:max(s1 - 1, 1)); base = mean(pre(isfinite(pre)));
    end

    m = struct();
    m.base       = base;
    m.resolvable = hfs_pulses_resolvable(config);

    % --- Train-level --------------------------------------------------------
    post = max(1, round(p.auc_post_ms / fd));
    m.q_train = trapz_ms(trace, s1:min(tend + post, nF), base, fd);
    pw = max(1, round(p.plateau_window_ms / fd));
    seg = trace(max(tend - pw, s1):tend - 1) - base;
    m.plateau = mean(seg(isfinite(seg)));
    seg = trace(s1:tend - 1) - base;
    m.peak_train = max(seg(isfinite(seg)));
    if isempty(m.peak_train), m.peak_train = NaN; end

    % --- Per-pulse ----------------------------------------------------------
    m.amp = nan(1, nP); m.tonic = nan(1, nP); m.lat_ms = nan(1, nP);
    if m.resolvable
        nres   = max(1, round(p.residual_window_ms / fd));
        search = max(1, round(p.peak_search_ms / fd));
        for n = 1:nP
            s = pf(n);
            if n < nP, nxt = pf(n + 1); else, nxt = tend; end
            w = s:min([s + search, nxt - 1, nF]);
            r = max(s - nres, 1):s - 1;
            if isempty(w) || isempty(r), continue; end
            pre = trace(r); pre = pre(isfinite(pre));
            seg = trace(w);
            if isempty(pre) || ~any(isfinite(seg)), continue; end
            lvl = mean(pre);
            [pk, k] = max(seg);
            m.tonic(n)  = lvl - base;
            m.amp(n)    = pk - lvl;
            m.lat_ms(n) = (w(k) - s) * fd;
        end
    end
    m.a1       = m.amp(1);
    m.norm_amp = m.amp / m.a1;
    if ~(m.a1 > 0), m.norm_amp(:) = NaN; end
    nss = min(p.steady_state_pulses, nP);
    m.ss_amp   = mean(m.amp(end - nss + 1:end), 'omitnan');
    m.ss_tonic = mean(m.tonic(end - nss + 1:end), 'omitnan');
    m.ss_ratio = m.ss_amp / m.a1;
    if ~(m.a1 > 0), m.ss_ratio = NaN; end
    m.cum_amp  = cumsum(m.amp);
    m.pool_est = NaN; m.repl_per_pulse = NaN;
    nfit = min(p.cumulative_fit_pulses, nP);
    k = (nP - nfit + 1):nP;
    if nfit >= 3 && all(isfinite(m.cum_amp(k)))
        c = polyfit(k, m.cum_amp(k), 1);
        m.repl_per_pulse = c(1);
        m.pool_est       = c(2);
    end
    m.median_lat_ms = median(m.lat_ms, 'omitnan');

    % --- Post-train decay ---------------------------------------------------
    m.decay_tau_ms = NaN; m.decay_r2 = NaN; m.t_half_ms = NaN;
    if do_decay
        [m.decay_tau_ms, m.decay_r2, m.t_half_ms] = post_train_decay(trace, base, pf(end), tend, config);
    end
end

function [tau, r2, t_half] = post_train_decay(trace, base, last_pulse, tend, config)
    % Decay after the last pulse: starts at the largest point between the last
    % pulse and the train end, and runs for decay_window_ms (clipped to the trace).
    p = config.hfs; fd = config.frame_dur_ms; nF = numel(trace);
    tau = NaN; r2 = NaN; t_half = NaN;
    seg = trace(last_pulse:min(tend, nF));
    if ~any(isfinite(seg)), return; end
    [~, k] = max(seg);
    d0 = last_pulse + k - 1;
    d1 = min(d0 + round(p.decay_window_ms / fd), nF);
    y = trace(d0:d1) - base; t = (0:numel(y) - 1)' * fd;
    ok = isfinite(y); y = y(ok); t = t(ok);
    if numel(y) < 6 || ~(y(1) > 0), return; end

    % Model-free half-decay: first time a 3-frame running mean drops to y(1)/2.
    ys = movmean(y, 3);
    h  = find(ys <= y(1) / 2, 1);
    if ~isempty(h), t_half = t(h); end

    % Single exponential + offset (Curve Fitting Toolbox, as fit_decay_double).
    try
        ft = fittype('A*exp(-x/tau) + C', 'independent', 'x', 'coefficients', {'A', 'tau', 'C'});
        opts = fitoptions(ft);
        opts.StartPoint = [y(1), max(t(end) / 5, fd), 0];
        opts.Lower      = [0, fd / 2, -Inf];
        opts.Upper      = [Inf, 10 * t(end), Inf];
        [f, gof] = fit(t, y, ft, opts);
        r2 = gof.rsquare;
        if r2 >= p.decay_min_r2, tau = f.tau; end
    catch
    end
end

function a = trapz_ms(trace, win, base, fd)
    a = NaN;
    if numel(win) < 2, return; end
    seg = trace(win) - base; seg = seg(isfinite(seg));
    if numel(seg) >= 2, a = trapz(seg) * fd; end
end
