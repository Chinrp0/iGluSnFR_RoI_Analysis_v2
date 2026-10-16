function amp = twoap_pulse_amplitudes(trace, stim1, stim2, base, config, template)
    % TWOAP_PULSE_AMPLITUDES - Pulse-1 / pulse-2 amplitudes of one dF/F trace, five ways
    %
    % Shared by twoap_analyzer (one ROI) and twoap_batch_analysis (averaged
    % traces), so a per-ROI amplitude and an average-trace amplitude are always
    % measured with the same frames, baselines and peak rule. A1 is above the
    % pre-stim1 baseline (base); A2 is INCREMENTAL, above the local pre-stim2
    % residual (mean of the twoap_residual_window frames), except decay_sub.
    %
    %   peak_max  - largest local maximum in [stim, stim + peak_search_ms] (pulse-1
    %               window capped one frame before stim2). Window max if the
    %               window has no interior local maximum.
    %   win_mean  - mean dF/F over the fixed window amp_window_ms after each stim
    %               (pulse-1 window capped before stim2). No max-of-noise inflation.
    %   sync      - single frame at sync_latency_ms after each stim.
    %   p1lock    - A1 = peak_max A1; A2 read at pulse-1's own peak latency after stim2.
    %   decay_sub - A1 = peak_max A1; A2 = peak_max of (trace - base - pulse-1 tail),
    %               where the tail is the pulse-1 template scaled to the trace's
    %               residual level. Corrects for pulse 1 still decaying under
    %               pulse 2. NaN without a template.
    %
    % Inputs:
    %   trace    - [nF x 1] dF/F trace
    %   stim1/2  - stimulus frames
    %   base     - pre-stim1 baseline level of this trace
    %   config   - tracenorm_config() struct (config.twoAP, config.frame_dur_ms)
    %   template - pulse-1 waveform, baseline-subtracted and peak-normalised,
    %              template(k) = frame stim1 + k - 1 ([] -> decay_sub is NaN)
    %
    % Output:
    %   amp - struct: amp.<method>.a1 / .a2 for the five methods, plus
    %         resid (residual level), a2_abs (pulse-2 peak above base),
    %         t1 / t2 (peak_max peak latencies in frames after each stim),
    %         t2_decay (decay_sub pulse-2 peak latency, frames after stim2),
    %         decay_scale (template scale s; NaN when the flat residual was used)
    if nargin < 6, template = []; end
    p  = config.twoAP;
    fd = config.frame_dur_ms;
    trace = double(trace(:));
    nF = numel(trace);

    search = max(1, round(p.peak_search_ms / fd));
    isi_f  = stim2 - stim1;
    s1_end = min(search, max(isi_f - 1, 1));     % pulse-1 window stops before stim2
    w1 = stim1:min(stim1 + s1_end, nF);
    w2 = stim2:min(stim2 + search, nF);

    [r_lo, r_hi] = twoap_residual_window(stim1, stim2, p.residual_window_ms, fd);
    rseg  = trace(r_lo:r_hi);
    resid = mean(rseg(isfinite(rseg)));
    if ~isfinite(resid), resid = base; end

    nanpair = struct('a1', NaN, 'a2', NaN);
    amp = struct('peak_max', nanpair, 'win_mean', nanpair, 'sync', nanpair, ...
        'p1lock', nanpair, 'decay_sub', nanpair, 'resid', resid, 'a2_abs', NaN, ...
        't1', NaN, 't2', NaN, 't2_decay', NaN, 'decay_scale', NaN);

    % --- peak_max ---
    % When stim2 cuts the pulse-1 window short (short ISI), pulse 1 can still be
    % rising on the window's last frame; that frame then counts as a peak if the
    % trace rose into it (see local_peak).
    [p1v, i1] = local_peak(trace(w1), s1_end < search);
    [p2v, i2] = local_peak(trace(w2));
    amp.peak_max.a1 = p1v - base;
    amp.peak_max.a2 = p2v - resid;
    amp.a2_abs      = p2v - base;
    amp.t1 = i1 - 1; amp.t2 = i2 - 1;

    % --- win_mean ---
    aw = getfielddef(p, 'amp_window_ms', [5 20]);
    k0 = max(0, round(aw(1) / fd)); k1 = max(k0, round(aw(2) / fd));
    k1a = min(k1, isi_f - 1);
    if k1a >= k0 && stim1 + k1a <= nF
        amp.win_mean.a1 = mean_finite(trace(stim1 + (k0:k1a))) - base;
    end
    if stim2 + k1 <= nF
        amp.win_mean.a2 = mean_finite(trace(stim2 + (k0:k1))) - resid;
    end

    % --- sync ---
    ks = max(1, round(p.sync_latency_ms / fd));
    if ks < isi_f && stim1 + ks <= nF, amp.sync.a1 = trace(stim1 + ks) - base; end
    if stim2 + ks <= nF, amp.sync.a2 = trace(stim2 + ks) - resid; end

    % --- p1lock ---
    amp.p1lock.a1 = amp.peak_max.a1;
    if isfinite(amp.t1) && stim2 + amp.t1 <= nF
        amp.p1lock.a2 = trace(stim2 + amp.t1) - resid;
    end

    % --- decay_sub ---
    amp.decay_sub.a1 = amp.peak_max.a1;
    if ~isempty(template)
        T = double(template(:));
        tr_idx = (r_lo:r_hi) - stim1 + 1;
        w2_idx = w2 - stim1 + 1;
        t_r = NaN;
        if tr_idx(end) <= numel(T), t_r = mean(T(tr_idx)); end
        if t_r > getfielddef(p, 'decay_template_min_frac', 0.05) ...
                && w2_idx(end) <= numel(T) && all(isfinite(T(w2_idx)))
            % Scale the template so its value at the residual window equals this
            % trace's residual above base, then subtract the projected pulse-1
            % tail from the pulse-2 window.
            s   = (resid - base) / t_r;
            sub = trace(w2) - base - s * T(w2_idx);
            amp.decay_scale = s;
        else
            % Pulse 1 has decayed (long ISI, or past the template's end): flat
            % residual, same as peak_max.
            sub = trace(w2) - resid;
        end
        [amp.decay_sub.a2, i2d] = local_peak(sub);
        amp.t2_decay = i2d - 1;
    end
end

function [v, i] = local_peak(seg, allow_end)
    % Largest interior local maximum of seg (same rule as findpeaks: endpoints
    % never count); window max when there is none. With allow_end, the last
    % sample is also a candidate when it is above the one before it (a window
    % cut short while the response is still rising).
    if nargin < 2, allow_end = false; end
    v = NaN; i = NaN;
    fin = isfinite(seg);
    if ~any(fin), return; end
    seg(~fin) = min(seg(fin));
    n = numel(seg);
    cand = [];
    if n >= 3
        c = 2:n-1;
        cand = c(seg(c) > seg(c-1) & seg(c) >= seg(c+1));
    end
    if allow_end && n >= 2 && seg(n) > seg(n-1), cand(end+1) = n; end
    if ~isempty(cand)
        [v, k] = max(seg(cand)); i = cand(k);
        return;
    end
    [v, i] = max(seg);
end

function m = mean_finite(x)
    x = x(isfinite(x));
    if isempty(x), m = NaN; else, m = mean(x); end
end

function v = getfielddef(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
