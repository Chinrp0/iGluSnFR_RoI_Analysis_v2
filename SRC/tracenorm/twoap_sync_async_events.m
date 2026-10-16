function ev = twoap_sync_async_events(trace, lo, hi, base, sd, config, allow_end)
    % TWOAP_SYNC_ASYNC_EVENTS - Detect every release event in one window, binned sync / async
    %
    % An event is a local maximum of the dF/F trace inside frames [lo, hi] that
    %   - clears the UPPER threshold (base + upper_threshold_sigma * sd),
    %   - passes the same Schmitt run as the responder gates (>= min_frames_above_lower
    %     consecutive frames >= base + lower_threshold_sigma * sd), and
    %   - RISES at least events.min_prominence_sigma * sd above the trough before it
    %     (the lowest frame since the previous event, or since the frame before lo),
    %     so noise riding on a decaying response is not counted as a new event.
    %     Only the trough BEFORE the peak counts: a window cut short by the next
    %     stim must not make a real peak fail for lack of a decay after it.
    % Events closer than min_peak_distance_ms keep the larger one.
    %
    % Binning is on event ONSET, not peak: the onset is the first frame of the
    % final climb above onset_frac of the way from the trough (lowest frame in
    % the onset_trough_ms before the peak) to the peak. Events whose onset is <= events.sync_cutoff_ms
    % after lo (rounded to frames) are SYNC, later ones ASYNC; the peak can land
    % later (3fast takes ~12 ms to rise, peaking 15-20 ms after release). lo is the
    % stim frame for a pulse window; for a spontaneous (pre-stim) window the bins
    % are ignored.
    %
    % The frames on either side of [lo, hi] decide whether its first / last frame
    % is a local maximum; only peaks inside [lo, hi] count. With allow_end (pulse-1
    % window cut short by stim2) the last frame is a candidate whenever it is above
    % the frame before it (pulse 1 still rising when the window ends), the same
    % rule as the pulse-1 responder gate.
    %
    % Inputs:
    %   trace     - [nF x 1] dF/F trace
    %   lo, hi    - window frames (lo = stim frame for a pulse window)
    %   base      - pre-stim1 baseline level of this trace
    %   sd        - pre-stim noise SD of this trace
    %   config    - tracenorm_config() struct (config.twoAP, config.twoAP.events,
    %               config.frame_dur_ms)
    %   allow_end - optional, default false
    %
    % Output:
    %   ev - struct:
    %        lat_ms    [1 x K] event peak latency after lo (ms)
    %        onset_ms  [1 x K] event onset latency after lo (ms, >= 0)
    %        height    [1 x K] event peak above base (dF/F)
    %        prom      [1 x K] rise above the trough before it (dF/F)
    %        is_sync   [1 x K] logical, onset <= sync cutoff
    %        n_sync, n_async   event counts per bin
    %        sync_ms, async_ms window length (ms) covered by each bin
    if nargin < 7, allow_end = false; end
    p  = config.twoAP;
    e  = p.events;
    fd = config.frame_dur_ms;
    trace = double(trace(:));
    nF = numel(trace);
    lo = max(1, lo); hi = min(nF, hi);

    ks = max(0, round(e.sync_cutoff_ms / fd));      % last sync onset frame after lo
    nWin = max(0, hi - lo + 1);
    ev = struct('lat_ms', zeros(1, 0), 'onset_ms', zeros(1, 0), 'height', zeros(1, 0), 'prom', zeros(1, 0), ...
        'is_sync', false(1, 0), 'n_sync', 0, 'n_async', 0, ...
        'sync_ms', min(ks + 1, nWin) * fd, 'async_ms', max(nWin - ks - 1, 0) * fd);
    if nWin < 1 || ~(isfinite(sd) && sd > 0) || ~isfinite(base), return; end

    % Window plus one frame either side (or -Inf past the trace / with allow_end).
    seg = [-Inf; trace(lo:hi); -Inf];
    if lo > 1, seg(1) = trace(lo - 1); end
    if hi < nF && ~allow_end, seg(end) = trace(hi + 1); end
    fin = isfinite(seg(2:end-1));
    if ~any(fin), return; end
    inner = seg(2:end-1); inner(~fin) = min(inner(fin)); seg(2:end-1) = inner;   % NaNs can't be peaks

    u    = base + p.upper_threshold_sigma * sd;
    l    = base + p.lower_threshold_sigma * sd;
    prom = e.min_prominence_sigma * sd;
    mpd  = max(1, round(p.min_peak_distance_ms / fd));

    c    = 2:numel(seg) - 1;                         % window frames, in seg indices
    cand = c(seg(c) > seg(c - 1) & seg(c) >= seg(c + 1) & seg(c) >= u);
    % Walk the candidates in time: rise above the lowest frame since the last
    % accepted event; within mpd of it, the higher of the two is kept.
    acc = []; rise = []; from = 1;
    for k = cand
        tr0 = min(seg(max(from, 1):k));
        if ~isfinite(tr0), tr0 = min(seg(max(from, 2):k)); end
        r = seg(k) - tr0;
        if ~isempty(acc) && k - acc(end) < mpd
            if seg(k) > seg(acc(end)), acc(end) = k; rise(end) = max(rise(end), r); from = k; end
            continue;
        end
        if r >= prom
            acc(end+1) = k; rise(end+1) = r; from = k; %#ok<AGROW>
        end
    end
    absloc = lo + acc - 2;
    % Onset: trough = lowest frame in the onset_trough_ms before the peak (not
    % reaching back past the previous event); onset = first frame of the final
    % climb above trough + onset_frac * (peak - trough). A one-frame noise dip
    % part-way up the rise only moves the onset if it falls back near the trough.
    % (seg index j -> absolute frame lo + j - 2), never before lo.
    nTr   = max(1, round(e.onset_trough_ms / fd));
    onset = zeros(size(acc));
    prev  = 1;
    for k = 1:numel(acc)
        a  = acc(k);
        w  = max([prev, a - nTr, 1]):a;
        sw = seg(w); sw(~isfinite(sw)) = Inf;
        [tv, it] = min(sw);
        t0 = w(it);
        lvl = tv + e.onset_frac * (seg(a) - tv);
        j = find(seg(t0:a) <= lvl, 1, 'last') + t0 - 1;   % last frame at/below the level
        onset(k) = max(lo + j - 1, lo);                    % the frame after it
        prev = a;
    end

    % Schmitt run around each peak (absolute frames), as in the responder gates.
    run_ok = false(size(absloc));
    for k = 1:numel(absloc)
        li = absloc(k); ri = absloc(k);
        while li > 1  && isfinite(trace(li-1)) && trace(li-1) >= l, li = li - 1; end
        while ri < nF && isfinite(trace(ri+1)) && trace(ri+1) >= l, ri = ri + 1; end
        run_ok(k) = (ri - li + 1) >= p.min_frames_above_lower;
    end
    absloc = absloc(run_ok); rise = rise(run_ok); onset = onset(run_ok);

    ev.lat_ms   = (absloc - lo) * fd;
    ev.onset_ms = (onset - lo) * fd;
    ev.height   = trace(absloc)' - base;
    ev.prom     = rise;
    ev.is_sync  = (onset - lo) <= ks;
    ev.n_sync  = sum(ev.is_sync);
    ev.n_async = sum(~ev.is_sync);
end
