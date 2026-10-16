function test_twoap_sync_async_events()
    % TEST_TWOAP_SYNC_ASYNC_EVENTS - Synthetic checks of the 2AP sync / async event binning
    %
    % 5 ms frames, iGluSnFR3fast-like events (peak ~7 ms after release):
    %   - noise-free: a sync event after each stim and one async event released at
    %     50 ms are found and binned (onset <= sync_cutoff_ms = sync, later = async)
    %   - the pulse-1 window stops before stim2 (a short-ISI pulse-2 response is
    %     not counted as pulse-1 async); a peak at the cut window's end is found
    %   - binning uses ONSET: release at +12 ms is async though it peaks ~20 ms
    %   - noise: the false async rate on a noisy decaying sync response, and the
    %     false event rate in a noise-only baseline, stay small
    % Run: test_twoap_sync_async_events  (errors on failure)

    evalc('cfg = tracenorm_config();');
    cfg.frame_dur_ms = 5;
    fd = 5; nF = 600; stim1 = 201; base = 0; sd = 0.002;
    tms  = (0:nF - 1)' * fd;
    kern = @(t) (t > 0) .* (1 - exp(-max(t, 0) / 6)) .* exp(-max(t, 0) / 15);
    kpk  = max(kern(0:0.1:200));
    ev_at = @(t0, a) a * kern(tms - t0) / kpk;     % event released at t0 (ms), peak a
    aw_f = round(cfg.twoAP.events.async_window_ms / fd);
    t1 = tms(stim1);

    % --- noise-free: sync + async after pulse 1, sync after pulse 2 ---
    stim2 = stim1 + round(250 / fd);
    tr = base + ev_at(t1, 0.02) + ev_at(t1 + 50, 0.015) + ev_at(tms(stim2), 0.02);
    e1 = twoap_sync_async_events(tr, stim1, min(stim2 - 1, stim1 + aw_f), base, sd, cfg);
    e2 = twoap_sync_async_events(tr, stim2, stim2 + aw_f, base, sd, cfg);
    assert(e1.n_sync == 1 && e1.n_async == 1, 'pulse 1: expected 1 sync + 1 async, got %d + %d', ...
        e1.n_sync, e1.n_async);
    assert(e1.lat_ms(1) <= cfg.twoAP.events.sync_cutoff_ms && e1.lat_ms(2) > cfg.twoAP.events.sync_cutoff_ms, ...
        'pulse 1 latencies %s', mat2str(e1.lat_ms));
    assert(e2.n_sync == 1 && e2.n_async == 0, 'pulse 2: expected 1 sync, got %d + %d', e2.n_sync, e2.n_async);
    assert(abs(e1.async_ms - (aw_f - round(cfg.twoAP.events.sync_cutoff_ms / fd)) * fd) < 1e-9, 'async bin length');

    % --- short ISI: pulse-1 window ends before stim2 ---
    stim2 = stim1 + round(60 / fd);
    tr = base + ev_at(t1, 0.02) + ev_at(tms(stim2), 0.02);
    e1 = twoap_sync_async_events(tr, stim1, min(stim2 - 1, stim1 + aw_f), base, sd, cfg);
    assert(e1.n_sync == 1 && e1.n_async == 0, 'short ISI pulse 1: %d + %d', e1.n_sync, e1.n_async);
    e2 = twoap_sync_async_events(tr, stim2, stim2 + aw_f, base, sd, cfg);
    assert(e2.n_sync == 1 && e2.n_async == 0, 'short ISI pulse 2 (on pulse-1 decay): %d + %d', ...
        e2.n_sync, e2.n_async);

    % --- 30 ms ISI: pulse-1 window (0-25 ms) cut short by stim2 ---
    % A peak near the window end still counts (its rise is measured from the trough
    % BEFORE it, not limited by the still-elevated stim2 frame after it) ...
    stim2 = stim1 + round(30 / fd);
    w1 = [stim1, stim2 - 1];
    tr = base + ev_at(t1 + 12, 0.02) + ev_at(tms(stim2), 0.02);
    e1 = twoap_sync_async_events(tr, w1(1), w1(2), base, sd, cfg, true);
    assert(e1.n_sync + e1.n_async == 1, '30 ms ISI: pulse-1 peak near the window end missed');
    % ... and a pulse 1 still rising on the last frame counts only with allow_end.
    tr = base + ev_at(t1 + 22, 0.05);
    e1 = twoap_sync_async_events(tr, w1(1), w1(2), base, sd, cfg, true);
    assert(e1.n_sync + e1.n_async == 1, '30 ms ISI: rising last frame not counted with allow_end');
    e1 = twoap_sync_async_events(tr, w1(1), w1(2), base, sd, cfg, false);
    assert(e1.n_sync + e1.n_async == 0, '30 ms ISI: rising last frame counted without allow_end');

    % --- binning is on ONSET: release at +3 ms is sync, release at +12 ms is
    % async even though its peak (~20 ms) lands where sync peaks usually do ---
    for rel = [3 12]
        tr = base + ev_at(t1 + rel, 0.03);
        e1 = twoap_sync_async_events(tr, stim1, stim1 + aw_f, base, sd, cfg);
        assert(isscalar(e1.lat_ms),'release +%g ms: expected 1 event', rel);
        assert(e1.onset_ms <= e1.lat_ms, 'onset after peak');
        assert(e1.is_sync == (rel < cfg.twoAP.events.sync_cutoff_ms), ...
            'release +%g ms: onset %g ms, peak %g ms, is_sync %d', rel, e1.onset_ms, e1.lat_ms, e1.is_sync);
    end

    % --- noise: false async on a decaying sync response, false spont events ---
    rng(1);
    nRep = 2000; fa = 0; fs = 0; ns = 0;
    for r = 1:nRep
        tr = base + ev_at(t1, 0.03) + sd * randn(nF, 1);
        e1 = twoap_sync_async_events(tr, stim1, stim1 + aw_f, base, sd, cfg);
        es = twoap_sync_async_events(tr, 1, stim1 - 1, base, sd, cfg);
        fa = fa + e1.n_async; ns = ns + (e1.n_sync > 0); fs = fs + numel(es.lat_ms);
    end
    fprintf('Noise (sd %.3f, 15-sigma sync): P(sync) %.3f, false async %.4f per ROI, false spont %.4f per s\n', ...
        sd, ns / nRep, fa / nRep, fs / nRep / ((stim1 - 1) * fd) * 1000);
    assert(ns / nRep > 0.99, 'sync event missed in noise');
    assert(fa / nRep < 0.01, 'false async rate %.4f per ROI', fa / nRep);
    assert(fs / nRep / ((stim1 - 1) * fd) * 1000 < 0.05, 'false spont rate');
    fprintf('test_twoap_sync_async_events: all checks passed\n');
end
