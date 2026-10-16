function test_twoap_pulse_amplitudes()
    % TEST_TWOAP_PULSE_AMPLITUDES - Synthetic checks of the five 2AP A1/A2 definitions
    %
    % Builds noise-free paired-pulse traces from a known pulse waveform
    % (fast rise, single-exponential decay) at 5 ms and 10 ms frames and checks:
    %   - peak_max A1 = the sampled pulse-1 peak; sync / win_mean match the
    %     waveform sampled at their frames
    %   - decay_sub A2 recovers the true pulse-2 amplitude at every ISI (pulse 1's
    %     tail is subtracted with the exact template)
    %   - peak_max A2 (flat residual) UNDER-estimates A2 at short ISI, where pulse 1
    %     is still decaying under pulse 2, and converges at long ISI
    %   - twoap_analyzer's per-ROI amplitudes equal the shared function's
    % Run: test_twoap_pulse_amplitudes  (errors on failure)

    evalc('cfg = tracenorm_config();');
    for fd = [5 10]
        cfg.frame_dur_ms = fd;
        nF = 700; stim1 = 201; base = 0.002;
        A1 = 0.010; A2 = 0.015;                         % true amplitudes (PPR 1.5)
        tms = (0:nF - 1)' * fd;
        kern = @(t) (t > 0) .* (1 - exp(-max(t, 0) / 3)) .* exp(-max(t, 0) / 40);
        % Template = the sampled pulse-1 waveform, peak-normalised
        k1 = kern(tms(stim1:end) - tms(stim1));
        kpk = max(k1(1:max(1, round(cfg.twoAP.peak_search_ms / fd)) + 1));
        T = k1 / kpk;

        for isi = [30 60 90 250 1000]
            isi_f = round(isi / fd); stim2 = stim1 + isi_f;
            if stim2 + 50 > nF, continue; end
            tr = base + A1 * kern(tms - tms(stim1)) / kpk + A2 * kern(tms - tms(stim2)) / kpk;
            amp = twoap_pulse_amplitudes(tr, stim1, stim2, base, cfg, T);

            assert(abs(amp.peak_max.a1 - A1) < 1e-9, 'peak_max A1 (fd %g, ISI %g)', fd, isi);
            ks = max(1, round(cfg.twoAP.sync_latency_ms / fd));
            if ks < isi_f
                assert(abs(amp.sync.a1 - (tr(stim1 + ks) - base)) < 1e-12, 'sync A1');
            end
            aw = cfg.twoAP.amp_window_ms;
            k0 = max(0, round(aw(1) / fd)); k1w = min(round(aw(2) / fd), isi_f - 1);
            if k1w >= k0
                assert(abs(amp.win_mean.a1 - (mean(tr(stim1 + (k0:k1w))) - base)) < 1e-12, 'win_mean A1');
            end
            % decay_sub recovers A2 (to within the residual-window sampling)
            assert(abs(amp.decay_sub.a2 - A2) / A2 < 0.02, ...
                'decay_sub A2 %.5f vs true %.5f (fd %g, ISI %g)', amp.decay_sub.a2, A2, fd, isi);
            if isi <= 60
                assert(amp.peak_max.a2 < 0.97 * A2, ...
                    'peak_max A2 should underestimate at ISI %g (got %.5f)', isi, amp.peak_max.a2);
            end
            if isi >= 250
                assert(abs(amp.peak_max.a2 - A2) / A2 < 0.02, 'peak_max A2 at long ISI');
            end
            fprintf('  fd %2g ms ISI %4g: PPR peak_max %.3f | decay_sub %.3f | win_mean %.3f | sync %.3f | p1lock %.3f (true 1.500)\n', ...
                fd, isi, amp.peak_max.a2 / amp.peak_max.a1, amp.decay_sub.a2 / amp.decay_sub.a1, ...
                amp.win_mean.a2 / amp.win_mean.a1, amp.sync.a2 / amp.sync.a1, amp.p1lock.a2 / amp.p1lock.a1);
        end
    end

    % Analyzer amplitudes == shared function (noisy multi-ROI trace)
    cfg.frame_dur_ms = 5; cfg.frame_rate = 200;
    rng(1);
    nR = 20; stim1 = 201; stim2 = stim1 + 12;
    cfg.twoAP.stim1_frame = stim1; cfg.twoAP.stim2_frame = stim2; cfg.twoAP.isi_ms = 60;
    tms = (0:nF - 1)' * 5;
    kern = @(t) (t > 0) .* (1 - exp(-max(t, 0) / 3)) .* exp(-max(t, 0) / 40);
    D = 0.002 * randn(nF, nR) + kern(tms - tms(stim1)) * (0.03 * rand(1, nR)) ...
        + kern(tms - tms(stim2)) * (0.03 * rand(1, nR));
    res = twoap_analyzer(D, struct(), cfg);
    for r = 1:nR
        if ~isfinite(res.prestim_base(r)), continue; end
        amp = twoap_pulse_amplitudes(D(:, r), stim1, stim2, res.prestim_base(r), cfg, []);
        assert(abs(res.peak1_amp(r) - amp.peak_max.a1) < 1e-12, 'analyzer peak1_amp ROI %d', r);
        if amp.peak_max.a2 >= 0
            assert(abs(res.peak2_amp(r) - amp.peak_max.a2) < 1e-12, 'analyzer peak2_amp ROI %d', r);
        end
        assert(isequaln(res.a1_sync(r), amp.sync.a1) && isequaln(res.a2_sync(r), amp.sync.a2), 'analyzer sync ROI %d', r);
    end

    % Short-ISI edge case: 30 ms ISI at 10 ms frames leaves a pulse-1 window of
    % frames 0/10/20 ms; a slower-rising pulse peaks at 20 ms = the window's last
    % frame. It must still count as the pulse-1 peak (gate + amplitude).
    cfg.frame_dur_ms = 10; cfg.frame_rate = 100;
    stim1 = 101; stim2 = stim1 + 3; nF = 300; nR = 10;
    cfg.twoAP.stim1_frame = stim1; cfg.twoAP.stim2_frame = stim2; cfg.twoAP.isi_ms = 30;
    tms = (0:nF - 1)' * 10;
    kslow = @(t) (t > 0) .* (1 - exp(-max(t, 0) / 12)) .* exp(-max(t, 0) / 40);   % peak ~18 ms
    rng(2);
    D = 0.0002 * randn(nF, nR) + kslow(tms - tms(stim1)) * (0.05 * ones(1, nR)) ...
        + kslow(tms - tms(stim2)) * (0.05 * ones(1, nR));
    res = twoap_analyzer(D, struct(), cfg);
    edgeVal = D(stim1 + 2, :) - res.prestim_base;
    assert(all(res.is_responder1), 'edge pulse-1 peak should pass the gate (%d/%d)', sum(res.is_responder1), nR);
    assert(all(abs(res.peak1_amp - edgeVal) < 1e-12), 'A1 should be the last-frame (20 ms) value');
    fprintf('  edge case (30 ms ISI, 10 ms frames): %d/%d pulse-1 responders, A1 = 20 ms frame\n', ...
        sum(res.is_responder1), nR);
    fprintf('test_twoap_pulse_amplitudes: all checks passed.\n');
end
