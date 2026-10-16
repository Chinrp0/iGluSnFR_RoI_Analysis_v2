function out = deconv_events(trace, frame_ms, tau_ms, p, tau_rise_ms)
    % DECONV_EVENTS - Deconvolution-based release-event detection and noise floor
    %
    % After Nat Commun 2022, doi:10.1038/s41467-022-31070-4 (Methods): band-pass the
    % trace, deconvolve the unitary iGluSnFR response (instantaneous rise,
    % single-exponential decay tau), band-pass again, fit one Gaussian to the
    % core of the deconvolved all-points histogram to get the baseline-noise
    % SD sigma_BN, and call every local maximum above theta = k * sigma_BN an
    % event.
    %
    % Adapted to this pipeline's recordings:
    %   - The deconvolution is the exact inverse of the SAMPLED kernel, done
    %     in the time domain. This equals the paper's Fourier division
    %     F^-1[F(x)/F(kernel)] but without its circular wrap-around at the
    %     epoch edges. d is in dF/F units: a single event of peak amplitude A
    %     gives d = A at its onset frame.
    %       tau_rise_ms = 0 (paper): instantaneous rise, kernel a^n,
    %         a = exp(-frame_ms / tau_ms)          -> AR(1) inverse
    %           d[n] = x[n] - a * x[n-1]
    %       tau_rise_ms > 0: difference of exponentials (a^n - b^n) / kpk,
    %         b = exp(-frame_ms / tau_rise_ms), kpk = its continuous peak
    %                                               -> AR(2) inverse
    %           d[n] = kpk / (a - b) * (x[n+1] - (a + b) x[n] + a b x[n-1])
    %       The rise matters when it spans more than a frame (SF-iGluSnFR4f at
    %       10 ms): with an instantaneous-rise kernel the slow rise is spread
    %       over 2-3 frames of d and its peak can fall under threshold.
    %   - The band-pass edges are in Hz and the upper edge is capped at
    %     p.max_hi_frac_nyquist of the Nyquist rate (a 30 Hz low-pass means
    %     nothing at 25 ms frames, where Nyquist is 20 Hz).
    %   - The band-pass is a Gaussian frequency response applied by FFT to a
    %     mirror-padded trace (no wrap-around). Each edge is the -3 dB point.
    %   - NaN frames (epoch padding) are trimmed from the ends; interior NaNs
    %     are linearly interpolated for filtering and NaN again in the output.
    %
    % USE FOR TIMING, DETECTION AND NOISE ONLY. Both filters change amplitudes
    % (the 0.5 Hz high-pass also removes slow tonic / plateau signal), so keep
    % measuring amplitudes and charge on the dF/F trace.
    %
    % Inputs:
    %   trace    - [nF x 1] dF/F trace
    %   frame_ms - frame interval (ms), from Recordings.csv
    %   tau_ms   - unitary decay time constant (ms) for this sensor / prep
    %   tau_rise_ms - unitary rise time constant (ms); 0 or omitted = instantaneous
    %   p        - struct (config.deconv): band_hz [lo hi], max_hi_frac_nyquist,
    %              theta_sigma, hist_core_sigma, hist_bins, min_peak_distance_ms
    %
    % Output struct:
    %   filtered   [nF x 1] band-passed dF/F
    %   deconv     [nF x 1] filtered deconvolved trace (dF/F jump units)
    %   sigma_bn   baseline-noise SD of deconv (Gaussian fit to histogram core)
    %   mu_bn      centre of that Gaussian (should be ~0)
    %   sigma_mad  MAD-based SD of deconv (cross-check of sigma_bn)
    %   theta      detection threshold = mu_bn + theta_sigma * sigma_bn
    %   ev_frames  [1 x K] frames of the events (local maxima above theta)
    %   ev_amp     [1 x K] deconv value at those frames
    %   band_hz    band actually used (after the Nyquist cap)
    %   a, b       decay / rise coefficients exp(-frame_ms / tau) (b = 0: AR(1))

    trace = double(trace(:));
    nF = numel(trace);
    out = struct('filtered', nan(nF, 1), 'deconv', nan(nF, 1), 'sigma_bn', NaN, ...
        'mu_bn', NaN, 'sigma_mad', NaN, 'theta', NaN, 'ev_frames', zeros(1, 0), ...
        'ev_amp', zeros(1, 0), 'band_hz', [NaN NaN], 'a', NaN, 'b', NaN);
    if nargin < 5 || isempty(tau_rise_ms), tau_rise_ms = 0; end

    if ~(isfinite(frame_ms) && frame_ms > 0 && isfinite(tau_ms) && tau_ms > 0)
        error('deconv_events:timing', 'frame_ms and tau_ms must be positive (got %g, %g).', ...
            frame_ms, tau_ms);
    end
    if ~(isfinite(tau_rise_ms) && tau_rise_ms >= 0 && tau_rise_ms < tau_ms)
        error('deconv_events:rise', 'tau_rise_ms must be in [0, tau_ms) (got %g, tau %g).', ...
            tau_rise_ms, tau_ms);
    end

    fs   = 1000 / frame_ms;
    nyq  = fs / 2;
    band = p.band_hz;
    band(2) = min(band(2), p.max_hi_frac_nyquist * nyq);
    if band(1) >= band(2)
        error('deconv_events:band', 'Band [%g %g] Hz is empty at %g ms frames (Nyquist %g Hz).', ...
            band(1), band(2), frame_ms, nyq);
    end
    out.band_hz = band;
    a = exp(-frame_ms / tau_ms);
    b = 0;
    if tau_rise_ms > 0, b = exp(-frame_ms / tau_rise_ms); end
    out.a = a; out.b = b;

    % --- Finite span (epoch NaN padding sits at the ends) -------------------
    fin = isfinite(trace);
    i0 = find(fin, 1, 'first'); i1 = find(fin, 1, 'last');
    if isempty(i0) || (i1 - i0 + 1) < 20, return; end
    idx = (i0:i1)';
    x = trace(idx);
    gaps = ~isfinite(x);
    if any(gaps)
        k = find(~gaps);
        x(gaps) = interp1(k, x(k), find(gaps), 'linear', 'extrap');
    end

    % --- Band-pass -> deconvolve -> band-pass -------------------------------
    xf = gauss_bandpass(x, fs, band);
    if b == 0
        d = [xf(1) * (1 - a); xf(2:end) - a * xf(1:end-1)];
    else
        tpk = tau_ms * tau_rise_ms / (tau_ms - tau_rise_ms) * log(tau_ms / tau_rise_ms);
        kpk = exp(-tpk / tau_ms) - exp(-tpk / tau_rise_ms);
        xe  = [xf(1); xf; xf(end)];                    % hold the ends for x[n-1], x[n+1]
        d   = kpk / (a - b) * (xe(3:end) - (a + b) * xe(2:end-1) + a * b * xe(1:end-2));
    end
    df = gauss_bandpass(d, fs, band);

    xf(gaps) = NaN; df(gaps) = NaN;
    out.filtered(idx) = xf;
    out.deconv(idx)   = df;

    % --- sigma_BN: Gaussian fit to the core of the all-points histogram -----
    v = df(isfinite(df));
    [out.mu_bn, out.sigma_bn, out.sigma_mad] = hist_gauss_sigma(v, p.hist_core_sigma, p.hist_bins);
    if ~isfinite(out.sigma_bn) || out.sigma_bn <= 0, return; end
    out.theta = out.mu_bn + p.theta_sigma * out.sigma_bn;

    % --- Events: local maxima above theta -----------------------------------
    dpk = df; dpk(~isfinite(dpk)) = -Inf;
    if ~any(dpk > out.theta), return; end
    mind = max(1, round(p.min_peak_distance_ms / frame_ms));
    [pk, loc] = findpeaks([-Inf; dpk; -Inf], 'MinPeakHeight', out.theta, ...
        'MinPeakDistance', mind);
    out.ev_frames = (loc(:)' - 1) + (i0 - 1);
    out.ev_amp    = pk(:)';
end

function y = gauss_bandpass(x, fs, band)
    % Gaussian-response band-pass via FFT on a mirror-padded copy. Each edge
    % is the -3 dB frequency: low-pass exp(-f^2/(2 s^2)) with s = fc/sqrt(ln 2);
    % the high-pass is 1 minus the same shape at the low edge.
    n   = numel(x);
    pad = min(n - 1, max(32, round(2 * fs / band(1))));   % ~2 cycles of the low edge
    xp  = [flipud(x(2:pad+1)); x; flipud(x(end-pad:end-1))];
    m   = numel(xp);
    f   = (0:m-1)' * fs / m;
    f   = min(f, fs - f);                                  % two-sided |f|
    s_lo = band(1) / sqrt(log(2));
    s_hi = band(2) / sqrt(log(2));
    H = (1 - exp(-f.^2 / (2 * s_lo^2))) .* exp(-f.^2 / (2 * s_hi^2));
    yp = real(ifft(fft(xp) .* H));
    y  = yp(pad + 1 : pad + n);
end

function [mu, sigma, sigma_mad] = hist_gauss_sigma(v, core_k, nbins)
    % Fit a Gaussian to the all-points histogram restricted to its core
    % (|x - median| <= core_k * MAD-sigma), as a parabola in log(counts).
    % The core keeps the positive event tail from widening the fit.
    mu = NaN; sigma = NaN;
    med = median(v);
    sigma_mad = mad(v, 1) * 1.4826;
    if ~(isfinite(sigma_mad) && sigma_mad > 0) || numel(v) < 50, return; end
    edges = linspace(med - core_k * sigma_mad, med + core_k * sigma_mad, nbins + 1);
    c = histcounts(v, edges)';
    xc = (edges(1:end-1) + edges(2:end))' / 2;
    ok = c > 0;
    if sum(ok) < 5, return; end
    % log c = b0 + b1 x + b2 x^2 ; weighted by counts (Poisson variance of log c ~ 1/c)
    X = [ones(sum(ok), 1), xc(ok), xc(ok).^2];
    w = sqrt(c(ok));
    b = (X .* w) \ (log(c(ok)) .* w);
    if b(3) >= 0, return; end                              % not a downward parabola
    sigma = sqrt(-1 / (2 * b(3)));
    mu    = -b(2) / (2 * b(3));
end
