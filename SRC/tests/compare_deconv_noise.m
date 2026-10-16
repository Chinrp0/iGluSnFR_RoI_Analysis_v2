function cmp = compare_deconv_noise(folder, interval_ms, out_dir, kernel, p_override)
    % COMPARE_DECONV_NOISE - Deconvolution noise floor / events vs the current 2AP gates
    %
    % Runs every file of one PPF folder through twoap_pipeline (same epochs,
    % F0 and dF/F as the real run), then deconv_events on each ROI's dF/F
    % epoch, and compares per ROI:
    %   - noise: sigma_BN (deconvolved histogram fit) vs prestim MAD SD (dF/F)
    %   - pulse-1 responder: deconvolved event in [stim1, stim1 + peak_search_ms]
    %     vs the Schmitt gate (is_responder1)
    %   - pulse-2 responder: deconvolved event in [stim2, stim2 + peak_search_ms]
    %     vs is_responder2, by ISI (short ISIs are where deconvolution should help)
    %   - spontaneous events in the pre-stim1 baseline vs twoap n_spont
    %   - chance level: a deconvolved "response" in a SHAM window of the same
    %     length as the pulse-1 window, placed in the pre-stim baseline
    %     (sham_ms before stim1). Its rate is the deconv false-positive rate.
    % The unitary kernel comes from config.deconv.tau_ms for the folder's
    % sensor; NaN values are measured here from the grand-average pulse-1
    % response of the longest-ISI files:
    %   'ar1' - instantaneous rise, decay tau = single exponential from the peak
    %   'ar2' - rise + decay: A*(exp(-t/tau_d) - exp(-t/tau_r)) fitted from the
    %           stim (onset delay free), both taus from that one fit. Default 'ar1':
    %           on the 2026-10-15 PPF data AR(2) detected no more (usually fewer)
    %           responders; the rise is <= 2 frames and its inverse adds noise.
    %
    % Usage:
    %   cmp = compare_deconv_noise(folder)            % main frame-interval set
    %   cmp = compare_deconv_noise(folder, 10)        % only the 10 ms files
    %   cmp = compare_deconv_noise(folder, 10, dir)   % also write CSV + figures to dir
    %   cmp = compare_deconv_noise(folder, 10, dir, 'ar1')   % instantaneous-rise kernel
    %   cmp = compare_deconv_noise(folder, 10, '', 'ar1', struct('band_hz', [0.5 15]))
    %                                     % override config.deconv fields for a sweep
    %
    % Output: cmp struct with .sensor .frame_ms .kernel .tau_ms .tau_rise_ms .tau_fit .table (one
    % row per analysed ROI) .summary (printed) .example (traces for figures)

    if nargin < 2, interval_ms = []; end
    if nargin < 3, out_dir = ''; end
    if nargin < 4 || isempty(kernel), kernel = 'ar1'; end
    kernel = lower(kernel);
    if ~any(strcmp(kernel, {'ar1', 'ar2'}))
        error('compare_deconv_noise:kernel', 'kernel must be ''ar1'' or ''ar2'' (got %s).', kernel);
    end
    evalc('config = tracenorm_config();');
    config.dataType = '2AP';
    p = config.deconv;
    if nargin >= 5 && ~isempty(p_override)
        for fn = fieldnames(p_override)'
            p.(fn{1}) = p_override.(fn{1});
        end
    end

    [~, ds_all] = evalc('load_dataset_folder(folder, config)');
    parts = split_dataset_by_interval(ds_all, config);
    if isempty(interval_ms)
        ds = parts(end);                                  % main set (most files)
    else
        fis = [parts.frame_interval_ms];
        [dmin, k] = min(abs(fis - interval_ms));
        if dmin > 0.05 * interval_ms
            error('compare_deconv_noise:interval', 'No %g ms set (have %s ms).', ...
                interval_ms, num2str(fis, '%.4g '));
        end
        ds = parts(k);
    end
    info = experiment_info(folder, ds.files, config);
    frame_ms = ds.frame_interval_ms;
    fprintf('\n=== Deconvolution comparison: %s | %s | %.4g ms frames | %d files | %s kernel ===\n', ...
        info.sensor, info.date, frame_ms, numel(ds.files), upper(kernel));

    %% 1. Real pipeline per file (epoch -> F0 -> dF/F -> twoap_analyzer)
    fr = cell(numel(ds.files), 1);
    for i = 1:numel(ds.files)
        try
            fr{i} = twoap_pipeline(ds.files(i), struct('verbose', false), config);
        catch ME
            fprintf('  skip %s: %s\n', ds.files(i).name, ME.message);
        end
    end
    fr = fr(~cellfun(@isempty, fr));
    if isempty(fr), error('compare_deconv_noise:noFiles', 'No file processed.'); end
    s1 = fr{1}.twoap.stim1_frame;

    %% 2. Unitary kernel for this sensor (config, else measured)
    tc = sensor_tau(p.tau_ms, [info.sensor ' ' ds.files(1).name]);
    tau_fit = measure_tau(fr, s1, frame_ms, config);
    if strcmp(kernel, 'ar1')
        tau_ms = tc.tau_ms; tau_rise_ms = 0; src = 'config';
        if ~isfinite(tau_ms), tau_ms = tau_fit.tau_ms; src = 'MEASURED'; end
        fprintf('  decay tau %.1f ms (%s; single exp from peak, R^2 %.2f), rise instantaneous\n', ...
            tau_ms, src, tau_fit.r2);
    else
        tau_ms = tc.tau_ms; tau_rise_ms = tc.tau_rise_ms; src = 'config';
        if ~isfinite(tau_ms) || ~isfinite(tau_rise_ms)
            tau_ms = tau_fit.rd_decay_ms; tau_rise_ms = tau_fit.rd_rise_ms; src = 'MEASURED';
        end
        fprintf('  rise tau %.1f ms, decay tau %.1f ms (%s; rise-decay fit from stim, onset %.1f ms, R^2 %.2f)\n', ...
            tau_rise_ms, tau_ms, src, tau_fit.rd_onset_ms, tau_fit.rd_r2);
    end
    fprintf('  (grand-average pulse 1, ISI >= %g ms, %d ROIs; peak-only single exp %.1f ms)\n', ...
        tau_fit.min_isi_ms, tau_fit.n_rois, tau_fit.tau_ms);
    if strcmp(src, 'MEASURED')
        fprintf('       -> set config.deconv.tau_ms for this sensor to pin it.\n');
    end
    if ~isfinite(tau_ms) || ~isfinite(tau_rise_ms)
        error('compare_deconv_noise:tau', 'Could not measure the kernel; set config.deconv.tau_ms.');
    end
    a = exp(-frame_ms / tau_ms);
    b = 0; if tau_rise_ms > 0, b = exp(-frame_ms / tau_rise_ms); end
    fprintf('  frames %.4g ms (Nyquist %.0f Hz) -> band %.2g-%.3g Hz, a = %.3f, b = %.3f\n', ...
        frame_ms, 500 / frame_ms, p.band_hz(1), min(p.band_hz(2), p.max_hi_frac_nyquist * 500 / frame_ms), a, b);

    %% 3. Deconvolve every analysed ROI
    win_f = round(config.twoAP.peak_search_ms / frame_ms);
    sham_ms = 500;                                        % sham window this far before stim1
    rows = {};
    ex = struct('trace', {}, 'deconv', {}, 'theta', {}, 'ev', {}, 'name', {}, 's1', {}, 's2', {}, ...
        'resp_schmitt', {}, 'resp_deconv', {});
    for i = 1:numel(fr)
        t = fr{i}.twoap;
        s2 = t.stim2_frame;
        pre_end = t.prestim_idx(end);
        w1 = [s1, min(s1 + win_f, s2 - 1)];
        w2 = [s2, s2 + win_f];
        ws = w1 - round(sham_ms / frame_ms);              % sham: same length, in the baseline
        for r = 1:fr{i}.numROIs
            if ~isfinite(t.prestim_sd(r)) || t.noise_excluded(r), continue; end
            d = deconv_events(fr{i}.dfof_data(:, r), frame_ms, tau_ms, p, tau_rise_ms);
            if ~isfinite(d.sigma_bn), continue; end
            ev = d.ev_frames;
            in1 = ev >= w1(1) & ev <= w1(2);
            in2 = ev >= w2(1) & ev <= w2(2);
            ins = ev >= ws(1) & ev <= ws(2);
            e1 = max([d.ev_amp(in1), NaN]);               % NaN when none
            pre = d.deconv(1:pre_end); pre = pre(isfinite(pre));
            neg_tail = mean(pre < d.mu_bn - p.theta_sigma * d.sigma_bn);
            rows(end+1, :) = {string(fr{i}.filename), string(fr{i}.genotype), fr{i}.isi_ms, r, ...
                t.prestim_sd(r), d.sigma_bn, d.sigma_mad, mad(pre, 1) * 1.4826, ...
                t.peak1_amp(r), t.peak1_amp(r) / t.prestim_sd(r), e1, e1 / d.sigma_bn, ...
                t.is_responder1(r), any(in1), any(ins), t.is_responder2(r), any(in2), ...
                sum(ev <= pre_end), t.n_spont(r), neg_tail}; %#ok<AGROW>
            if numel(ex) < 400
                ex(end+1) = struct('trace', fr{i}.dfof_data(:, r), 'deconv', d.deconv, ...
                    'theta', d.theta, 'ev', ev, 'name', sprintf('%s r%d', fr{i}.filename, r), ...
                    's1', s1, 's2', s2, 'resp_schmitt', t.is_responder1(r), 'resp_deconv', any(in1)); %#ok<AGROW>
            end
        end
    end
    T = cell2table(rows, 'VariableNames', {'file', 'genotype', 'isi_ms', 'roi', ...
        'prestim_sd', 'sigma_bn', 'sigma_mad_deconv', 'sigma_mad_deconv_pre', ...
        'A1', 'snr_A1', 'deconv_e1', 'snr_deconv_e1', ...
        'resp1_schmitt', 'resp1_deconv', 'sham_deconv', 'resp2_schmitt', 'resp2_deconv', ...
        'n_spont_deconv', 'n_spont_schmitt', 'neg_tail_frac'});

    %% 4. Summary
    S = struct();
    S.n_rois = height(T);
    rr = corr(T.prestim_sd, T.sigma_bn, 'type', 'Spearman');
    ratio = T.sigma_bn ./ T.prestim_sd;
    S.noise_spearman = rr;
    S.sigma_ratio_median = median(ratio);
    S.sigma_ratio_iqr = iqr(ratio);
    S.sigma_bn_vs_mad = median(T.sigma_bn ./ T.sigma_mad_deconv);
    S.sigma_bn_full_vs_pre = median(T.sigma_bn ./ T.sigma_mad_deconv_pre);
    S.neg_tail_mean = mean(T.neg_tail_frac);
    S.resp1 = crosstab_counts(T.resp1_schmitt, T.resp1_deconv);
    S.spont_schmitt = sum(T.n_spont_schmitt);      % NaN when twoAP events are off (frames > max_frame_ms)
    S.spont_deconv  = sum(T.n_spont_deconv);

    fprintf('\n  %d ROIs analysed\n', S.n_rois);
    fprintf('  NOISE  sigma_BN vs prestim MAD SD: Spearman %.2f | ratio median %.2f (IQR %.2f)\n', ...
        rr, S.sigma_ratio_median, S.sigma_ratio_iqr);
    fprintf('         (a constant ratio = same ranking of ROIs, no new QC information)\n');
    fprintf('         sigma_BN / MAD of deconv: %.2f (1 = Gaussian core) | full trace / pre-stim only: %.2f\n', ...
        S.sigma_bn_vs_mad, S.sigma_bn_full_vs_pre);
    fprintf('         negative tail beyond -%g sigma_BN: %.2g%% of pre-stim frames (Gaussian: %.2g%%)\n', ...
        p.theta_sigma, 100 * S.neg_tail_mean, 100 * normcdf(-p.theta_sigma));
    fprintf('  PULSE 1 responders (Schmitt rows x deconv cols):\n');
    fprintf('                 deconv-   deconv+\n');
    fprintf('    Schmitt-    %6d    %6d\n', S.resp1(1, 1), S.resp1(1, 2));
    fprintf('    Schmitt+    %6d    %6d\n', S.resp1(2, 1), S.resp1(2, 2));
    fprintf('    agreement %.1f%% | Schmitt %d vs deconv %d responders\n', ...
        100 * sum(diag(S.resp1)) / max(S.n_rois, 1), sum(S.resp1(2, :)), sum(S.resp1(:, 2)));
    S.sham_rate = mean(T.sham_deconv);
    fprintf('    deconv chance level (sham window %g ms before stim1): %.2f%% of ROIs -> ~%.0f expected false responders\n', ...
        sham_ms, 100 * S.sham_rate, S.sham_rate * S.n_rois);
    only_d = T.resp1_deconv & ~T.resp1_schmitt;
    fprintf('    median A1 / prestim SD: deconv-only %.2f | Schmitt-only %.2f | both %.2f\n', ...
        median(T.snr_A1(only_d)), median(T.snr_A1(T.resp1_schmitt & ~T.resp1_deconv)), ...
        median(T.snr_A1(T.resp1_schmitt & T.resp1_deconv)));

    isis = unique(T.isi_ms)';
    fprintf('  PULSE 2 responders by ISI (fraction of analysed ROIs):\n');
    fprintf('    ISI(ms)   n    Schmitt   deconv   P1 Schmitt   P1 deconv\n');
    S.by_isi = zeros(numel(isis), 6);
    for k = 1:numel(isis)
        m = T.isi_ms == isis(k);
        S.by_isi(k, :) = [isis(k), sum(m), mean(T.resp2_schmitt(m)), mean(T.resp2_deconv(m)), ...
            mean(T.resp1_schmitt(m)), mean(T.resp1_deconv(m))];
        fprintf('    %6g  %4d    %6.2f   %6.2f     %6.2f      %6.2f\n', S.by_isi(k, :));
    end
    pre_s = (fr{1}.twoap.prestim_idx(end)) * frame_ms / 1000;
    fprintf('  SPONT events in pre-stim1 (%.2f s/ROI): Schmitt %d (%.3f Hz/ROI) | deconv %d (%.3f Hz/ROI)\n', ...
        pre_s, S.spont_schmitt, S.spont_schmitt / (S.n_rois * pre_s), ...
        S.spont_deconv, S.spont_deconv / (S.n_rois * pre_s));

    cmp = struct('sensor', info.sensor, 'date', info.date, 'frame_ms', frame_ms, ...
        'kernel', kernel, 'tau_ms', tau_ms, 'tau_rise_ms', tau_rise_ms, 'tau_fit', tau_fit, 'table', T, 'summary', S, 'example', ex);

    %% 5. Figures (+ files when out_dir is given)
    figs = make_figures(cmp, config, p);
    if ~isempty(out_dir)
        if ~isfolder(out_dir), mkdir(out_dir); end
        tag = sprintf('%s_%gms_%s', regexprep(info.sensor, '\W', ''), round(frame_ms), kernel);
        writetable(T, fullfile(out_dir, sprintf('deconv_compare_%s.csv', tag)));
        for k = 1:numel(figs)
            exportgraphics(figs(k), fullfile(out_dir, sprintf('deconv_%s_%d.png', tag, k)), ...
                'Resolution', 150);
        end
        fprintf('  wrote table + %d figures to %s\n', numel(figs), out_dir);
    end
end

function tc = sensor_tau(tab, name)
    tc = struct('tau_ms', NaN, 'tau_rise_ms', NaN);
    for k = 1:numel(tab)
        if ~isempty(regexpi(name, tab(k).pattern, 'once'))
            tc.tau_ms = tab(k).tau_ms; tc.tau_rise_ms = tab(k).tau_rise_ms; return;
        end
    end
end

function f = measure_tau(fr, s1, frame_ms, config)
    % Single exponential A*exp(-t/tau) + C fitted from the grand-average
    % pulse-1 peak, using only files whose ISI leaves >= 300 ms of clean decay.
    isis = cellfun(@(x) x.isi_ms, fr);
    min_isi = min(500, max(isis));
    tr = [];
    for i = find(isis >= min_isi)'
        t = fr{i}.twoap;
        ok = isfinite(t.prestim_sd) & ~t.noise_excluded;
        y = fr{i}.dfof_data(:, ok) - t.prestim_base(ok);
        tr = [tr, y]; %#ok<AGROW>
    end
    f = struct('tau_ms', NaN, 'r2', NaN, 'n_rois', size(tr, 2), 'min_isi_ms', min_isi, ...
        'tau_fast', NaN, 'tau_slow', NaN);
    if isempty(tr), return; end
    g = double(mean(tr, 2, 'omitnan'));
    pk_f = round(config.twoAP.peak_search_ms / frame_ms);
    [~, k] = max(g(s1 : s1 + pk_f));
    pk = s1 + k - 1;
    n = min(round(min(min_isi - 20, 400) / frame_ms), numel(g) - pk);
    t = (0:n)' * frame_ms;
    y = g(pk : pk + n);
    sse = @(b) sum((b(1) * exp(-t / b(2)) + b(3) - y).^2);
    b = fminsearch(sse, [y(1), 30, 0], optimset('MaxFunEvals', 4000, 'MaxIter', 4000));
    f.tau_ms = b(2);
    f.r2 = 1 - sse(b) / sum((y - mean(y)).^2);
    f.t = t; f.y = y; f.fit = b(1) * exp(-t / b(2)) + b(3);
    f.peak_latency_ms = (pk - s1) * frame_ms;
    [f.tau_fast, f.tau_slow] = fit_decay_double(t, y);

    % Rise + decay from the stim: C for t < t0, else
    %   A * (exp(-(t-t0)/tau_d) - exp(-(t-t0)/tau_r)) / kpk + C   (A = peak height)
    % tau_d > tau_r enforced by tau_d = tau_r + exp(q); t0 in (0, peak latency).
    tt = (0:(pk - s1 + n))' * frame_ms;
    yy = g(s1 : pk + n);
    lat = max(f.peak_latency_ms, frame_ms);
    model = @(q, x) rd_model(q, x, lat);
    sse2 = @(q) sum((model(q, tt) - yy).^2);
    q0 = [max(yy) - min(yy), log(max(frame_ms, 3)), log(max(f.tau_ms, 5)), 0, 0];
    q = fminsearch(sse2, q0, optimset('MaxFunEvals', 20000, 'MaxIter', 20000));
    [~, tr, td, t0] = model(q, tt);
    f.rd_rise_ms = tr; f.rd_decay_ms = td; f.rd_onset_ms = t0;
    f.rd_r2 = 1 - sse2(q) / sum((yy - mean(yy)).^2);
    f.rd_t = tt; f.rd_y = yy; f.rd_fit = model(q, tt);
end

function [y, tr, td, t0] = rd_model(q, t, lat)
    tr = exp(q(2));
    td = tr + exp(q(3));
    t0 = lat / (1 + exp(-q(4)));                          % logistic -> (0, lat)
    tpk = td * tr / (td - tr) * log(td / tr);
    kpk = exp(-tpk / td) - exp(-tpk / tr);
    s = max(t - t0, 0);
    y = q(1) * (exp(-s / td) - exp(-s / tr)) / kpk + q(5);
end

function c = crosstab_counts(a, b)
    a = logical(a); b = logical(b);
    c = [sum(~a & ~b), sum(~a & b); sum(a & ~b), sum(a & b)];
end

function figs = make_figures(cmp, config, p)
    T = cmp.table; fd = cmp.frame_ms;
    ttl = sprintf('%s %s | %.4g ms frames | %s: \\tau_{rise} %.1f, \\tau_{decay} %.1f ms', ...
        cmp.sensor, cmp.date, fd, upper(cmp.kernel), cmp.tau_rise_ms, cmp.tau_ms);

    % 1: tau fit + noise comparison + an example histogram
    figs(1) = figure('Name', 'Deconv: tau + noise', 'Color', 'w', 'Position', [80 80 1300 380]);
    tl = tiledlayout(figs(1), 1, 3); title(tl, ttl);
    nexttile; f = cmp.tau_fit;
    if isfield(f, 'rd_t')
        plot(f.rd_t, f.rd_y, 'k.'); hold on;
        plot(f.rd_t, f.rd_fit, 'r-', 'LineWidth', 1.5);
        plot(f.t + f.peak_latency_ms, f.fit, 'b--', 'LineWidth', 1);
        legend({'grand avg', sprintf('rise-decay: \\tau_r %.1f, \\tau_d %.1f (R^2 %.2f)', ...
            f.rd_rise_ms, f.rd_decay_ms, f.rd_r2), sprintf('peak exp: \\tau %.1f', f.tau_ms)}, ...
            'Location', 'northeast', 'Box', 'off');
        title(sprintf('Pulse-1 kernel (onset %.1f ms)', f.rd_onset_ms));
        xlim([0, min(f.rd_t(end), 250)]);
    end
    xlabel('ms from stim1'); ylabel('grand-average dF/F');
    nexttile;
    scatter(T.prestim_sd, T.sigma_bn, 8, double(T.resp1_schmitt), 'filled');
    xlabel('pre-stim MAD SD (dF/F)'); ylabel('\sigma_{BN} (deconvolved)');
    title(sprintf('Spearman %.2f | colour = Schmitt responder', cmp.summary.noise_spearman));
    colormap(gca, [0.6 0.6 0.6; 0.85 0.2 0.2]);
    nexttile;
    [~, k] = min(abs(T.sigma_bn - median(T.sigma_bn)));
    e = cmp.example(min(k, numel(cmp.example)));
    v = e.deconv(isfinite(e.deconv));
    histogram(v, 80, 'Normalization', 'pdf', 'EdgeColor', 'none'); hold on;
    xs = linspace(min(v), max(v), 200);
    sig = (e.theta) / p.theta_sigma;                       % approx (mu ~ 0)
    plot(xs, normpdf(xs, 0, sig), 'r-', 'LineWidth', 1.5);
    xline(e.theta, 'r--', sprintf('%g\\sigma_{BN}', p.theta_sigma));
    set(gca, 'YScale', 'log'); xlabel('deconvolved value'); ylabel('pdf');
    title('All-points histogram (median-noise ROI)', 'Interpreter', 'none');

    % 2: example traces, where the two gates disagree first
    dis = find([cmp.example.resp_schmitt] ~= [cmp.example.resp_deconv]);
    agr = find([cmp.example.resp_schmitt] & [cmp.example.resp_deconv]);
    pick = [dis(1:min(6, end)), agr(1:min(6 - min(6, numel(dis)) + 2, end))];
    pick = pick(1:min(8, end));
    figs(2) = figure('Name', 'Deconv: examples', 'Color', 'w', 'Position', [80 80 1400 800]);
    tl = tiledlayout(figs(2), numel(pick), 1, 'TileSpacing', 'tight'); title(tl, [ttl ' | examples (disagreements first)']);
    for j = 1:numel(pick)
        e = cmp.example(pick(j));
        x = ((1:numel(e.trace))' - e.s1) * fd;
        keep = x >= -300 & x <= (e.s2 - e.s1) * fd + 400;
        nexttile;
        yyaxis left; plot(x(keep), e.trace(keep), 'Color', [0.2 0.2 0.2]); ylabel('dF/F');
        set(gca, 'YColor', [0.2 0.2 0.2]);
        yyaxis right; plot(x(keep), e.deconv(keep), 'Color', [0.1 0.45 0.85]); hold on;
        set(gca, 'YColor', [0.1 0.45 0.85]); ylabel('deconv');
        yline(e.theta, ':', 'Color', [0.1 0.45 0.85]);
        evx = (e.ev - e.s1) * fd; evk = evx >= -300 & evx <= max(x(keep));
        plot(evx(evk), e.deconv(e.ev(evk)), 'v', 'Color', [0.85 0.2 0.2], 'MarkerFaceColor', [0.85 0.2 0.2]);
        xline([0, (e.s2 - e.s1) * fd], 'k:');
        title(sprintf('%s | Schmitt %d, deconv %d', e.name, e.resp_schmitt, e.resp_deconv), ...
            'Interpreter', 'none', 'FontSize', 8);
        if j < numel(pick), set(gca, 'XTickLabel', []); end
    end
    xlabel('ms from stim1');
end
