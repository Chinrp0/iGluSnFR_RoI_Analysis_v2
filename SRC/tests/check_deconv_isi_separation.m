function res = check_deconv_isi_separation(folder, interval_ms, tau_ms, bands_hz, test_isis_ms, out_dir)
    % CHECK_DECONV_ISI_SEPARATION - Does deconv_events keep pulse 1 and pulse 2 apart at short ISI?
    %
    % A lower band-pass edge smooths the deconvolved trace, so pulse 1's event
    % can spill into the pulse-2 window. This measures that spill with a
    % control that has NO second stimulus: the long-ISI files (ISI >= 500 ms).
    % For each test ISI and each low-pass edge (AR(1) kernel):
    %   real  - ISI files: P(event in [stim2, stim2 + peak_search_ms] | event in
    %           [stim1, stim2 - 1])
    %   leak  - long-ISI files, same windows placed at stim1 + ISI although no
    %           stim2 comes: P(event in pseudo pulse-2 window | event in pulse-1
    %           window) = pulse-1 spill + chance
    %   sep   - real - leak (what stim2 adds). A band that separates the pulses
    %           keeps leak near the chance level and sep high.
    %   valley - grand-average deconvolved trace (in sigma_BN units): the
    %           minimum between the pulse-1 and pulse-2 peaks / the smaller peak.
    %           0 = fully resolved, 1 = merged into one hump.
    %   iei   - median interval between the first pulse-1 and first pulse-2
    %           event in ROIs with both (should track the ISI)
    %
    % Usage:
    %   res = check_deconv_isi_separation(folder, 10, 42.8)
    %   res = check_deconv_isi_separation(folder, 5, 30.1, [10 15 20 30], [30 60], out_dir)
    %
    % Inputs:
    %   folder       - PPF raw-mean folder (with Recordings.csv)
    %   interval_ms  - frame-interval set to use ([] = main set)
    %   tau_ms       - unitary decay tau (ms); [] = config.deconv for the sensor
    %   bands_hz     - low-pass edges to compare (default [10 15 20 30])
    %   test_isis_ms - ISIs to test (default [30 60])
    %   out_dir      - optional: write the table + figure here

    if nargin < 2, interval_ms = []; end
    if nargin < 3, tau_ms = []; end
    if nargin < 4 || isempty(bands_hz), bands_hz = [10 15 20 30]; end
    if nargin < 5 || isempty(test_isis_ms), test_isis_ms = [30 60]; end
    if nargin < 6, out_dir = ''; end
    evalc('config = tracenorm_config();');
    config.dataType = '2AP';
    p = config.deconv;

    [~, ds_all] = evalc('load_dataset_folder(folder, config)');
    parts = split_dataset_by_interval(ds_all, config);
    if isempty(interval_ms)
        ds = parts(end);
    else
        [~, k] = min(abs([parts.frame_interval_ms] - interval_ms));
        ds = parts(k);
    end
    info = experiment_info(folder, ds.files, config);
    fd = ds.frame_interval_ms;
    if isempty(tau_ms)
        for k = 1:numel(p.tau_ms)
            if ~isempty(regexpi([info.sensor ' ' ds.files(1).name], p.tau_ms(k).pattern, 'once'))
                tau_ms = p.tau_ms(k).tau_ms; break;
            end
        end
    end
    if isempty(tau_ms) || ~isfinite(tau_ms)
        error('check_deconv_isi_separation:tau', 'Pass tau_ms or set config.deconv.tau_ms for %s.', info.sensor);
    end
    fprintf('\n=== ISI separation: %s | %.4g ms frames | tau %.1f ms | AR(1) | low-pass %s Hz ===\n', ...
        info.sensor, fd, tau_ms, num2str(bands_hz));

    %% Real pipeline per file (same epochs / dF/F as the analysis)
    fr = {};
    for i = 1:numel(ds.files)
        try
            fr{end+1} = twoap_pipeline(ds.files(i), struct('verbose', false), config); %#ok<AGROW>
        catch ME
            fprintf('  skip %s: %s\n', ds.files(i).name, ME.message);
        end
    end
    isi  = cellfun(@(x) x.isi_ms, fr);
    s1   = fr{1}.twoap.stim1_frame;
    win_f = round(config.twoAP.peak_search_ms / fd);
    long = find(isi >= 500);
    span = round((max(test_isis_ms) + 150) / fd);         % grand-average window after stim1
    pre_show = round(50 / fd);

    rows = {};
    avg = struct('band', {}, 'isi', {}, 'real', {}, 'ctrl', {});
    for hi = bands_hz
        pb = p; pb.band_hz = [p.band_hz(1) hi];
        % Deconvolve every analysed ROI once per band
        D = cell(numel(fr), 1);                           % {file} -> [nF x ROIs] in sigma_BN units
        E = cell(numel(fr), 1);                           % {file}{roi} -> event frames
        for i = 1:numel(fr)
            t = fr{i}.twoap;
            nF = size(fr{i}.dfof_data, 1);
            D{i} = nan(nF, fr{i}.numROIs); E{i} = cell(1, fr{i}.numROIs);
            for r = 1:fr{i}.numROIs
                if ~isfinite(t.prestim_sd(r)) || t.noise_excluded(r), continue; end
                d = deconv_events(fr{i}.dfof_data(:, r), fd, tau_ms, pb, 0);
                if ~isfinite(d.sigma_bn), continue; end
                D{i}(:, r) = d.deconv / d.sigma_bn;
                E{i}{r} = d.ev_frames;
            end
        end
        for isi_ms = test_isis_ms
            f_isi = find(abs(isi - isi_ms) < 1);
            if isempty(f_isi), continue; end
            isi_f = round(isi_ms / fd);
            [real_c, n_real, iei_real, gR] = window_stats(fr(f_isi), E(f_isi), D(f_isi), s1, isi_f, win_f, span, pre_show);
            [leak_c, n_ctrl, ~, gC]        = window_stats(fr(long), E(long), D(long), s1, isi_f, win_f, span, pre_show);
            vR = valley(gR, pre_show, isi_f, win_f);
            vC = valley(gC, pre_show, isi_f, win_f);
            rows(end+1, :) = {hi, isi_ms, n_real.p1, real_c, n_ctrl.p1, leak_c, real_c - leak_c, ...
                vR, vC, iei_real * fd, n_real.p2_all, n_ctrl.p2_all}; %#ok<AGROW>
            avg(end+1) = struct('band', hi, 'isi', isi_ms, 'real', gR, 'ctrl', gC); %#ok<AGROW>
        end
    end
    T = cell2table(rows, 'VariableNames', {'lowpass_hz', 'isi_ms', 'n_p1_real', 'p2_given_p1_real', ...
        'n_p1_ctrl', 'p2_given_p1_leak', 'sep', 'valley_real', 'valley_ctrl', 'iei_median_ms', ...
        'p2_all_real', 'p2_all_ctrl'});
    disp(T);
    fprintf(['  sep = what stim2 adds over pulse-1 spill (higher = better separated)\n' ...
             '  valley: 0 = resolved, 1 = one hump; valley_ctrl is the no-stim2 baseline\n']);

    res = struct('sensor', info.sensor, 'frame_ms', fd, 'tau_ms', tau_ms, 'table', T, 'avg', avg);

    %% Figure: grand-average deconvolved traces, real vs no-stim2 control
    fig = figure('Name', 'Deconv ISI separation', 'Color', 'w', 'Position', [60 60 1500 650]);
    tl = tiledlayout(fig, numel(test_isis_ms), numel(bands_hz), 'TileSpacing', 'compact');
    title(tl, sprintf('%s | %.4g ms frames | \\tau %.1f ms | mean deconv / \\sigma_{BN}: ISI files (colour) vs ISI >= 500 ms control (grey)', ...
        info.sensor, fd, tau_ms));
    for a = 1:numel(test_isis_ms)
        for b = 1:numel(bands_hz)
            nexttile;
            k = find([avg.band] == bands_hz(b) & [avg.isi] == test_isis_ms(a), 1);
            if isempty(k), axis off; continue; end
            x = ((1:numel(avg(k).real)) - pre_show - 1) * fd;
            plot(x, avg(k).ctrl, 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2); hold on;
            plot(x, avg(k).real, 'Color', [0.1 0.45 0.85], 'LineWidth', 1.5);
            xline([0 test_isis_ms(a)], 'k:');
            r = T(T.lowpass_hz == bands_hz(b) & T.isi_ms == test_isis_ms(a), :);
            title(sprintf('ISI %g | LP %g Hz | sep %.2f | valley %.2f', test_isis_ms(a), bands_hz(b), ...
                r.sep, r.valley_real), 'FontSize', 9);
            if b == 1, ylabel('mean deconv (\sigma_{BN})'); end
            if a == numel(test_isis_ms), xlabel('ms from stim1'); end
            xlim([x(1), x(end)]);
        end
    end
    if ~isempty(out_dir)
        if ~isfolder(out_dir), mkdir(out_dir); end
        tag = sprintf('%s_%gms', regexprep(info.sensor, '\W', ''), round(fd));
        writetable(T, fullfile(out_dir, sprintf('deconv_isi_sep_%s.csv', tag)));
        exportgraphics(fig, fullfile(out_dir, sprintf('deconv_isi_sep_%s.png', tag)), 'Resolution', 150);
    end
end

function [p2_given_p1, n, iei_med, g] = window_stats(fr, E, D, s1, isi_f, win_f, span, pre_show)
    % Pulse-1 window [s1, s1+isi_f-1], pulse-2 window [s1+isi_f, s1+isi_f+win_f]
    w1 = [s1, s1 + isi_f - 1];
    w2 = [s1 + isi_f, s1 + isi_f + win_f];
    has1 = []; has2 = []; iei = []; G = [];
    for i = 1:numel(fr)
        for r = 1:numel(E{i})
            if all(isnan(D{i}(:, r))), continue; end
            ev = E{i}{r};
            e1 = ev(ev >= w1(1) & ev <= w1(2));
            e2 = ev(ev >= w2(1) & ev <= w2(2));
            has1(end+1) = ~isempty(e1); has2(end+1) = ~isempty(e2); %#ok<AGROW>
            if ~isempty(e1) && ~isempty(e2), iei(end+1) = e2(1) - e1(1); end %#ok<AGROW>
        end
        G = [G, D{i}(s1 - pre_show : min(s1 + span, end), :)]; %#ok<AGROW>
    end
    has1 = logical(has1); has2 = logical(has2);
    n = struct('p1', sum(has1), 'p2_all', mean(has2));
    p2_given_p1 = mean(has2(has1));
    iei_med = median(iei);
    g = mean(G, 2, 'omitnan');
end

function v = valley(g, pre_show, isi_f, win_f)
    % min between the pulse-1 and pulse-2 peaks / the smaller peak
    o = pre_show + 1;                                     % stim1 index in g
    [p1, k1] = max(g(o : o + isi_f - 1));     k1 = o + k1 - 1;
    [p2, k2] = max(g(o + isi_f : min(o + isi_f + win_f, end))); k2 = o + isi_f + k2 - 1;
    if k2 <= k1 + 1 || min(p1, p2) <= 0, v = NaN; return; end
    v = max(0, min(g(k1:k2))) / min(p1, p2);
end
