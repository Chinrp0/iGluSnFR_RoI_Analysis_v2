function results = hfs_batch_analysis(dataset, options)
    % HFS_BATCH_ANALYSIS - Batch analysis of high-frequency train (HFS) data
    %
    % Takes a dataset from load_dataset_folder (one frame interval; run_pipeline
    % splits mixed folders first). Every file must share one train (Train_Pulses
    % x Train_Hz). Files are grouped by genotype x drug and processed through
    % hfs_pipeline. For each coverslip the QC-passed ROIs of all its files are
    % averaged into one grand-average dF/F trace (no responder gating), and the
    % PRIMARY metrics are measured on that trace with hfs_train_metrics:
    %   per pulse  - phasic amplitude A_n, A_n / A_1, tonic level, peak latency
    %   train      - charge, plateau, peak, steady-state amplitude / ratio,
    %                cumulative-amplitude pool estimate + replenishment per pulse,
    %                post-train decay tau and half-decay time
    % The coverslip is the statistical unit (regions / files / ROIs within it are
    % pseudoreplicates). Per-ROI values are kept as diagnostics in the ROI table.
    %
    % When single pulses are not resolvable at this frame interval (fewer than
    % config.hfs.min_frames_per_pulse frames per pulse, e.g. 20 Hz at 25 ms) only
    % the train-level metrics are reported and the per-pulse figure is skipped.
    %
    % Inputs:
    %   dataset - struct from load_dataset_folder (one frame interval)
    %   options - struct: verbose, createPlots
    %
    % Output:
    %   results - struct with .groups, .plot_handles, .hfs_roi_table,
    %             .hfs_pulse_table, .hfs_coverslip_table, .file_checks,
    %             .processing_info

    if nargin < 2, options = struct(); end
    if ~isfield(options, 'verbose'),     options.verbose = true; end
    if ~isfield(options, 'createPlots'), options.createPlots = true; end
    if ~isstruct(dataset) || ~isfield(dataset, 'files')
        error('hfs_batch_analysis:input', 'hfs_batch_analysis takes a dataset from load_dataset_folder.');
    end
    evalc('config = tracenorm_config();');
    config.dataType = 'HFS';

    fprintf('=== HFS (high-frequency train) Batch Analysis ===\nFolder: %s\n', dataset.folder);

    %% Step 1: one frame interval and one train for the whole set
    files = dataset.files;
    if isempty(files)
        error('hfs_batch_analysis:noFiles', 'The dataset has no loadable files.');
    end
    fi = arrayfun(@(f) f.meta.frame_interval_ms, files);
    if (max(fi) - min(fi)) > config.metadata.interval_tol_frac * median(fi)
        error('hfs_batch_analysis:mixedInterval', ...
            'Files have different Frame_Interval_ms; split them first (run_pipeline does).');
    end
    nP = arrayfun(@(f) f.meta.train_pulses, files);
    hz = arrayfun(@(f) f.meta.train_hz, files);
    if numel(unique(nP)) > 1 || numel(unique(hz)) > 1
        error('hfs_batch_analysis:mixedTrain', ...
            'Files have different trains (%s); analyse each train from its own folder.', ...
            strjoin(unique(arrayfun(@(a, b) sprintf('%gAP-%gHz', a, b), nP, hz, 'UniformOutput', false)), ', '));
    end
    pre_f  = round(config.hfs.epoch_pre_ms / median(fi));
    post_f = round((nP(1) * 1000 / hz(1) + config.hfs.post_train_ms) / median(fi));
    config = apply_recording_timing(config, median(fi), pre_f + post_f + 1);
    config = apply_train_timing(config, pre_f + 1, nP(1), hz(1));
    resolvable = hfs_pulses_resolvable(config);
    fprintf('  %g pulses at %g Hz | frame interval %.4g ms (%.1f frames/pulse) | epoch %d frames, stim1 at %d\n', ...
        nP(1), hz(1), config.frame_dur_ms, config.hfs.isi_ms / config.frame_dur_ms, ...
        config.expected_frames, config.hfs.stim1_frame);
    if ~resolvable
        fprintf(['  Single pulses NOT resolvable (< %d frames/pulse): train-level metrics only ' ...
            '(charge, plateau, peak, decay).\n'], config.hfs.min_frames_per_pulse);
    end

    %% Step 2: process every file
    fprintf('\nProcessing files...\n');
    nFiles = numel(files);
    fileResults = cell(nFiles, 1);
    tic;
    for i = 1:nFiles
        try
            fr = hfs_pipeline(files(i), struct('verbose', false), config);
            fileResults{i} = fr;
            h = fr.hfs;
            fprintf('  [%d/%d] %s (%s %s): %d ROIs, %d pass QC, %d responders (%.0f%%)\n', ...
                i, nFiles, fr.filename, fr.genotype, fr.drug, fr.numROIs, h.n_analyzable, ...
                h.n_responders, 100 * h.responder_fraction);
            if ~isempty(fr.train_check), fprintf('        CHECK: %s\n', fr.train_check); end
            if ~isempty(fr.meta.flags), fprintf('        Step1 flags: %s\n', strjoin(fr.meta.flags, ', ')); end
        catch ME
            fprintf('  [%d/%d] FAILED %s: %s\n', i, nFiles, files(i).name, ME.message);
        end
    end
    process_time = toc;
    ok = ~cellfun(@isempty, fileResults(:)');

    file_checks = table('Size', [0 5], 'VariableTypes', {'string', 'string', 'string', 'double', 'string'}, ...
        'VariableNames', {'filename', 'genotype', 'drug', 'median_peak_latency_ms', 'check'});
    for i = find(ok)
        fr = fileResults{i};
        c = strjoin([{fr.train_check}, fr.meta.flags(:)'], '; ');
        c = regexprep(c, '^(; )+|(; )+$', '');
        file_checks(end+1, :) = {string(fr.filename), string(fr.genotype), string(fr.drug), ...
            fr.median_lat_ms, string(c)}; %#ok<AGROW>
    end

    %% Step 3: groups (genotype x drug) with per-coverslip grand averages
    fprintf('\nAggregating by genotype x drug (per coverslip)...\n');
    genotypes = {files.genotype}; drugs = {files.drug};
    [geno_list, drug_list] = genotype_drug_order(genotypes(ok), drugs(ok), config);
    groups = [];
    for d = 1:numel(drug_list)
        for g = 1:numel(geno_list)
            mask = ok & strcmp(genotypes, geno_list{g}) & strcmp(drugs, drug_list{d});
            if ~any(mask), continue; end
            gd = aggregate_group(fileResults(mask), config);
            gd.name     = sprintf('%s %s', geno_list{g}, drug_list{d});
            gd.genotype = geno_list{g};
            gd.drug     = drug_list{d};
            gd.color    = group_color(geno_list{g}, d, numel(drug_list), config);
            if isempty(groups), groups = gd; else, groups(end+1) = gd; end %#ok<AGROW>
            fprintf('  %s: %d files, %d coverslips, %d QC ROIs, %d responders\n', ...
                gd.name, gd.num_files, numel(gd.cs), gd.n_analyzable, gd.n_responders);
        end
    end

    %% Step 4: figures
    plot_handles = struct();
    if options.createPlots && ~isempty(groups)
        fprintf('\nCreating HFS figures...\n');
        plot_handles.traces  = plot_hfs_traces(groups, config);
        if resolvable
            plot_handles.pulses = plot_hfs_pulses(groups, config);
        end
        plot_handles.summary = plot_hfs_summary(groups, config, resolvable);
        plot_handles.files   = plot_hfs_file_traces(fileResults(ok), groups, config);
    end

    %% Step 5: tables
    results = struct();
    results.dataType            = 'HFS';
    results.groups              = groups;
    results.plot_handles        = plot_handles;
    results.hfs_roi_table       = build_roi_table(fileResults(ok));
    results.hfs_pulse_table     = build_pulse_table(groups, resolvable);
    results.hfs_coverslip_table = build_coverslip_table(groups);
    results.file_checks         = file_checks;
    results.processing_info     = struct('folder_path', dataset.folder, ...
        'metadata_file', dataset.metadata_file, 'process_time', process_time, ...
        'timestamp', datetime('now'), 'frame_interval_ms', config.frame_dur_ms, ...
        'n_pulses', nP(1), 'train_hz', hz(1), 'pulses_resolvable', resolvable, ...
        'epoch_frames', config.expected_frames, 'epoch_stim1_frame', config.hfs.stim1_frame);

    print_summary(groups, resolvable);
    flagged = strlength(file_checks.check) > 0;
    if any(flagged)
        fprintf('\n*** %d file(s) with checks / Step1 flags (see hfs_file_checks.csv) ***\n', sum(flagged));
        disp(file_checks(flagged, :));
    end
end

%% ======================================================================
%% Aggregation
%% ======================================================================
function gd = aggregate_group(frs, config)
    % Coverslip grand averages (QC-passed ROIs of every file on that coverslip)
    % and their metrics. Each region (c#) mean is kept for the overlay.
    gd = struct();
    gd.num_files = numel(frs);
    gd.filenames = cellfun(@(f) f.filename, frs, 'UniformOutput', false);
    gd.n_analyzable = sum(cellfun(@(f) f.hfs.n_analyzable, frs));
    gd.n_responders = sum(cellfun(@(f) f.hfs.n_responders, frs));
    cs_names = unique(cellfun(@(f) f.coverslip, frs, 'UniformOutput', false), 'stable');
    cs = struct('cs', {}, 'trace', {}, 'sem', {}, 'n', {}, 'regions', {}, 'm', {}, ...
        'responder_fraction', {});
    for c = 1:numel(cs_names)
        sel = frs(cellfun(@(f) strcmp(f.coverslip, cs_names{c}), frs));
        D = []; regions = struct('label', {}, 'trace', {}, 'n', {});
        nA = 0; nR = 0;
        for k = 1:numel(sel)
            Dk = double(sel{k}.dfof_data(:, sel{k}.hfs.analyzed));
            D = [D, Dk]; %#ok<AGROW>
            regions(end+1) = struct('label', sel{k}.region, ...
                'trace', mean(Dk, 2, 'omitnan')', 'n', size(Dk, 2)); %#ok<AGROW>
            nA = nA + sel{k}.hfs.n_analyzable; nR = nR + sel{k}.hfs.n_responders;
        end
        n   = sum(isfinite(D), 2);
        tr  = mean(D, 2, 'omitnan');
        sem = std(D, 0, 2, 'omitnan') ./ sqrt(max(n, 1));
        cs(end+1) = struct('cs', cs_names{c}, 'trace', tr', 'sem', sem', 'n', size(D, 2), ...
            'regions', regions, 'm', hfs_train_metrics(tr, [], config, true), ...
            'responder_fraction', nR / max(nA, 1)); %#ok<AGROW>
    end
    gd.cs = cs;
end

function [mn, se, vals] = across_cs(gd, field)
    % Mean +/- SEM across coverslips of a scalar metric (or a per-pulse row).
    if strcmp(field, 'responder_fraction')
        vals = [gd.cs.responder_fraction]';
    else
        vals = cell2mat(arrayfun(@(c) c.m.(field)(:)', gd.cs(:), 'UniformOutput', false));
    end
    mn = mean(vals, 1, 'omitnan');
    nn = sum(isfinite(vals), 1);
    se = std(vals, 0, 1, 'omitnan') ./ sqrt(max(nn, 1));
    se(nn < 2) = NaN;
end

%% ======================================================================
%% Figures
%% ======================================================================
function fig = plot_hfs_traces(groups, config)
    % Group average traces (mean of coverslip grand averages), one column per
    % drug, genotypes overlaid. Rows: whole epoch, first 300 ms of the train,
    % last 300 ms of the train. Thin lines = individual coverslips.
    p = config.hfs; fd = config.frame_dur_ms; s1 = p.stim1_frame;
    drugs = unique({groups.drug}, 'stable'); nD = numel(drugs);
    fig = figure('Name', 'HFS average traces', 'Position', [60 60 max(700, 560*nD) 900], 'Color', 'w');
    tl = tiledlayout(fig, 3, nD, 'TileSpacing', 'compact', 'Padding', 'compact');
    train_ms = p.n_pulses * p.isi_ms;
    xl = {[-p.epoch_pre_ms / 2, train_ms + p.post_train_ms], [-50 300], [train_ms - 300, train_ms + 250]};
    rowname = {'whole train', 'train start', 'train end'};
    axs = gobjects(1, 0);
    for r = 1:3
        for d = 1:nD
            ax = nexttile(tl, (r - 1) * nD + d); hold(ax, 'on'); axs(end+1) = ax; %#ok<AGROW>
            setappdata(ax, 'ymin', config.visualization.ymin_mean_trace);
            gsel = find(strcmp({groups.drug}, drugs{d}));
            for g = gsel
                gd = groups(g);
                T = vertcat(gd.cs.trace); t = ((1:size(T, 2)) - s1) * fd;
                if size(T, 1) > 1
                    plot(ax, t, T', '-', 'Color', [gd.color 0.35], 'LineWidth', 0.5, 'HandleVisibility', 'off');
                end
                plot(ax, t, mean(T, 1, 'omitnan'), '-', 'Color', gd.color, 'LineWidth', 1.2, ...
                    'DisplayName', sprintf('%s (%d cs)', gd.genotype, size(T, 1)));
            end
            if r > 1
                pt = (p.pulse_frames - s1) * fd;
                pt = pt(pt >= xl{r}(1) & pt <= xl{r}(2));
                if ~isempty(pt), xline(ax, pt, ':', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off'); end
            else
                xregion(ax, 0, train_ms, 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.08, 'HandleVisibility', 'off');
            end
            xlim(ax, xl{r}); grid(ax, 'on'); box(ax, 'on');
            xlabel(ax, 'ms from first pulse');
            if d == 1, ylabel(ax, sprintf('\\DeltaF/F (%s)', rowname{r}), 'Interpreter', 'tex'); end
            if r == 1, title(ax, drugs{d}, 'Interpreter', 'none'); legend(ax, 'Location', 'northeast'); end
        end
    end
    share_ylim(axs, config);   % one y-scale for all panels (largest trace)
    title(tl, {sprintf('HFS %g AP @ %g Hz: coverslip grand-average dF/F (all QC ROIs, no gating)', ...
        p.n_pulses, p.train_hz), 'Thick = group mean, thin = coverslips; dotted = pulses'}, ...
        'FontWeight', 'bold');
end

function fig = plot_hfs_pulses(groups, config)
    % Per-pulse curves from the coverslip averages (mean +/- SEM across coverslips).
    p = config.hfs; n = 1:p.n_pulses;
    fig = figure('Name', 'HFS per-pulse', 'Position', [80 80 1400 900], 'Color', 'w');
    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    spec = {'amp', 'Phasic amplitude A_n (\DeltaF/F above pre-pulse level)'; ...
            'norm_amp', 'A_n / A_1'; ...
            'tonic', 'Tonic level before pulse n (\DeltaF/F above baseline)'; ...
            'cum_amp', 'Cumulative phasic amplitude \Sigma A (\DeltaF/F)'};
    for s = 1:size(spec, 1)
        ax = nexttile(tl); hold(ax, 'on');
        for g = 1:numel(groups)
            gd = groups(g);
            [mn, se] = across_cs(gd, spec{s, 1});
            ok = isfinite(se);
            if any(ok)
                fill(ax, [n(ok), fliplr(n(ok))], [mn(ok) + se(ok), fliplr(mn(ok) - se(ok))], gd.color, ...
                    'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            end
            plot(ax, n, mn, '.-', 'Color', gd.color, 'LineWidth', 1, 'MarkerSize', 8, 'DisplayName', gd.name);
            if strcmp(spec{s, 1}, 'cum_amp')
                pool = mean(arrayfun(@(c) c.m.pool_est, gd.cs), 'omitnan');
                repl = mean(arrayfun(@(c) c.m.repl_per_pulse, gd.cs), 'omitnan');
                plot(ax, [0 n(end)], pool + repl * [0 n(end)], '--', 'Color', gd.color, 'HandleVisibility', 'off');
            end
        end
        if strcmp(spec{s, 1}, 'norm_amp'), yline(ax, 1, 'k:', 'HandleVisibility', 'off'); end
        if strcmp(spec{s, 1}, 'cum_amp'), xlim(ax, [0 n(end)]); end
        grid(ax, 'on'); box(ax, 'on'); xlabel(ax, 'Pulse number');
        ylabel(ax, spec{s, 2}, 'Interpreter', 'tex');
        if s == 1, legend(ax, 'Location', 'best', 'Interpreter', 'none'); end
    end
    title(tl, {sprintf('HFS %g AP @ %g Hz per pulse (coverslip average traces; mean \\pm SEM across coverslips)', ...
        p.n_pulses, p.train_hz), sprintf(['Dashed: line through the last %d cumulative points, back-extrapolated ' ...
        'to pulse 0 (intercept ~ releasable pool, slope ~ replenishment per pulse)'], p.cumulative_fit_pulses)}, ...
        'FontWeight', 'bold', 'Interpreter', 'tex');
end

function fig = plot_hfs_summary(groups, config, resolvable)
    % Train-level metrics: one dot per coverslip, bar = mean across coverslips.
    spec = {'q_train', 'Train charge (\DeltaF/F\cdotms)'; ...
            'plateau', 'Plateau (\DeltaF/F, last 500 ms)'; ...
            'peak_train', 'Peak during train (\DeltaF/F)'; ...
            'decay_tau_ms', 'Post-train decay \tau (ms)'; ...
            't_half_ms', 'Post-train half-decay (ms)'; ...
            'responder_fraction', 'Responder fraction'};
    if resolvable
        spec = [spec; {'a1', 'A_1 (\DeltaF/F)'; 'ss_amp', 'Steady-state A_n (\DeltaF/F)'; ...
            'ss_ratio', 'Steady-state A_n / A_1'; 'ss_tonic', 'Steady-state tonic (\DeltaF/F)'; ...
            'pool_est', 'Pool estimate (\DeltaF/F, intercept)'; 'repl_per_pulse', 'Replenishment (\DeltaF/F per pulse)'}];
    end
    spec{2, 2} = sprintf('Plateau (\\DeltaF/F, last %g ms)', config.hfs.plateau_window_ms);
    nS = size(spec, 1); nG = numel(groups);
    nc = min(nS, 6); nr = ceil(nS / nc);
    fig = figure('Name', 'HFS summary', 'Position', [100 100 300*nc 380*nr], 'Color', 'w');
    tl = tiledlayout(fig, nr, nc, 'TileSpacing', 'compact', 'Padding', 'compact');
    for s = 1:nS
        ax = nexttile(tl); hold(ax, 'on');
        lo = 0;
        for g = 1:nG
            [mn, ~, vals] = across_cs(groups(g), spec{s, 1});
            lo = min([lo; vals(isfinite(vals))]);
            setappdata(ax, 'ymin', lo);   % keep negative values visible (apply_ymin_zero)
            bar(ax, g, mn, 0.6, 'FaceColor', groups(g).color, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
            v = vals(isfinite(vals));
            scatter(ax, g + linspace(-0.08, 0.08, numel(v)) * (numel(v) > 1), v, 18, 'k', 'filled');
        end
        set(ax, 'XTick', 1:nG, 'XTickLabel', {groups.name}, 'XTickLabelRotation', 35, 'TickLabelInterpreter', 'none');
        title(ax, spec{s, 2}, 'Interpreter', 'tex', 'FontSize', 9);
        grid(ax, 'on'); box(ax, 'on');
    end
    title(tl, sprintf(['HFS %g AP @ %g Hz train metrics (bar = mean across coverslips, dot = coverslip; ' ...
        'measured on coverslip average traces)'], config.hfs.n_pulses, config.hfs.train_hz), 'FontWeight', 'bold');
end

function fig = plot_hfs_file_traces(frs, groups, config)
    % QC: each file's mean trace (QC-passed ROIs), with any check in the title.
    p = config.hfs; fd = config.frame_dur_ms; s1 = p.stim1_frame;
    nF = numel(frs); nc = min(nF, 2); nr = ceil(nF / nc);
    fig = figure('Name', 'HFS file traces', 'Position', [120 40 1500 max(300, 200*nr)], 'Color', 'w');
    tl = tiledlayout(fig, nr, nc, 'TileSpacing', 'compact', 'Padding', 'compact');
    axs = gobjects(1, 0);
    for i = 1:nF
        fr = frs{i};
        g = find(strcmp({groups.genotype}, fr.genotype) & strcmp({groups.drug}, fr.drug), 1);
        ax = nexttile(tl); hold(ax, 'on'); axs(end+1) = ax; %#ok<AGROW>
        setappdata(ax, 'ymin', config.visualization.ymin_mean_trace);
        tr = mean(fr.dfof_data(:, fr.hfs.analyzed), 2, 'omitnan');
        t = ((1:numel(tr)) - s1) * fd;
        plot(ax, t, tr, '-', 'Color', groups(g).color, 'LineWidth', 0.8);
        xline(ax, [0, p.n_pulses * p.isi_ms], 'r-');
        xlim(ax, [-500, p.n_pulses * p.isi_ms + p.post_train_ms]);
        ttl = sprintf('%s | %d QC ROIs, %d resp', fr.filename, fr.hfs.n_analyzable, fr.hfs.n_responders);
        if ~isempty(fr.train_check) || ~isempty(fr.meta.flags)
            ttl = {ttl, ['CHECK: ' strjoin([{fr.train_check}, fr.meta.flags(:)'], ' ')]};
        end
        title(ax, ttl, 'Interpreter', 'none', 'FontSize', 7);
        grid(ax, 'on'); box(ax, 'on');
    end
    share_ylim(axs, config);   % one y-scale for all panels (largest trace)
    title(tl, 'HFS per-file mean dF/F (QC-passed ROIs); red = train start / end', 'FontWeight', 'bold');
end

%% ======================================================================
%% Tables / summary
%% ======================================================================
function T = build_roi_table(frs)
    rows = {};
    for i = 1:numel(frs)
        fr = frs{i}; h = fr.hfs; n = fr.numROIs;
        rows{end+1} = table(repmat(string(fr.genotype), n, 1), repmat(string(fr.drug), n, 1), ...
            repmat(string(fr.filename), n, 1), repmat(string(fr.coverslip), n, 1), (1:n)', ...
            h.analyzed(:), h.noise_excluded(:), h.is_responder(:), h.prestim_sd(:), ...
            h.q_train(:), h.plateau(:), h.peak_train(:), h.a1(:), h.ss_amp(:), h.ss_ratio(:), ...
            h.ss_tonic(:), h.pool_est(:), h.repl_per_pulse(:), h.median_lat_ms(:), ...
            'VariableNames', {'genotype', 'drug', 'filename', 'coverslip', 'roi', 'qc_pass', ...
            'noise_excluded', 'is_responder', 'prestim_sd', 'q_train', 'plateau', 'peak_train', ...
            'a1', 'ss_amp', 'ss_ratio', 'ss_tonic', 'pool_est', 'repl_per_pulse', 'median_peak_latency_ms'}); %#ok<AGROW>
    end
    T = vertcat(rows{:});
end

function T = build_pulse_table(groups, resolvable)
    rows = {};
    if ~resolvable, T = table(); return; end
    for g = 1:numel(groups)
        for c = 1:numel(groups(g).cs)
            e = groups(g).cs(c); n = numel(e.m.amp);
            rows{end+1} = table(repmat(string(groups(g).genotype), n, 1), repmat(string(groups(g).drug), n, 1), ...
                repmat(string(e.cs), n, 1), (1:n)', e.m.amp(:), e.m.norm_amp(:), e.m.tonic(:), ...
                e.m.cum_amp(:), e.m.lat_ms(:), ...
                'VariableNames', {'genotype', 'drug', 'coverslip', 'pulse', 'phasic_amp', ...
                'phasic_amp_norm', 'tonic', 'cumulative_amp', 'peak_latency_ms'}); %#ok<AGROW>
        end
    end
    T = vertcat(rows{:});
end

function T = build_coverslip_table(groups)
    rows = {};
    for g = 1:numel(groups)
        for c = 1:numel(groups(g).cs)
            e = groups(g).cs(c); m = e.m;
            rows{end+1} = table(string(groups(g).genotype), string(groups(g).drug), string(e.cs), ...
                numel(e.regions), e.n, e.responder_fraction, m.q_train, m.plateau, m.peak_train, ...
                m.a1, m.ss_amp, m.ss_ratio, m.ss_tonic, m.pool_est, m.repl_per_pulse, ...
                m.decay_tau_ms, m.decay_r2, m.t_half_ms, m.median_lat_ms, ...
                'VariableNames', {'genotype', 'drug', 'coverslip', 'n_regions', 'n_qc_rois', ...
                'responder_fraction', 'q_train', 'plateau', 'peak_train', 'a1', 'ss_amp', 'ss_ratio', ...
                'ss_tonic', 'pool_est', 'repl_per_pulse', 'decay_tau_ms', 'decay_r2', 't_half_ms', ...
                'median_peak_latency_ms'}); %#ok<AGROW>
        end
    end
    T = vertcat(rows{:});
end

function print_summary(groups, resolvable)
    fprintf('\n=== HFS SUMMARY (coverslip average traces; mean across coverslips) ===\n');
    fprintf('(Statistical unit = coverslip; regions / ROIs within it are pseudoreplicates.)\n');
    for g = 1:numel(groups)
        gd = groups(g);
        q = across_cs(gd, 'q_train'); pl = across_cs(gd, 'plateau'); tau = across_cs(gd, 'decay_tau_ms');
        fprintf('%-22s %d cs | Q %.1f | plateau %.4f | decay tau %.0f ms', gd.name, numel(gd.cs), q, pl, tau);
        if resolvable
            fprintf(' | A1 %.4f | SS %.4f (x%.2f of A1) | pool %.3f | repl %.4f/pulse', ...
                across_cs(gd, 'a1'), across_cs(gd, 'ss_amp'), across_cs(gd, 'ss_ratio'), ...
                across_cs(gd, 'pool_est'), across_cs(gd, 'repl_per_pulse'));
        end
        fprintf('\n');
    end
    fprintf('=====================================================================\n');
end
