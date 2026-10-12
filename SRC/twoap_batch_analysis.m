function results = twoap_batch_analysis(dataset, options)
    % TWOAP_BATCH_ANALYSIS - Batch analysis of 2AP (paired-pulse / PPF) data
    %
    % Counterpart to oneap_batch_analysis for paired-pulse recordings. Takes a
    % dataset from load_dataset_folder: each file is one region at one ISI, with
    % stim1, the ISI and the frame interval from Recordings.csv. All files must
    % share one frame interval (run_pipeline splits mixed folders first). Files
    % are grouped by genotype x drug, processed through twoap_pipeline, and the
    % paired-pulse ratio (PPR = A2/A1) is aggregated PER COVERSLIP at each ISI, for
    % both the peak-max and time-locked methods. Produces the PPR-vs-ISI figure, a
    % method-comparison montage, a per-ROI metrics table and a summary. Files
    % whose stims don't match the reported ISI are listed (isi_check).
    %
    % Inputs:
    %   dataset - struct from load_dataset_folder (one frame interval)
    %   options - struct: verbose, createPlots, createMontage
    %
    % Output:
    %   results - struct with .groups, .ppr_by_isi, .plot_handles,
    %             .twoap_metrics_table, .isi_checks, .processing_info

    if nargin < 2, options = struct(); end
    if ~isfield(options, 'verbose'),      options.verbose = true; end
    if ~isfield(options, 'createPlots'),  options.createPlots = true; end
    if ~isfield(options, 'createMontage'),options.createMontage = true; end
    % saveFigures defaults OFF: run_pipeline saves all open figures after
    % this returns, so leave it off there to avoid double-saving. Set true when
    % calling twoap_batch_analysis directly to get a dated figures folder.
    if ~isfield(options, 'saveFigures'),  options.saveFigures = false; end

    if ~isstruct(dataset) || ~isfield(dataset, 'files')
        error('twoap_batch_analysis:input', ...
            'twoap_batch_analysis takes a dataset from load_dataset_folder.');
    end
    folder_path = dataset.folder;
    config = local_config();

    fprintf('=== 2AP Paired-Pulse (PPF) Batch Analysis ===\n');
    fprintf('Folder: %s\n', folder_path);

    %% Step 1: files + one shared time base (stim1 at the same epoch index)
    files = dataset.files;
    if isempty(files)
        error('twoap_batch_analysis:noFiles', 'The dataset has no loadable files.');
    end
    fi = arrayfun(@(f) f.meta.frame_interval_ms, files);
    if (max(fi) - min(fi)) > config.metadata.interval_tol_frac * median(fi)
        error('twoap_batch_analysis:mixedInterval', ...
            'Files have different Frame_Interval_ms (%s); split them first (run_pipeline does).', ...
            strjoin(unique(cellstr(num2str(fi(:), '%.4g'))), ', '));
    end
    pre_f  = round(config.twoAP.epoch_pre_ms  / median(fi));
    post_f = round(config.twoAP.epoch_post_ms / median(fi));
    config = apply_recording_timing(config, median(fi), pre_f + post_f + 1);
    config.twoAP.stim1_frame = pre_f + 1;            % stim1 index within each epoch
    fprintf('  Frame interval %.4g ms (metadata) | epoch %d frames, stim1 at %d\n', ...
        config.frame_dur_ms, config.expected_frames, config.twoAP.stim1_frame);

    %% Step 2: process every file
    fprintf('\nProcessing files...\n');
    nFiles = numel(files);
    fileResults = cell(nFiles, 1);
    tic;
    for i = 1:nFiles
        try
            fileResults{i} = twoap_pipeline(files(i), struct('verbose', false), config);
            fr = fileResults{i};
            t  = fr.twoap;
            fprintf('  [%d/%d] %s (%s %s, ISI %g ms, %s): %d ROIs, %d pulse-1 responders, PPR n=%d/%d/%d (max/sync/p1lock)\n', ...
                i, nFiles, fr.filename, fr.genotype, fr.drug, fr.isi_ms, fr.ppf_mode, fr.numROIs, ...
                t.n_responders1, t.n_ppr_valid_max, t.n_ppr_valid_sync, t.n_ppr_valid_p1lock);
            if ~isempty(fr.isi_check)
                fprintf('        ISI CHECK: %s\n', fr.isi_check);
            end
            if ~isempty(fr.meta.flags)
                fprintf('        Step1 flags: %s\n', strjoin(fr.meta.flags, ', '));
            end
        catch ME
            fprintf('  [%d/%d] FAILED %s: %s\n', i, nFiles, files(i).name, ME.message);
            fileResults{i} = [];
        end
    end
    process_time = toc;
    fprintf('  Processing complete in %.2f s\n', process_time);

    % Files whose stims don't match the reported ISI (look at these when
    % interpreting PPR; the ISI-derived stim2 was used for them).
    ok = ~cellfun(@isempty, fileResults(:)');
    isi_checks = table('Size', [0 4], 'VariableTypes', {'string', 'double', 'string', 'string'}, ...
        'VariableNames', {'filename', 'isi_ms', 'ppf_mode', 'isi_check'});
    for i = find(ok)
        fr = fileResults{i};
        if ~isempty(fr.isi_check)
            isi_checks(end+1, :) = {string(fr.filename), fr.isi_ms, string(fr.ppf_mode), string(fr.isi_check)}; %#ok<AGROW>
        end
    end

    %% Step 3: aggregate per group (genotype x drug)
    fprintf('\nAggregating by genotype x drug / ISI (per coverslip)...\n');
    genotypes = {files.genotype};
    drugs     = {files.drug};
    [geno_list, drug_list] = genotype_drug_order(genotypes(ok), drugs(ok), config);
    groups = [];
    for d = 1:numel(drug_list)
        for g = 1:numel(geno_list)
            mask = ok & strcmp(genotypes, geno_list{g}) & strcmp(drugs, drug_list{d});
            if ~any(mask), continue; end
            gname = sprintf('%s %s', geno_list{g}, drug_list{d});
            gd = aggregate_group_2ap(fileResults(mask), gname, ...
                group_color(geno_list{g}, d, numel(drug_list), config));
            gd.genotype = geno_list{g};
            gd.drug     = drug_list{d};
            if isempty(groups), groups = gd; else, groups(end+1) = gd; end %#ok<AGROW>
            fprintf('  %s: %d files, %d coverslips, %d ISIs (%s ms), %d ROIs, %d pulse-1 responders\n', ...
                gd.name, gd.num_files, gd.num_coverslips, numel(gd.isi_list), ...
                strjoin(arrayfun(@(x) num2str(x), gd.isi_list, 'UniformOutput', false), ','), ...
                gd.total_rois, gd.n_responders1);
        end
    end

    %% Step 4: PPR summaries (per coverslip, by ISI)
    ppr_by_isi = struct();
    for g = 1:numel(groups)
        ppr_by_isi.(matlab.lang.makeValidName(groups(g).name)) = struct( ...
            'ratio_of_means', compute_rom_by_isi(groups(g), config), ...  % PRIMARY (unbiased)
            'peak_max',   compute_ppr_by_isi(groups(g), 'ppr_max'), ...   % per-ROI diagnostic
            'p1lock',     compute_ppr_by_isi(groups(g), 'ppr_p1lock'), ...
            'sync',       compute_ppr_by_isi(groups(g), 'ppr_sync'));
    end

    %% Step 5: plots
    plot_handles = struct();
    if options.createPlots
        fprintf('\nCreating 2AP figures...\n');
        plot_handles.ppr_vs_isi   = plot_twoap_ppr(groups, config);
        plot_handles.amplitudes   = plot_twoap_amplitudes(groups, config);
        plot_handles.gating_bias  = plot_twoap_gating_bias(groups, config);
        plot_handles.cs_average   = plot_twoap_cs_average(groups, config);
        if options.createMontage
            plot_handles.montage = plot_twoap_montage(groups, config);
        end
        fprintf('  Figures created.\n');
    end

    %% Step 5b: optionally save the figures this run created (opt-in)
    output_dir = '';
    if options.createPlots && options.saveFigures
        output_dir = make_figures_dir(folder_path);
        figs = gobjects(1, 0);
        fn = fieldnames(plot_handles);
        for i = 1:numel(fn)
            h = plot_handles.(fn{i});
            figs = [figs, reshape(h(isgraphics(h, 'figure')), 1, [])]; %#ok<AGROW>
        end
        n_saved = save_open_figures(output_dir, figs);
        fprintf('  Saved %d figures to:\n    %s\n', n_saved, output_dir);
    end

    %% Step 6: per-ROI metrics table
    twoap_metrics_table = build_twoap_table(groups);

    %% Assemble results
    results = struct();
    results.dataType            = '2AP';
    results.groups              = groups;
    results.ppr_by_isi          = ppr_by_isi;
    results.plot_handles        = plot_handles;
    results.twoap_metrics_table = twoap_metrics_table;
    results.isi_checks          = isi_checks;
    results.output_dir          = output_dir;   % '' unless options.saveFigures
    results.processing_info     = struct('folder_path', folder_path, ...
        'metadata_file', dataset.metadata_file, ...
        'process_time', process_time, 'timestamp', datetime('now'), ...
        'frame_interval_ms', config.frame_dur_ms, 'epoch_frames', config.expected_frames, ...
        'epoch_stim1_frame', config.twoAP.stim1_frame);

    print_summary_2ap(groups, config);
    print_gating_bias(groups, config);
    if height(isi_checks) > 0
        fprintf('\n*** %d file(s) whose stims do not match the reported ISI (stim2 set from ISI) ***\n', ...
            height(isi_checks));
        disp(isi_checks);
    else
        fprintf('\nISI check: all files'' stims match the reported ISI.\n');
    end
end

%% ======================================================================
%% Aggregation
%% ======================================================================
function gd = aggregate_group_2ap(fileResults, name, color)
    gd = struct();
    gd.name = name; gd.color = color;
    gd.num_files = numel(fileResults);

    ppr_max = []; ppr_sync = []; ppr_p1lock = []; peak1_amp = []; peak2_amp = [];
    a1_sync = []; a2_sync = []; a2_p1lock = [];
    q_pair = [];
    is_responder1 = []; is_responder2 = []; prestim_sd = [];
    isi_of_roi = []; coverslip_of_roi = {}; file_of_roi = []; filenames = {};
    file_isi_check = {}; file_ppf_mode = {}; roi_index = [];

    % Reference file for the method-comparison montage: prefer the SHORTEST ISI
    % with enough pulse-1 responders (methods diverge most there), else the file
    % with the most responders.
    montage_ref = struct('dfof', [], 'twoap', [], 'isi', NaN, 'filename', '');
    best_short_isi = Inf; best_resp = -1;

    % Coverslip x ISI grand-average accumulator (running nanmean over QC-passed ROIs
    % of all files sharing a coverslip and ISI). Unbiased: peak2/peak1 of this mean
    % trace = ratio-of-means PPR. Also accumulates responder-subset sums for the
    % gating diagnostic. Keys 'coverslip|isi'.
    csMap = containers.Map('KeyType', 'char', 'ValueType', 'any');

    for i = 1:numel(fileResults)
        fr = fileResults{i};
        if isempty(fr), continue; end
        t = fr.twoap;
        filenames{end+1} = fr.filename; %#ok<AGROW>
        file_isi_check{end+1} = fr.isi_check; %#ok<AGROW>
        file_ppf_mode{end+1}  = fr.ppf_mode;  %#ok<AGROW>
        nR = fr.numROIs;

        ppr_max       = [ppr_max,       t.ppr_max];        %#ok<AGROW>
        ppr_sync      = [ppr_sync,      t.ppr_sync];       %#ok<AGROW>
        ppr_p1lock    = [ppr_p1lock,    t.ppr_p1lock];     %#ok<AGROW>
        peak1_amp     = [peak1_amp,     t.peak1_amp];      %#ok<AGROW>
        peak2_amp     = [peak2_amp,     t.peak2_amp];      %#ok<AGROW>
        a1_sync       = [a1_sync,       t.a1_sync];        %#ok<AGROW>
        a2_sync       = [a2_sync,       t.a2_sync];        %#ok<AGROW>
        a2_p1lock     = [a2_p1lock,     t.a2_p1lock];      %#ok<AGROW>
        q_pair        = [q_pair,        t.q_pair];         %#ok<AGROW>
        is_responder1 = [is_responder1, t.is_responder1];  %#ok<AGROW>
        is_responder2 = [is_responder2, t.is_responder2];  %#ok<AGROW>
        prestim_sd    = [prestim_sd,    t.prestim_sd];     %#ok<AGROW>
        isi_of_roi    = [isi_of_roi,    repmat(fr.isi_ms, 1, nR)]; %#ok<AGROW>
        coverslip_of_roi = [coverslip_of_roi, repmat({parse_coverslip_id(fr.filename)}, 1, nR)]; %#ok<AGROW>
        file_of_roi   = [file_of_roi,   repmat(i, 1, nR)]; %#ok<AGROW>
        roi_index     = [roi_index,     1:nR];             %#ok<AGROW>

        % Accumulate coverslip x ISI average traces, nan-aware, over QC-passed ROIs
        % (analyzable + not noise-excluded; NOT filtered on response, so no selection
        % bias). Also keep the same sums restricted to responder subsets, so the
        % gating-bias diagnostic is built on the SAME average-trace estimator as the
        % primary ROM: 'all' -> no-gate (== primary ROM); gate_A1 -> P1 responders;
        % gate_A2 -> P2 responders (peak2 >= upper_sigma*sd); either / both.
        D     = fr.dfof_data;                        % [nF x nROI]
        qc    = isfinite(t.upper_threshold);         % passed noise QC (analyzable, not noise-excluded)
        thr   = t.upper_threshold - t.prestim_base;  % = upper_sigma*sd per ROI (pulse-2 responder cut)
        r1raw = t.is_responder1;
        r2raw = isfinite(t.peak2_amp) & (t.peak2_amp >= thr);
        cond.all         = qc;
        cond.gate_A1     = qc & r1raw;
        cond.gate_A2     = qc & r2raw;
        cond.gate_either = qc & (r1raw | r2raw);
        cond.gate_both   = qc & (r1raw & r2raw);
        cnames = fieldnames(cond);
        key = sprintf('%s|%g', parse_coverslip_id(fr.filename), fr.isi_ms);
        sq  = sum(D(:, qc).^2, 2, 'omitnan');        % [nF x 1] sum of squares over QC ROIs (for SEM)
        % This file's own region (c#) QC-mean numerator/denominator, so each
        % region can later be overlaid on the coverslip grand-average trace.
        rlabel   = parse_region_id(fr.filename);
        [rs, rc] = colsum(D, qc);
        if isKey(csMap, key)
            e = csMap(key);
            m = min(numel(e.sq), numel(sq));         % guard uneven frame counts
            e.sq(1:m) = e.sq(1:m) + sq(1:m);
            for q = 1:numel(cnames)
                [sc, cc] = colsum(D, cond.(cnames{q}));
                e.s.(cnames{q})(1:m) = e.s.(cnames{q})(1:m) + sc(1:m);
                e.c.(cnames{q})(1:m) = e.c.(cnames{q})(1:m) + cc(1:m);
            end
            % Accumulate the per-region sum/count (a region can span >1 file).
            ri = find(strcmp({e.regs.label}, rlabel), 1);
            if isempty(ri)
                e.regs(end+1) = struct('label', rlabel, 's', rs, 'c', rc);
            else
                mr = min(numel(e.regs(ri).s), numel(rs));
                e.regs(ri).s(1:mr) = e.regs(ri).s(1:mr) + rs(1:mr);
                e.regs(ri).c(1:mr) = e.regs(ri).c(1:mr) + rc(1:mr);
            end
            csMap(key) = e;
        else
            sS = struct(); cS = struct();
            for q = 1:numel(cnames)
                [sc, cc] = colsum(D, cond.(cnames{q}));
                sS.(cnames{q}) = sc; cS.(cnames{q}) = cc;
            end
            csMap(key) = struct('cs', parse_coverslip_id(fr.filename), 'isi', fr.isi_ms, ...
                'stim2', t.stim2_frame, 'sq', sq, 's', sS, 'c', cS, ...
                'regs', struct('label', rlabel, 's', rs, 'c', rc));
        end

        % Montage reference selection
        pick = false;
        if t.n_responders1 >= 5 && fr.isi_ms < best_short_isi
            pick = true; best_short_isi = fr.isi_ms;
        elseif isinf(best_short_isi) && t.n_responders1 > best_resp
            pick = true; best_resp = t.n_responders1;
        end
        if pick
            montage_ref.dfof = fr.dfof_data; montage_ref.twoap = t;
            montage_ref.isi = fr.isi_ms; montage_ref.filename = fr.filename;
        end
    end

    gd.ppr_max          = ppr_max;
    gd.ppr_sync         = ppr_sync;
    gd.ppr_p1lock       = ppr_p1lock;
    gd.peak1_amp        = peak1_amp;
    gd.peak2_amp        = peak2_amp;
    gd.a1_sync          = a1_sync;
    gd.a2_sync          = a2_sync;
    gd.a2_p1lock        = a2_p1lock;
    gd.q_pair           = q_pair;
    gd.is_responder1    = logical(is_responder1);
    gd.is_responder2    = logical(is_responder2);
    gd.prestim_sd       = prestim_sd;
    gd.isi_of_roi       = isi_of_roi;
    gd.coverslip_of_roi = coverslip_of_roi;
    gd.file_of_roi      = file_of_roi;
    gd.filenames        = filenames;
    gd.file_isi_check   = file_isi_check;
    gd.roi_index        = roi_index;              % ROI number within its file
    gd.file_ppf_mode    = file_ppf_mode;
    gd.montage_ref      = montage_ref;

    % Finalize coverslip x ISI average traces. cs_avg = the QC-passed grand average
    % (mean +/- SEM), used by the primary ROM + amplitude figures. cs_avg_gate holds
    % the same average restricted to each responder subset, for the gating diagnostic.
    cs_avg      = struct('cs', {}, 'isi', {}, 'trace', {}, 'sem', {}, 'n', {}, 'stim2', {}, 'regions', {});
    cs_avg_gate = struct('cs', {}, 'isi', {}, 'stim2', {}, 'traces', {}, 'n', {});
    ks = csMap.keys;
    for kk = 1:numel(ks)
        e   = csMap(ks{kk});
        ca  = e.c.all; na = max(ca, 1);
        mn  = (e.s.all ./ na)';                     % 1 x nF mean over QC ROIs
        v   = e.sq ./ na - (e.s.all ./ na).^2; v(v < 0) = 0;
        sem = (sqrt(v) ./ sqrt(na))';               % 1 x nF SEM across QC ROIs
        % Per-region (c#) QC-mean traces for the within-coverslip overlay.
        regions = struct('label', {}, 'trace', {}, 'n', {});
        for rr = 1:numel(e.regs)
            rn = max(e.regs(rr).c, 1);
            regions(end+1) = struct('label', e.regs(rr).label, ...      %#ok<AGROW>
                'trace', (e.regs(rr).s ./ rn)', 'n', round(median(e.regs(rr).c)));
        end
        cs_avg(end+1) = struct('cs', e.cs, 'isi', e.isi, 'trace', mn, ...    %#ok<AGROW>
            'sem', sem, 'n', round(median(ca)), 'stim2', e.stim2, 'regions', regions);
        cnames = fieldnames(e.s);
        traces = struct(); ncnt = struct();
        for q = 1:numel(cnames)
            cc = e.c.(cnames{q});
            traces.(cnames{q}) = (e.s.(cnames{q}) ./ max(cc, 1))';   % 1 x nF subset mean
            ncnt.(cnames{q})   = round(median(cc));
        end
        cs_avg_gate(end+1) = struct('cs', e.cs, 'isi', e.isi, 'stim2', e.stim2, ...  %#ok<AGROW>
            'traces', traces, 'n', ncnt);
    end
    gd.cs_avg           = cs_avg;
    gd.cs_avg_gate      = cs_avg_gate;

    gd.total_rois       = numel(ppr_max);
    gd.n_responders1    = sum(gd.is_responder1);
    gd.num_coverslips   = numel(unique(coverslip_of_roi));
    gd.isi_list         = unique(isi_of_roi);
end

function [s, c] = colsum(D, mask)
    % Column-subset running sum + finite-count for a [frames x ROIs] block.
    if any(mask)
        Dm = D(:, mask);
        s = sum(Dm, 2, 'omitnan'); c = sum(isfinite(Dm), 2);
    else
        s = zeros(size(D, 1), 1); c = zeros(size(D, 1), 1);
    end
end

function S = compute_ppr_by_isi(gd, field)
    % Per-coverslip mean PPR at each ISI (responder-gated: PPR is already NaN for
    % non-responders). Returns sorted ISIs, the per-coverslip value vector at
    % each, and the mean +/- SEM across coverslips.
    ppr = gd.(field);
    isi = gd.isi_of_roi;
    cs  = gd.coverslip_of_roi;
    keep = isfinite(ppr);

    uISI = unique(isi(keep));
    S = struct('isi', uISI, 'cs_vals', {cell(1, numel(uISI))}, ...
        'mean', nan(1, numel(uISI)), 'sem', nan(1, numel(uISI)), 'n_cs', zeros(1, numel(uISI)));
    for k = 1:numel(uISI)
        m   = keep & (isi == uISI(k));
        csk = cs(m); pprk = ppr(m);
        ucs = unique(csk, 'stable');
        vals = nan(1, numel(ucs));
        for j = 1:numel(ucs)
            vals(j) = mean(pprk(strcmp(csk, ucs{j})));
        end
        vals = vals(isfinite(vals));
        S.cs_vals{k} = vals;
        S.n_cs(k)    = numel(vals);
        if ~isempty(vals)
            S.mean(k) = mean(vals);
            S.sem(k)  = std(vals) / sqrt(numel(vals));
        end
    end
end

%% ======================================================================
%% Plots
%% ======================================================================
function fig = plot_twoap_ppr(groups, config)
    % PPR vs ISI. The PRIMARY estimator is the unbiased ratio-of-means (peak2/peak1
    % and AUC2/AUC1 of each coverslip's grand-average trace); the per-ROI
    % responder-gated mean is shown alongside as a DIAGNOSTIC because gating on the
    % pulse-1 responders inflates the A1 denominator (regression to the mean),
    % biasing PPR below 1 even when the population facilitates.
    fig = figure('Name', '2AP PPR vs ISI', 'Position', [90 120 1600 520]);

    % Ratio-of-means per group (one value per coverslip per ISI)
    rom = arrayfun(@(gd) compute_rom_by_isi(gd, config), groups, 'UniformOutput', false);

    % --- Panel 1: ratio-of-means, peak PPR (PRIMARY) ---
    ax1 = subplot(1, 3, 1); hold on;
    for g = 1:numel(groups)
        draw_ppr_series(ax1, rom{g}.isi, rom{g}.ppr_pk, groups(g).color, groups(g).name, 'o');
    end
    finish_ppr_axis(ax1, 'PPR (A2 / A1)');
    title(ax1, {'\bfRatio-of-means \rm(primary)', 'peak2 / peak1 of avg trace \bullet all ROIs'}, ...
        'Interpreter', 'tex');

    % --- Panel 2: total charge over the encompassing window (absolute, PRIMARY) ---
    % Charge in dF/F*ms, its own y-scale (NOT a ratio, so no PPR=1 line and not
    % shared with the ratio panels). Grows with ISI by construction.
    ax2 = subplot(1, 3, 2); hold on;
    for g = 1:numel(groups)
        draw_ppr_series(ax2, rom{g}.isi, rom{g}.q_pair, groups(g).color, groups(g).name, 'o');
    end
    set(ax2, 'XScale', 'log');
    xlabel(ax2, 'Inter-stimulus interval (ms)');
    ylabel(ax2, 'Total charge Q (dF/F\cdotms)');
    grid(ax2, 'on'); box(ax2, 'on'); legend(ax2, 'Location', 'best');
    title(ax2, {'\bfTotal charge \rm(primary)', ...
        sprintf('\\intdF/F over [stim1, stim2+%gms] of avg trace \\bullet all ROIs', ...
        config.twoAP.auc_post_ms)}, 'Interpreter', 'tex');

    % --- Panel 3: per-ROI responder-gated peak-max PPR (DIAGNOSTIC) ---
    ax3 = subplot(1, 3, 3); hold on;
    for g = 1:numel(groups)
        S = compute_ppr_by_isi(groups(g), 'ppr_max');
        draw_ppr_series(ax3, S.isi, S, groups(g).color, groups(g).name, 's');
    end
    finish_ppr_axis(ax3, 'PPR (A2 / A1)');
    title(ax3, {'\itPer-ROI responder-gated \rm(diagnostic)', ...
        'P1 responders only \bullet selection-biased'}, 'Interpreter', 'tex');

    % Consistent y across the two PPR-ratio panels (1 & 3) with y-min pinned to 0.
    % Panel 2 is absolute charge on its own scale, so it is excluded.
    axs  = [ax1 ax3];
    ymax = max(arrayfun(@(a) a.YLim(2), axs));
    for a = axs, ylim(a, [0 ymax]); end

    sgtitle({'2AP Paired-Pulse Ratio (PPR = A2 / A1) vs ISI', ...
        'per-coverslip mean \pm SEM (dots = individual coverslips); dashed line = PPR 1 (no change)'}, ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'tex');
end

function draw_ppr_series(ax, isi, M, color, name, marker)
    % Faint per-coverslip points (jittered in log-x) + mean +/- SEM across
    % coverslips. M must have fields cs_vals (cell), mean, sem aligned to isi.
    if nargin < 6, marker = 'o'; end
    for k = 1:numel(isi)
        vals = M.cs_vals{k};
        if isempty(vals), continue; end
        jx = isi(k) .* (1 + (rand(1, numel(vals)) - 0.5) * 0.06);
        scatter(ax, jx, vals, 28, color, 'filled', 'MarkerFaceAlpha', 0.30, ...
            'MarkerEdgeColor', 'none', 'HandleVisibility', 'off');
    end
    valid = isfinite(M.mean);
    if any(valid)
        errorbar(ax, isi(valid), M.mean(valid), M.sem(valid), ['-' marker], 'Color', color, ...
            'MarkerFaceColor', color, 'LineWidth', 1.8, 'MarkerSize', 6, ...
            'CapSize', 4, 'DisplayName', name);
    end
end

function finish_ppr_axis(ax, ylab)
    set(ax, 'XScale', 'log');
    plot(ax, xlim(ax), [1 1], 'k--', 'HandleVisibility', 'off');   % PPR = 1 (no change)
    xlabel(ax, 'Inter-stimulus interval (ms)'); ylabel(ax, ylab);
    grid(ax, 'on'); box(ax, 'on'); legend(ax, 'Location', 'best');
end

function fig = plot_twoap_amplitudes(groups, config)
    % Pulse-1 vs pulse-2 amplitude at each ISI, from the coverslip grand-average
    % traces (ratio-of-means basis; unbiased, no responder gating). Grouped bars,
    % one panel per group, three bars per ISI:
    %   Pulse 1                 - A1, peak above the pre-stim1 baseline
    %   Pulse 2 (incremental)   - A2, peak above the LOCAL pre-stim2 residual
    %                             (the extra dF/F pulse 2 adds on pulse-1's decay)
    %   Pulse 2 (from baseline) - peak above the ORIGINAL pre-stim1 baseline
    % Bars = mean across coverslips; error bars = SEM. Incremental < from-baseline
    % at short ISI (pulse 2 rides high on pulse-1 decay) and they converge at long
    % ISI once the residual settles.
    rom = arrayfun(@(gd) compute_rom_by_isi(gd, config), groups, 'UniformOutput', false);
    nG  = numel(groups);
    fig = figure('Name', '2AP Pulse amplitudes vs ISI', ...
        'Position', [110 140 max(820, 620*nG) 500], 'Color', 'w');
    c1 = [0.35 0.35 0.35];   % pulse 1 = neutral reference

    % Shared y-limit across panels so groups are directly comparable
    ymax = 0;
    for g = 1:nG
        for f = {'a1', 'a2', 'a2_abs'}
            m = rom{g}.(f{1}).mean; s = rom{g}.(f{1}).sem;
            ymax = max([ymax, m(isfinite(m)) + s(isfinite(m))]);
        end
    end
    if ~(ymax > 0), ymax = 1; end

    for g = 1:nG
        gd = groups(g); ax = subplot(1, nG, g); hold on;
        isi = rom{g}.isi(:);
        x   = 1:numel(isi);
        Y = [rom{g}.a1.mean(:), rom{g}.a2.mean(:), rom{g}.a2_abs.mean(:)];
        E = [rom{g}.a1.sem(:),  rom{g}.a2.sem(:),  rom{g}.a2_abs.sem(:)];
        b = bar(ax, x, Y, 'grouped', 'EdgeColor', 'none');
        c2 = gd.color; c2l = c2 + (1 - c2) * 0.55;   % lighter shade of group color
        b(1).FaceColor = c1; b(2).FaceColor = c2; b(3).FaceColor = c2l;
        b(1).DisplayName = 'Pulse 1';
        b(2).DisplayName = 'Pulse 2 (incremental)';
        b(3).DisplayName = 'Pulse 2 (from baseline)';
        % SEM whiskers centred on each grouped bar
        for k = 1:numel(b)
            xk = b(k).XEndPoints; yk = Y(:, k)'; ek = E(:, k)';
            v  = isfinite(yk) & isfinite(ek) & (ek > 0);
            if any(v)
                errorbar(ax, xk(v), yk(v), ek(v), 'k', 'LineStyle', 'none', ...
                    'LineWidth', 0.8, 'CapSize', 3, 'HandleVisibility', 'off');
            end
        end
        set(ax, 'XTick', x, 'XTickLabel', compose('%g', isi));
        xlabel(ax, 'Inter-stimulus interval (ms)');
        ylabel(ax, 'Peak amplitude (\DeltaF/F)', 'Interpreter', 'tex');
        title(ax, gd.name, 'Interpreter', 'none');
        grid(ax, 'on'); box(ax, 'on'); legend(ax, 'Location', 'best');
        ylim(ax, [0 ymax * 1.05]);
    end
    sgtitle({'2AP Pulse-1 vs Pulse-2 amplitude vs ISI (all ROIs; per-coverslip mean \pm SEM)', ...
        'from coverslip avg trace. Pulse 2: incremental (above pre-stim2 residual) vs absolute (from pre-stim1 baseline)'}, ...
        'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'tex');
end

function G = compute_gating_bias(gd, config)
    % PPR vs ISI under responder gates, all on the SAME average-trace estimator as
    % the primary ROM figure. For each coverslip x ISI the group average TRACE is
    % restricted to an ROI subset, then peak2/peak1 of that trace is the PPR:
    %   no_gate     - all QC-passed ROIs  (== the primary ROM curve, "PPR vs ISI")
    %   gate_A1     - pulse-1 responders  (selecting large A1 inflates the mean
    %                 trace's pulse-1 peak but not pulse-2 -> biases PPR DOWN)
    %   gate_A2     - pulse-2 responders  (inflates pulse-2 -> biases PPR UP)
    %   gate_either / gate_both
    % Values are the per-coverslip mean across coverslips, so no_gate matches the
    % primary figure EXACTLY. gate_A1 < no_gate < gate_A2 is the regression-to-the-
    % mean signature that justifies the ungated ratio-of-means as the primary.
    conds   = {'no_gate', 'gate_A1', 'gate_A2', 'gate_either', 'gate_both'};
    srcName = {'all', 'gate_A1', 'gate_A2', 'gate_either', 'gate_both'};   % cs_avg_gate fields
    stim1   = config.twoAP.stim1_frame;
    G = struct('isi', [], 'labels', {conds});
    for c = 1:numel(conds), G.(conds{c}) = struct('ppr', [], 'n', []); end
    if ~isfield(gd, 'cs_avg_gate') || isempty(gd.cs_avg_gate), return; end

    ca   = gd.cs_avg_gate;
    isiE = [ca.isi];
    uISI = unique(isiE);
    G.isi = uISI;
    for c = 1:numel(conds)
        pv = nan(1, numel(uISI)); nv = zeros(1, numel(uISI));
        for k = 1:numel(uISI)
            idxs = find(isiE == uISI(k));
            vals = []; ntot = 0;
            for ii = idxs
                tr = ca(ii).traces.(srcName{c});
                nn = ca(ii).n.(srcName{c});
                if nn > 0 && any(isfinite(tr))
                    pk = trace_ppr(tr, stim1, ca(ii).stim2, config);   % peak2/peak1 of subset avg trace
                    if isfinite(pk), vals(end+1) = pk; ntot = ntot + nn; end %#ok<AGROW>
                end
            end
            if ~isempty(vals), pv(k) = mean(vals); nv(k) = ntot; end
        end
        G.(conds{c}).ppr = pv; G.(conds{c}).n = nv;
    end
end

function fig = plot_twoap_gating_bias(groups, config)
    % Selection-bias diagnostic: population PPR vs ISI under 5 responder gates
    % (see compute_gating_bias). One panel per group. The fanning of the gated
    % curves ABOVE and BELOW the unbiased no-gate line is the teaching point that
    % justifies using ratio-of-means as the primary estimator.
    nG = numel(groups);
    fig = figure('Name', '2AP Gating selection-bias diagnostic', ...
        'Position', [120 130 max(820, 600*nG) 500], 'Color', 'w');
    conds = {'no_gate', 'gate_A1', 'gate_A2', 'gate_either', 'gate_both'};
    labels = {'no gate (all ROIs = primary ROM)', 'gate on A1 (P1 responders)', ...
        'gate on A2 (P2 responders)', 'gate on either', 'gate on both'};
    cols = [0 0 0; 0.20 0.45 0.95; 0.90 0.30 0.20; 0.20 0.65 0.30; 0.60 0.30 0.70];
    lw   = [2.4 1.6 1.6 1.4 1.4];
    axs = gobjects(1, nG);
    for g = 1:nG
        gd = groups(g); G = compute_gating_bias(gd, config);
        ax = subplot(1, nG, g); hold on; axs(g) = ax;
        for c = 1:numel(conds)
            y = G.(conds{c}).ppr; v = isfinite(y);
            if any(v)
                plot(ax, G.isi(v), y(v), '-o', 'Color', cols(c, :), 'LineWidth', lw(c), ...
                    'MarkerFaceColor', cols(c, :), 'MarkerSize', 4, 'DisplayName', labels{c});
            end
        end
        set(ax, 'XScale', 'log');
        plot(ax, xlim(ax), [1 1], 'k--', 'HandleVisibility', 'off');
        xlabel(ax, 'Inter-stimulus interval (ms)'); ylabel(ax, 'PPR (A2 / A1 of avg trace)', 'Interpreter', 'tex');
        title(ax, gd.name, 'Interpreter', 'none'); grid(ax, 'on'); box(ax, 'on');
        legend(ax, 'Location', 'best', 'FontSize', 8);
    end
    % Consistent y across group panels, y-min pinned to 0 (matches the PPR figure)
    ymax = max(arrayfun(@(a) a.YLim(2), axs));
    for a = axs, ylim(a, [0 ymax]); end
    sgtitle({'2AP gating selection-bias diagnostic (per-coverslip avg traces; no-gate \equiv primary ROM)', ...
        'gate-on-A1 < no-gate < gate-on-A2 \Rightarrow regression to the mean, not biology'}, ...
        'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'tex');
end

function print_gating_bias(groups, config)
    fprintf('\n=== 2AP GATING SELECTION-BIAS DIAGNOSTIC (PPR = peak2/peak1 of subset avg trace) ===\n');
    fprintf('(Per-coverslip avg traces; no-gate == primary ROM. gate-A1 < no-gate < gate-A2 = bias.)\n');
    ksig = config.twoAP.upper_threshold_sigma;
    for g = 1:numel(groups)
        gd = groups(g); G = compute_gating_bias(gd, config);
        fprintf('%s (P2 responder = peak >= %g*sigma):\n', gd.name, ksig);
        fprintf('   ISI      no-gate      gate-A1      gate-A2   gate-either    gate-both\n');
        for k = 1:numel(G.isi)
            fprintf('  %5g   %5.2f(%4d)  %5.2f(%4d)  %5.2f(%4d)  %5.2f(%4d)  %5.2f(%4d)\n', ...
                G.isi(k), ...
                G.no_gate.ppr(k),     G.no_gate.n(k), ...
                G.gate_A1.ppr(k),     G.gate_A1.n(k), ...
                G.gate_A2.ppr(k),     G.gate_A2.n(k), ...
                G.gate_either.ppr(k), G.gate_either.n(k), ...
                G.gate_both.ppr(k),   G.gate_both.n(k));
        end
    end
    fprintf('=========================================================================\n');
end

function v = getfield_default(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function figs = plot_twoap_montage(groups, config)
    % Per group, example pulse-1 responders from the shortest-ISI reference file,
    % showing both stimuli and both methods' picks (o = peak-max, x = time-locked)
    % so the peak-selection methods can be compared by eye.
    figs = [];
    frame_dur = config.frame_dur_ms;
    stim1     = config.twoAP.stim1_frame;
    sync_f    = max(1, round(config.twoAP.sync_latency_ms / frame_dur));
    maxROIs   = 20;
    tax       = make_time_axis(config);
    cv        = tax.from_ms;   % ms position -> display units (for plotting)

    for g = 1:numel(groups)
        gd  = groups(g);
        ref = gd.montage_ref;
        if isempty(ref.dfof) || isempty(ref.twoap), continue; end
        t2    = ref.twoap;
        stim2 = t2.stim2_frame;
        nCol  = size(ref.dfof, 2);
        resp  = find(t2.is_responder1); resp = resp(resp <= nCol);
        if isempty(resp), continue; end
        order = resp(1:min(maxROIs, numel(resp)));

        fig = figure('Name', sprintf('2AP Montage %s (ISI %g ms)', gd.name, ref.isi), ...
            'Position', [50 50 1850 1000], 'Color', 'w');
        figs = [figs, fig]; %#ok<AGROW>
        nR = numel(order); nc = 5; nr = ceil(nR / nc);
        tl = tiledlayout(fig, nr, nc, 'TileSpacing', 'tight', 'Padding', 'tight');

        pre_f  = round(config.twoAP.montage_pre_ms  / frame_dur);
        post_f = round(config.twoAP.montage_post_ms / frame_dur);
        % Default window spans pre-stim1 .. post-stim2; montage_range_frames
        % (if set in config) overrides it with an absolute frame interval.
        nF_full = size(ref.dfof, 1);
        if ~isempty(tax.montage_range_frames)
            r  = tax.montage_range_frames;
            f0 = max(round(min(r)), 1);
            f1 = min(round(max(r)), nF_full);
        else
            f0 = max(stim1 - pre_f, 1);
            f1 = min(stim2 + post_f, nF_full);
        end
        fidx  = f0:f1;                                    % initial visible window
        tfull = cv(((1:nF_full) - stim1) * frame_dur);    % FULL-trace time axis
        xwin  = cv(([f0 f1]    - stim1) * frame_dur);     % initial x-limits (scrollable)
        stim2_ms = (stim2 - stim1) * frame_dur;         % stim2 position (ms, for arithmetic)
        % Peak-search windows (the range the peak-max method scans) for shading
        isi_frames = stim2 - stim1;
        search_f   = max(1, round(config.twoAP.peak_search_ms / frame_dur));
        s1_end_ms  = min(search_f, max(isi_frames - 1, 1)) * frame_dur;   % pulse-1 (capped at ISI)
        s2_end_ms  = config.twoAP.peak_search_ms;                         % pulse-2

        for k = 1:nR
            roi = order(k);
            nexttile(tl, k); hold on;
            % Shade the peak-search windows behind the trace (blue=pulse1, red=pulse2)
            xregion(0, cv(s1_end_ms), 'FaceColor', [0.2 0.4 1.0], 'FaceAlpha', 0.06);
            xregion(cv(stim2_ms), cv(stim2_ms + s2_end_ms), 'FaceColor', [1.0 0.3 0.3], 'FaceAlpha', 0.06);
            % Plot the WHOLE trace, then focus the view on the stim window; the
            % rest of the trace is a horizontal drag away (see scroll_time_axis).
            plot(tfull, ref.dfof(:, roi), 'k-', 'LineWidth', 1);
            % yline/xline span the full axis height/width, so thresholds and stim
            % markers stay put while the user pans/zooms in x.
            if isfinite(t2.upper_threshold(roi))
                yline(t2.upper_threshold(roi), 'g:', 'LineWidth', 1);
            end
            if isfield(t2, 'lower_threshold') && isfinite(t2.lower_threshold(roi))
                yline(t2.lower_threshold(roi), ':', 'Color', [0.5 0.8 0.5], 'LineWidth', 0.8);
            end
            % Focus the vertical scale on the stim window (not the whole trace).
            seg = ref.dfof(fidx, roi); seg = seg(isfinite(seg));
            if ~isempty(seg)
                lo = min(seg); hi = max(seg); pad = 0.1 * max(hi - lo, eps);
                ylim([lo - pad, hi + pad]);
            end
            xline(0, '--', 'Color', 'b', 'LineWidth', 1);            % stim1
            xline(cv(stim2_ms), '--', 'Color', 'r', 'LineWidth', 1); % stim2
            nF = nF_full;
            % Method 1 peak-max: filled circles at each pulse's peak
            t1ms = t2.peak1_time_ms(roi);
            if isfinite(t1ms)
                f = stim1 + round(t1ms / frame_dur);
                if f <= nF, plot(cv(t1ms), ref.dfof(f, roi), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5); end
            end
            if isfinite(t2.peak2_time_ms(roi))
                f = stim2 + round(t2.peak2_time_ms(roi) / frame_dur);
                if f <= nF, plot(cv(stim2_ms + t2.peak2_time_ms(roi)), ref.dfof(f, roi), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5); end
            end
            % Method 2 sync time-lock: crosses at the fixed sync latency
            if stim1 + sync_f <= nF
                plot(cv(config.twoAP.sync_latency_ms), ref.dfof(stim1 + sync_f, roi), 'bx', 'MarkerSize', 7, 'LineWidth', 1.2);
            end
            if stim2 + sync_f <= nF
                plot(cv(stim2_ms + config.twoAP.sync_latency_ms), ref.dfof(stim2 + sync_f, roi), 'rx', 'MarkerSize', 7, 'LineWidth', 1.2);
            end
            % Method 3 p1-latency lock: pulse-2 read at pulse-1's own peak latency
            if isfinite(t1ms)
                f = stim2 + round(t1ms / frame_dur);
                if f <= nF, plot(cv(stim2_ms + t1ms), ref.dfof(f, roi), 'rd', 'MarkerSize', 7, 'LineWidth', 1.2); end
            end
            ttl = sprintf('ROI%d  PPR=%.2f Q=%.1f', roi, t2.ppr_max(roi), t2.q_pair(roi));
            title(ttl, 'FontSize', 8);
            set(gca, 'FontSize', 7); grid on; box on;
            scroll_time_axis(gca, xwin);   % show stim window; drag to scroll x
        end
        title(tl, sprintf(['2AP %s, ISI=%g ms \\bullet example pulse-1 responders (one reference file); ' ...
            'tile title = PPR / total charge Q (dF/F*ms) \\bullet drag to scroll x, whole trace plotted\n' ...
            'blue=stim1, red=stim2; shaded=peak-search window; o=peak-max, x=sync %gms, d=p1-lock; ' ...
            'green=upper %g\\sigma / lower %g\\sigma; x in %s'], gd.name, ref.isi, ...
            config.twoAP.sync_latency_ms, config.twoAP.upper_threshold_sigma, config.twoAP.lower_threshold_sigma, tax.unit), ...
            'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'tex');
    end
end

function figs = plot_twoap_cs_average(groups, config)
    % Coverslip x ISI grand-average dF/F (QC). One tile per coverslip (rows) and
    % ISI (cols); the trace is the mean over ALL ROIs so peak2/peak1 of the trace
    % is the unbiased ratio-of-means PPR (annotated PPR / total charge Q in each
    % title). This is the honest, selection-free view of the paired-pulse response.
    %
    % Split into TWO figures by ISI so the short- and long-interval panels each
    % get readable column widths: short ISI (30-120 ms) and long ISI (250-1000
    % ms). ISIs outside both bands (if any) are dropped from these figures.

    % Row order = coverslips (grouped by group); col order = sorted ISIs
    rows = {}; rowGroup = []; allISI = [];
    for g = 1:numel(groups)
        if ~isfield(groups(g), 'cs_avg') || isempty(groups(g).cs_avg), continue; end
        csn = unique({groups(g).cs_avg.cs}, 'stable');
        for c = 1:numel(csn), rows{end+1} = csn{c}; rowGroup(end+1) = g; end %#ok<AGROW>
        allISI = [allISI, [groups(g).cs_avg.isi]]; %#ok<AGROW>
    end
    isis = unique(allISI);
    figs = gobjects(1, 0);
    if isempty(rows), return; end

    bands = {[30 120], 'short ISI (30-120 ms)'; [250 1000], 'long ISI (250-1000 ms)'};
    for b = 1:size(bands, 1)
        lo = bands{b,1}(1); hi = bands{b,1}(2);
        sub = isis(isis >= lo & isis <= hi);
        if isempty(sub), continue; end
        f = render_cs_average_fig(groups, config, rows, rowGroup, sub, bands{b,2});
        if ~isempty(f) && isgraphics(f), figs(end+1) = f; end %#ok<AGROW>
    end
end

function fig = render_cs_average_fig(groups, config, rows, rowGroup, isis, bandLabel)
    % Render ONE coverslip grand-average figure for the given ISI subset (cols).
    frame_dur = config.frame_dur_ms;
    stim1     = config.twoAP.stim1_frame;
    pre_f     = round(config.twoAP.montage_pre_ms  / frame_dur);
    post_f    = round(config.twoAP.montage_post_ms / frame_dur);

    nR = numel(rows); nC = numel(isis);
    fig = figure('Name', sprintf('2AP Coverslip grand-average dF/F - %s', bandLabel), ...
        'Position', [40 40 max(1600, 240*nC) max(500, 230*nR)], 'Color', 'w');
    tl = tiledlayout(fig, nR, nC, 'TileSpacing', 'tight', 'Padding', 'tight');

    rowAxes = gobjects(nR, nC);
    for r = 1:nR
        g = rowGroup(r); gd = groups(g); csname = rows{r};
        for cix = 1:nC
            isi = isis(cix);
            axh = nexttile(tl, (r-1)*nC + cix); hold on;
            idx = find(strcmp({gd.cs_avg.cs}, csname) & [gd.cs_avg.isi] == isi, 1);
            if isempty(idx), axis off; continue; end
            rowAxes(r, cix) = axh;
            e = gd.cs_avg(idx);
            trace = e.trace; stim2 = e.stim2; nF = numel(trace);
            f0 = max(stim1 - pre_f, 1); f1 = min(stim2 + post_f, nF); fidx = f0:f1;
            t  = (fidx - stim1) * frame_dur; stim2_ms = (stim2 - stim1) * frame_dur;
            % Shade the AUC total-charge window as ONE continuous band from stim1
            % to auc_post_ms after stim2 -- the encompassing window q_pair actually
            % integrates over. No gap between the pulses, so the whole paired
            % response the Q value sums is visible (this is NOT the peak-search
            % window, which is per-pulse and drives only the peak PPR).
            auc_end_ms = stim2_ms + config.twoAP.auc_post_ms;
            xregion(0, auc_end_ms, 'FaceColor', [0.45 0.45 0.50], 'FaceAlpha', 0.13);
            yy = trace(fidx); ss = e.sem(fidx);
            fill([t, fliplr(t)], [yy + ss, fliplr(yy - ss)], gd.color, ...
                'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            % Individual region (c#) traces behind the mean, to show the
            % within-coverslip variability the grand-average averages over.
            regs = e.regions; nrg = numel(regs);
            for rr = 1:nrg
                if numel(regs(rr).trace) < f1, continue; end
                plot(t, regs(rr).trace(fidx), '-', 'Color', [0.6 0.6 0.6], ...
                    'LineWidth', 0.5, 'HandleVisibility', 'off');
            end
            plot(t, yy, '-', 'Color', gd.color, 'LineWidth', 1.3, 'HandleVisibility', 'off');
            % xline spans the full axis height, so the stim markers stay full-height
            % after the per-row y-axes are linked below.
            xline(0, '--', 'Color', 'b', 'LineWidth', 0.8);
            xline(stim2_ms, '--', 'Color', 'r', 'LineWidth', 0.8);

            [ppr_pk, q_pair] = trace_ppr(trace, stim1, stim2, config);
            title(sprintf('%s  ISI %g  PPR=%.2f Q=%.1f (%d c)', csname, isi, ppr_pk, q_pair, nrg), 'FontSize', 8);
            set(gca, 'FontSize', 7); grid on; box on;
            axh.YAxis.Exponent = 0;   % plain decimal y ticks (no 10^n exponent)
            ytickformat(axh, '%.3g');
            if cix == 1, ylabel(sprintf('%s / %s', gd.name, csname), 'FontSize', 8, 'FontWeight', 'bold'); end
            if r == nR, xlabel('ms from stim1', 'FontSize', 7); end
        end
    end
    % Share a y-scale within each coverslip ROW (coverslips differ ~2x in
    % amplitude, so a per-row scale keeps pulse-2 vs pulse-1 heights comparable
    % across ISIs without the small coverslip being squashed by a global scale).
    for r = 1:nR
        va = isgraphics(rowAxes(r, :), 'axes');
        if nnz(va) >= 2, linkaxes(rowAxes(r, va), 'y'); end
    end
    title(tl, ['2AP Coverslip grand-average dF/F \bullet ' bandLabel ' (mean of ALL ROIs \pm SEM; ' ...
        'grey = individual region (c#) traces, showing within-coverslip spread; ' ...
        'blue=stim1, red=stim2, shaded=AUC charge window [stim1, stim2+' ...
        num2str(config.twoAP.auc_post_ms) 'ms]; y-scale shared per coverslip row). ' ...
        'PPR = peak2/peak1 of the avg trace; Q = total charge (dF/F*ms) over that window'], ...
        'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'tex');
end

function [ppr_pk, q_pair, A1, A2, A2_abs] = trace_ppr(trace, stim1, stim2, config)
    % Paired-pulse metrics of a single (already-averaged) trace, on the same
    % windows/baselines as the per-ROI analyzer (so these equal the ratio-of-means
    % values). Returns:
    %   ppr_pk  - peak PPR = A2 / A1
    %   q_pair  - TOTAL CHARGE (dF/F*ms) over ONE window encompassing both pulses,
    %             [stim1, stim2 + auc_post_ms], above the single pre-stim1 baseline.
    %             Not a ratio; grows with ISI by construction (window = ISI + post).
    %   A1, A2  - pulse-1/2 peak amplitude (dF/F); A2 above the local pre-stim2 residual
    %   A2_abs  - pulse-2 peak above the ORIGINAL pre-stim1 baseline. A2_abs - A2 =
    %             how far pulse-1's decay has lifted the baseline by the time pulse 2
    %             fires; they converge once the residual settles (long ISI).
    fd = config.frame_dur_ms; nF = numel(trace);
    search = max(1, round(config.twoAP.peak_search_ms / fd));
    isi_f  = stim2 - stim1; s1 = min(search, max(isi_f - 1, 1));

    pre  = trace(1:max(stim1 - 1, 1)); pre = pre(isfinite(pre));
    base = mean(pre);
    seg  = trace(max(stim2 - 2, stim1 + 1):(stim2 - 1)); seg = seg(isfinite(seg));
    resid = mean(seg); if ~isfinite(resid), resid = base; end

    w1 = stim1:min(stim1 + s1, nF); w2 = stim2:min(stim2 + search, nF);
    peak2 = max(trace(w2));
    A1 = max(trace(w1)) - base; A2 = peak2 - resid; A2_abs = peak2 - base;
    ppr_pk = A2 / A1;

    % Total charge over the encompassing window (above the single pre-stim1 baseline)
    post = max(1, round(config.twoAP.auc_post_ms / fd));
    wq   = stim1:min(stim2 + post, nF);
    sq   = trace(wq) - base; sq = sq(isfinite(sq));
    if numel(sq) >= 2, q_pair = trapz(sq) * fd; else, q_pair = NaN; end
end

function S = compute_rom_by_isi(gd, config)
    % Ratio-of-means (unbiased) PPR + absolute pulse amplitudes at each ISI,
    % built from the coverslip grand-average traces (gd.cs_avg). Each cs_avg
    % entry is one coverslip at one ISI, so per ISI the values across entries are
    % the per-coverslip values. Returns S.isi and, for each metric, a sub-struct
    % with per-coverslip vectors (cs_vals) and mean/SEM/N across coverslips:
    %   ppr_pk           - peak-based ratio-of-means PPR (A2/A1 of the avg trace)
    %   q_pair           - total charge (dF/F*ms) over the encompassing window
    %   a1, a2           - pulse-1/2 peak amplitude (dF/F); a2 above pre-stim2 residual
    %   a2_abs           - pulse-2 peak above the ORIGINAL pre-stim1 baseline
    names = {'ppr_pk', 'q_pair', 'a1', 'a2', 'a2_abs'};
    empty = struct('cs_vals', {{}}, 'mean', [], 'sem', [], 'n_cs', []);
    S = struct('isi', []);
    for n = 1:numel(names), S.(names{n}) = empty; end
    if ~isfield(gd, 'cs_avg') || isempty(gd.cs_avg), return; end

    stim1 = config.twoAP.stim1_frame;
    ca    = gd.cs_avg;
    nE    = numel(ca);
    isiE  = [ca.isi];
    vals  = struct();
    for n = 1:numel(names), vals.(names{n}) = nan(1, nE); end
    for k = 1:nE
        [pk, q, A1, A2, A2abs] = trace_ppr(ca(k).trace, stim1, ca(k).stim2, config);
        vals.ppr_pk(k) = pk; vals.q_pair(k) = q;
        vals.a1(k) = A1; vals.a2(k) = A2; vals.a2_abs(k) = A2abs;
    end

    uISI  = unique(isiE);
    S.isi = uISI;
    for n = 1:numel(names)
        va = vals.(names{n});
        cv = cell(1, numel(uISI)); mn = nan(1, numel(uISI));
        se = nan(1, numel(uISI)); nc = zeros(1, numel(uISI));
        for k = 1:numel(uISI)
            v = va(isiE == uISI(k)); v = v(isfinite(v));
            cv{k} = v; nc(k) = numel(v);
            if ~isempty(v), mn(k) = mean(v); se(k) = std(v) / sqrt(numel(v)); end
        end
        S.(names{n}) = struct('cs_vals', {cv}, 'mean', mn, 'sem', se, 'n_cs', nc);
    end
end

%% ======================================================================
%% Output / utility
%% ======================================================================
function T = build_twoap_table(groups)
    grp = {}; drug = {}; fname = {}; cs = {}; roi = []; isi = []; resp1 = [];
    pmode = {}; ichk = {};
    a1 = []; a2 = []; pmax = []; a1s = []; a2s = []; psync = []; a2p = []; pp1 = []; sd = [];
    q = [];
    for g = 1:numel(groups)
        gd = groups(g);
        n = gd.total_rois;
        for j = 1:n
            grp{end+1,1} = gd.genotype; %#ok<AGROW>
            drug{end+1,1} = gd.drug; %#ok<AGROW>
            fidx = gd.file_of_roi(j);
            if fidx >= 1 && fidx <= numel(gd.filenames)
                fname{end+1,1} = gd.filenames{fidx}; %#ok<AGROW>
                pmode{end+1,1} = gd.file_ppf_mode{fidx}; %#ok<AGROW>
                ichk{end+1,1}  = gd.file_isi_check{fidx}; %#ok<AGROW>
            else
                fname{end+1,1} = ''; pmode{end+1,1} = ''; ichk{end+1,1} = ''; %#ok<AGROW>
            end
            cs{end+1,1} = gd.coverslip_of_roi{j}; %#ok<AGROW>
        end
        roi   = [roi;   gd.roi_index(:)];        %#ok<AGROW>
        isi   = [isi;   gd.isi_of_roi(:)];        %#ok<AGROW>
        resp1 = [resp1; gd.is_responder1(:)];     %#ok<AGROW>
        a1    = [a1;    gd.peak1_amp(:)];          %#ok<AGROW>
        a2    = [a2;    gd.peak2_amp(:)];          %#ok<AGROW>
        pmax  = [pmax;  gd.ppr_max(:)];            %#ok<AGROW>
        a1s   = [a1s;   gd.a1_sync(:)];            %#ok<AGROW>
        a2s   = [a2s;   gd.a2_sync(:)];            %#ok<AGROW>
        psync = [psync; gd.ppr_sync(:)];           %#ok<AGROW>
        a2p   = [a2p;   gd.a2_p1lock(:)];          %#ok<AGROW>
        pp1   = [pp1;   gd.ppr_p1lock(:)];         %#ok<AGROW>
        q     = [q;     gd.q_pair(:)];             %#ok<AGROW>
        sd    = [sd;    gd.prestim_sd(:)];         %#ok<AGROW>
    end
    T = table(grp, drug, fname, cs, roi, isi, pmode, logical(resp1), a1, a2, pmax, a1s, a2s, psync, a2p, pp1, q, sd, ichk, ...
        'VariableNames', {'genotype', 'drug', 'filename', 'coverslip', 'roi', 'isi_ms', 'ppf_mode', 'is_responder1', ...
        'peak1_amp', 'peak2_amp', 'ppr_max', 'a1_sync', 'a2_sync', 'ppr_sync', ...
        'a2_p1lock', 'ppr_p1lock', 'q_pair', 'prestim_sd', 'isi_check'});
end

function print_summary_2ap(groups, config)
    fprintf('\n=== 2AP SUMMARY (PPR = A2/A1; per-coverslip mean at each ISI) ===\n');
    fprintf('(Statistical unit = coverslip. Regions/files within a coverslip are\n');
    fprintf(' pseudoreplicates. With 1 animal/genotype these are technical replicates.)\n');
    fprintf('PRIMARY = ratio-of-means (peak2/peak1 of the coverslip average trace,\n');
    fprintf(' unbiased). gated = per-ROI responder-gated mean (selection-biased low).\n');
    for g = 1:numel(groups)
        gd   = groups(g);
        rom  = compute_rom_by_isi(gd, config);
        Smax = compute_ppr_by_isi(gd, 'ppr_max');
        fprintf('%s: %d files, %d coverslips, %d pulse-1 responders\n', ...
            gd.name, gd.num_files, gd.num_coverslips, gd.n_responders1);
        for k = 1:numel(gd.isi_list)
            isi = gd.isi_list(k);
            [mPk, nPk] = mean_at_rom(rom.ppr_pk, rom.isi, isi);
            mQ         = mean_at_rom(rom.q_pair, rom.isi, isi);
            [mA1, ~]   = mean_at_rom(rom.a1, rom.isi, isi);
            mA2        = mean_at_rom(rom.a2, rom.isi, isi);
            mGate      = mean_at_isi(Smax, isi);
            fprintf(['   ISI %5g ms: ROM peak %.2f (N=%d cs) | Q %.1f | ' ...
                'A1=%.4f A2=%.4f | gated %.2f\n'], ...
                isi, mPk, nPk, mQ, mA1, mA2, mGate);
        end
    end
    fprintf('================================================================\n');
end

function [m, n] = mean_at_isi(S, isi)
    % Look up a method's per-coverslip mean (and N coverslips) at a given ISI.
    m = NaN; n = 0;
    k = find(S.isi == isi, 1);
    if ~isempty(k), m = S.mean(k); n = S.n_cs(k); end
end

function [m, n] = mean_at_rom(M, isiList, isi)
    % Look up a ratio-of-means metric sub-struct's mean (and N coverslips) at ISI.
    m = NaN; n = 0;
    k = find(isiList == isi, 1);
    if ~isempty(k), m = M.mean(k); n = M.n_cs(k); end
end

function config = local_config()
    evalc('config = tracenorm_config();');
    config.dataType = '2AP';
end
