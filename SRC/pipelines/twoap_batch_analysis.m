function results = twoap_batch_analysis(dataset, options)
    % TWOAP_BATCH_ANALYSIS - Batch analysis of 2AP (paired-pulse / PPF) data
    %
    % Counterpart to oneap_batch_analysis for paired-pulse recordings. Takes a
    % dataset from load_dataset_folder: each file is one region at one ISI, with
    % stim1, the ISI and the frame interval from Recordings.csv. All files must
    % share one frame interval (run_pipeline splits mixed folders first). Files
    % are grouped by genotype x drug, processed through twoap_pipeline, and the
    % paired-pulse ratio (PPR = A2/A1) is aggregated PER COVERSLIP at each ISI.
    % PRIMARY PPR = pooled SumA2 / SumA1 over the gate-either ROIs of each coverslip
    % x ISI (pulse-1 responder OR incremental A2 >= upper_sigma*sd), with A2
    % decay-subtracted (pulse 1's tail, projected from the group's longest-ISI
    % average trace, removed under pulse 2). The ungated average-trace ratio is
    % kept as a QC reference. Produces the PPR-vs-ISI figure, estimator / gating
    % diagnostics, decay-subtraction examples, a per-ROI metrics table and a
    % summary. Files whose stims don't match the reported ISI are listed (isi_check).
    %
    % Inputs:
    %   dataset - struct from load_dataset_folder (one frame interval)
    %   options - struct: verbose, createPlots, createMontage, saveFigures,
    %             diagnosticPlots (default false)
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
    % diagnosticPlots: estimator / gating diagnostics used to choose the primary
    % PPR (gating-bias, estimator comparison, per-ROI distributions + their
    % printed tables). Off by default; the primary PPR is computed either way.
    if ~isfield(options, 'diagnosticPlots'), options.diagnosticPlots = false; end

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
            gd.ppr_diag = compute_ppr_diagnostic(gd, fileResults(mask), config);
            gd.sync_async = compute_sync_async(gd, config);
            if gd.sync_async.enabled
                gd.sync_async.examples = pick_sync_async_examples(fileResults(mask), config);
            end
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
            'primary',        groups(g).ppr_diag.decay_sub.pooled, ...    % PRIMARY: gate-either SumA2/SumA1, decay-sub
            'ref_no_gate',    groups(g).ppr_diag.decay_sub.trace_all, ... % QC: ungated avg trace, decay-sub
            'ratio_of_means', compute_rom_by_isi(groups(g), config), ...  % ungated avg trace, peak-max
            'peak_max',   compute_ppr_by_isi(groups(g), 'ppr_max'), ...   % per-ROI diagnostic
            'p1lock',     compute_ppr_by_isi(groups(g), 'ppr_p1lock'), ...
            'sync',       compute_ppr_by_isi(groups(g), 'ppr_sync'), ...
            'sync_async', groups(g).sync_async);                         % event binning (fast frames only)
    end

    %% Step 5: plots
    plot_handles = struct();
    if options.createPlots
        fprintf('\nCreating 2AP figures...\n');
        plot_handles.ppr_vs_isi   = plot_twoap_ppr(groups, config);
        plot_handles.amplitudes   = plot_twoap_amplitudes(groups, config);
        plot_handles.amplitude_hist = plot_twoap_amplitude_hist(groups, config);
        if options.diagnosticPlots
            plot_handles.gating_bias  = plot_twoap_gating_bias(groups, config);
            plot_handles.ppr_diag     = plot_twoap_ppr_diagnostic(groups, config);
            plot_handles.ppr_dist     = plot_twoap_ppr_distributions(groups, config);
        end
        plot_handles.decay_ex     = plot_twoap_decay_examples(groups, config);
        plot_handles.cs_average   = plot_twoap_cs_average(groups, config);
        plot_handles.sync_async   = plot_twoap_sync_async(groups, config);   % [] at slow frames
        plot_handles.sync_async_ex = plot_twoap_sync_async_examples(groups, config);
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
        n_saved = save_open_figures(output_dir, config, figs);
        fprintf('  Saved %d figures to:\n    %s\n', n_saved, output_dir);
    end

    %% Step 6: per-ROI metrics table
    twoap_metrics_table   = build_twoap_table(groups);
    twoap_coverslip_table = build_twoap_cs_table(groups);
    twoap_sync_async_table = build_twoap_sync_async_table(groups);

    %% Assemble results
    results = struct();
    results.dataType            = '2AP';
    results.groups              = groups;
    results.ppr_by_isi          = ppr_by_isi;
    results.plot_handles        = plot_handles;
    results.twoap_metrics_table = twoap_metrics_table;
    results.twoap_coverslip_table = twoap_coverslip_table;
    results.twoap_sync_async_table = twoap_sync_async_table;   % empty at slow frames
    results.isi_checks          = isi_checks;
    results.output_dir          = output_dir;   % '' unless options.saveFigures
    results.processing_info     = struct('folder_path', folder_path, ...
        'metadata_file', dataset.metadata_file, ...
        'process_time', process_time, 'timestamp', datetime('now'), ...
        'frame_interval_ms', config.frame_dur_ms, 'epoch_frames', config.expected_frames, ...
        'epoch_stim1_frame', config.twoAP.stim1_frame);

    print_summary_2ap(groups, config);
    print_sync_async(groups, config);
    if options.diagnosticPlots
        print_gating_bias(groups, config);
        print_ppr_diagnostic(groups, config);
    end
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
    is_responder1 = []; is_responder2 = []; prestim_sd = []; qc_roi = false(1, 0);
    % Sync / async events: per-ROI counts, per-ROI bin lengths (ms) and one flat
    % event list (ROI index into this group's ROI order).
    ev = struct('enabled', false, 'n_sync1', [], 'n_async1', [], 'n_sync2', [], ...
        'n_async2', [], 'n_spont', [], 'async1_ms', [], 'async2_ms', [], 'spont_ms', [], ...
        'list', struct('roi', [], 'pulse', [], 'lat_ms', [], 'onset_ms', [], 'prom', [], ...
                       'is_sync', false(1, 0)));
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
        qc_roi        = [qc_roi,        isfinite(t.upper_threshold)]; %#ok<AGROW>
        ev = append_events(ev, t, numel(isi_of_roi));
        isi_of_roi    = [isi_of_roi,    repmat(fr.isi_ms, 1, nR)]; %#ok<AGROW>
        coverslip_of_roi = [coverslip_of_roi, repmat({parse_coverslip_id(fr.filename)}, 1, nR)]; %#ok<AGROW>
        file_of_roi   = [file_of_roi,   repmat(i, 1, nR)]; %#ok<AGROW>
        roi_index     = [roi_index,     1:nR];             %#ok<AGROW>

        % Accumulate coverslip x ISI average traces, nan-aware, over QC-passed ROIs
        % (analyzable + not noise-excluded; NOT filtered on response, so no selection
        % bias). Also keep the same sums restricted to responder subsets, so the
        % gating-bias diagnostic and the PPR figure's QC reference use the SAME
        % average-trace estimator: 'all' -> no-gate (panel 3); gate_A1 -> P1 responders;
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
        cond.excluded    = qc & ~(r1raw | r2raw);        % sub-threshold: fail both gates
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
                'stim2', t.stim2_frame, 'pre', t.prestim_idx, 'sq', sq, 's', sS, 'c', cS, ...
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
    gd.qc               = qc_roi;                 % passed noise QC (analyzable, not noise-excluded)
    gd.events           = ev;
    gd.isi_of_roi       = isi_of_roi;
    gd.coverslip_of_roi = coverslip_of_roi;
    gd.file_of_roi      = file_of_roi;
    gd.filenames        = filenames;
    gd.file_isi_check   = file_isi_check;
    gd.roi_index        = roi_index;              % ROI number within its file
    gd.file_ppf_mode    = file_ppf_mode;
    gd.montage_ref      = montage_ref;

    % Finalize coverslip x ISI average traces. cs_avg = the QC-passed grand average
    % (mean; SEM across ROIs kept but not plotted), used by the charge panel, the
    % decay template and the coverslip-average figure. cs_avg_gate holds
    % the same average restricted to each responder subset, for the gating diagnostic.
    cs_avg      = struct('cs', {}, 'isi', {}, 'trace', {}, 'sem', {}, 'n', {}, 'stim2', {}, 'pre', {}, 'regions', {});
    cs_avg_gate = struct('cs', {}, 'isi', {}, 'stim2', {}, 'pre', {}, 'traces', {}, 'n', {});
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
            'sem', sem, 'n', round(median(ca)), 'stim2', e.stim2, 'pre', e.pre, 'regions', regions);
        cnames = fieldnames(e.s);
        traces = struct(); ncnt = struct();
        for q = 1:numel(cnames)
            cc = e.c.(cnames{q});
            traces.(cnames{q}) = (e.s.(cnames{q}) ./ max(cc, 1))';   % 1 x nF subset mean
            ncnt.(cnames{q})   = round(median(cc));
        end
        cs_avg_gate(end+1) = struct('cs', e.cs, 'isi', e.isi, 'stim2', e.stim2, 'pre', e.pre, ...  %#ok<AGROW>
            'traces', traces, 'n', ncnt);
    end
    gd.cs_avg           = cs_avg;
    gd.cs_avg_gate      = cs_avg_gate;

    gd.total_rois       = numel(ppr_max);
    gd.n_responders1    = sum(gd.is_responder1);
    gd.num_coverslips   = numel(unique(coverslip_of_roi));
    gd.isi_list         = unique(isi_of_roi);
end

function ev = append_events(ev, t, offset)
    % Add one file's sync / async event counts and event list to the group's
    % (offset = number of ROIs already in the group).
    nR = t.numROIs;
    ev.enabled  = ev.enabled || t.events_enabled;
    ev.n_sync1  = [ev.n_sync1,  t.n_sync1];
    ev.n_async1 = [ev.n_async1, t.n_async1];
    ev.n_sync2  = [ev.n_sync2,  t.n_sync2];
    ev.n_async2 = [ev.n_async2, t.n_async2];
    ev.n_spont  = [ev.n_spont,  t.n_spont];
    ev.async1_ms = [ev.async1_ms, repmat(t.events_ms.async1, 1, nR)];
    ev.async2_ms = [ev.async2_ms, repmat(t.events_ms.async2, 1, nR)];
    ev.spont_ms  = [ev.spont_ms,  repmat(t.events_ms.spont,  1, nR)];
    L = t.events;
    ev.list.roi     = [ev.list.roi,     L.roi + offset];
    ev.list.pulse   = [ev.list.pulse,   L.pulse];
    ev.list.lat_ms  = [ev.list.lat_ms,  L.lat_ms];
    ev.list.onset_ms = [ev.list.onset_ms, L.onset_ms];
    ev.list.prom    = [ev.list.prom,    L.prom];
    ev.list.is_sync = [ev.list.is_sync, L.is_sync];
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
    % PPR vs ISI. PRIMARY = pooled SumA2 / SumA1 over each coverslip's gate-either
    % ROIs (pulse-1 responder OR incremental A2 >= upper_sigma*sd), A2
    % decay-subtracted (gd.ppr_diag.decay_sub.pooled). Gating on either pulse is
    % close to symmetric, and dividing the sums once avoids the divide-by-small-A1
    % blow-up of a mean of per-ROI ratios. Panel 3 is the QC reference: the
    % ungated average-trace ratio on the same (decay-subtracted) amplitudes; a
    % large gap between panels 1 and 3 means the gate is biasing the PPR.
    fig = figure('Name', '2AP PPR vs ISI', 'Position', [90 120 1600 520], 'Color', 'w');

    % Ratio-of-means per group (one value per coverslip per ISI)
    rom = arrayfun(@(gd) compute_rom_by_isi(gd, config), groups, 'UniformOutput', false);

    % --- Panel 1: pooled gate-either, decay-subtracted (PRIMARY) ---
    ax1 = subplot(1, 3, 1); hold on;
    for g = 1:numel(groups)
        S = groups(g).ppr_diag.decay_sub.pooled;
        draw_ppr_series(ax1, S.isi, S, groups(g).color, groups(g).name, 'o');
    end
    % Sub-threshold ROIs (fail both 5-sigma gates): a separate population, PPR from
    % their average trace (single ROIs are too noisy for per-ROI peaks).
    for g = 1:numel(groups)
        S  = groups(g).ppr_diag.sub_threshold;
        cl = groups(g).color + (1 - groups(g).color) * 0.45;
        draw_ppr_series(ax1, S.isi, S, cl, [groups(g).name ' (sub-threshold)'], 'd', '--');
    end
    finish_ppr_axis(ax1, 'PPR (A2 / A1)');
    title(ax1, {'\bfGate-either \Sigma A2/\Sigma A1 \rm(primary)', ...
        'dashed: sub-threshold ROIs (avg trace)'}, 'Interpreter', 'tex');

    % --- Panel 2: total charge over the encompassing window (absolute, PRIMARY) ---
    % Charge in dF/F*ms, its own y-scale (NOT a ratio, so no PPR=1 line and not
    % shared with the ratio panels). The window lengthens with ISI, so baseline
    % drift accumulates into Q at long ISI (raw Q-vs-ISI is drift-confounded).
    ax2 = subplot(1, 3, 2); hold on;
    for g = 1:numel(groups)
        draw_ppr_series(ax2, rom{g}.isi, rom{g}.q_pair, groups(g).color, groups(g).name, 'o');
    end
    set(ax2, 'XScale', 'log');
    xlabel(ax2, 'Inter-stimulus interval (ms)');
    ylabel(ax2, 'Total charge Q (dF/F\cdotms)');
    grid(ax2, 'on'); box(ax2, 'on'); legend(ax2, 'Location', 'best');
    ylim(ax2, [0 max(ax2.YLim(2), eps)]);
    title(ax2, {'\bfTotal charge \rm(primary)', ...
        sprintf('\\intdF/F over [stim1, stim2+%gms] of avg trace \\bullet all ROIs', ...
        config.twoAP.auc_post_ms)}, 'Interpreter', 'tex');

    % --- Panel 3: ungated average-trace ratio, decay-subtracted (QC reference) ---
    ax3 = subplot(1, 3, 3); hold on;
    for g = 1:numel(groups)
        S = groups(g).ppr_diag.decay_sub.trace_all;
        draw_ppr_series(ax3, S.isi, S, groups(g).color, groups(g).name, 's');
    end
    finish_ppr_axis(ax3, 'PPR (A2 / A1 of avg trace)');
    title(ax3, {'\itNo gate, avg trace \rm(QC reference)', ...
        'all QC-passed ROIs \bullet A2 decay-subtracted'}, 'Interpreter', 'tex');

    % Consistent y across the two PPR-ratio panels (1 & 3) with y-min pinned to 0.
    % Panel 2 is absolute charge on its own scale, so it is excluded.
    axs  = [ax1 ax3];
    ymax = max(arrayfun(@(a) a.YLim(2), axs));
    for a = axs, ylim(a, [0 ymax]); end

    sgtitle({'2AP Paired-Pulse Ratio (PPR = A2 / A1) vs ISI', ...
        'per-coverslip mean \pm SEM (dots = individual coverslips); dashed line = PPR 1 (no change)'}, ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'tex');
end

function draw_ppr_series(ax, isi, M, color, name, marker, lineStyle)
    % Faint per-coverslip points (jittered in log-x) + mean +/- SEM across
    % coverslips. M must have fields cs_vals (cell), mean, sem aligned to isi.
    if nargin < 6, marker = 'o'; end
    if nargin < 7, lineStyle = '-'; end
    for k = 1:numel(isi)
        vals = M.cs_vals{k};
        if isempty(vals), continue; end
        jx = isi(k) .* (1 + (rand(1, numel(vals)) - 0.5) * 0.06);
        scatter(ax, jx, vals, 28, color, 'filled', 'MarkerFaceAlpha', 0.30, ...
            'MarkerEdgeColor', 'none', 'HandleVisibility', 'off');
    end
    valid = isfinite(M.mean);
    if any(valid)
        errorbar(ax, isi(valid), M.mean(valid), M.sem(valid), [lineStyle marker], 'Color', color, ...
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
    % Pulse-1 vs pulse-2 amplitude at each ISI on the SAME basis as the primary
    % PPR: per coverslip, the mean over the gate-either ROIs of A1 (pulse-1 peak
    % above the pre-stim1 baseline) and of the decay-subtracted A2 (pulse-2 peak
    % after pulse 1's projected tail is removed). Their ratio is the primary PPR.
    % A third bar shows the ABSOLUTE pulse-2 peak (above the pre-stim1 baseline, no
    % subtraction): absolute - decay-subtracted = how much of pulse 2's height is
    % pulse 1 still decaying underneath it (large at short ISI, ~0 at long ISI).
    % Grouped bars, one panel per group; bar = mean across coverslips, whisker = SEM.
    nG  = numel(groups);
    fig = figure('Name', '2AP Pulse amplitudes vs ISI', ...
        'Position', [110 140 max(820, 620*nG) 500], 'Color', 'w');
    c1 = [0.35 0.35 0.35];   % pulse 1 = neutral reference

    % Shared y-limit across panels so groups are directly comparable
    ymax = 0;
    for g = 1:nG
        D = groups(g).ppr_diag;
        for M = [D.amp_a1, D.amp_a2, D.amp_a2_abs]
            v = isfinite(M.mean);
            ymax = max([ymax, M.mean(v) + max(M.sem(v), 0)]);
        end
    end
    if ~(ymax > 0), ymax = 1; end

    for g = 1:nG
        gd = groups(g); D = gd.ppr_diag; ax = subplot(1, nG, g); hold on;
        isi = D.amp_a1.isi(:);
        x   = 1:numel(isi);
        Y = [D.amp_a1.mean(:), D.amp_a2.mean(:), D.amp_a2_abs.mean(:)];
        E = [D.amp_a1.sem(:),  D.amp_a2.sem(:),  D.amp_a2_abs.sem(:)];
        b = bar(ax, x, Y, 'grouped', 'EdgeColor', 'none');
        b(1).FaceColor = c1; b(2).FaceColor = gd.color;
        b(3).FaceColor = gd.color + (1 - gd.color) * 0.55;   % lighter shade
        b(1).DisplayName = 'Pulse 1 (A1)';
        b(2).DisplayName = 'Pulse 2 (A2, decay-subtracted)';
        b(3).DisplayName = 'Pulse 2 (absolute, from pre-stim1 baseline)';
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
        ylabel(ax, 'Mean peak amplitude per ROI (\DeltaF/F)', 'Interpreter', 'tex');
        title(ax, sprintf('%s (decay template %g ms ISI)', gd.name, D.template_isi), 'Interpreter', 'none');
        grid(ax, 'on'); box(ax, 'on'); legend(ax, 'Location', 'best');
        ylim(ax, [0 ymax * 1.05]);
    end
    sgtitle({sprintf(['2AP Pulse-1 vs Pulse-2 amplitude vs ISI (gate-either ROIs: P1 responder OR ' ...
        'incremental A2 \\geq %g\\sigma)'], config.twoAP.upper_threshold_sigma), ...
        ['per coverslip: mean over ROIs (decay-subtracted A2 / A1 = primary PPR); absolute = raw pulse-2 ' ...
        'peak height; bars = mean across coverslips \pm SEM']}, ...
        'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'tex');
end

function figs = plot_twoap_amplitude_hist(groups, config)
    % Per-ROI amplitude distributions of the gate-either ROIs (the PRIMARY PPR
    % population), one figure per pulse: pulse 1 (A1), then pulse 2 (decay-
    % subtracted A2). Each tile is one ISI with the groups side by side, so the
    % comparison is between genotypes at the same pulse, not pulse 1 vs pulse 2
    % within a genotype. ROIs pooled across coverslips: distribution shape only.
    % Bins are shared across both figures and all ISIs so every tile compares directly.
    figs = gobjects(1, 0);
    has = arrayfun(@(gd) isfield(gd, 'ppr_diag') && ~isempty(gd.ppr_diag) && ...
        ~isempty(gd.ppr_diag.isi), groups);
    groups = groups(has);
    if isempty(groups), return; end
    allv = [];
    for g = 1:numel(groups)
        r = groups(g).ppr_diag.roi;
        allv = [allv, r.a1.decay_sub, r.a2.decay_sub]; %#ok<AGROW>
    end
    allv = allv(isfinite(allv));
    if isempty(allv), return; end
    lim   = prctile(allv, [0.5 99.5]);               % keep a few outliers from setting the range
    edges = linspace(min(0, lim(1)), lim(2), 31);
    isiC   = arrayfun(@(gd) gd.ppr_diag.isi(:)', groups, 'UniformOutput', false);
    isiAll = unique([isiC{:}]);
    nI = numel(isiAll);
    nc = min(nI, 4); nr = ceil(nI / nc);
    colors = vertcat(groups.color);
    pulses = {'a1', 'Pulse 1 (A1)'; 'a2', 'Pulse 2 (A2, decay-subtracted)'};
    for p = 1:size(pulses, 1)
        fig = figure('Name', sprintf('2AP amplitude distributions - %s', pulses{p, 2}), ...
            'Position', [90 90 360*nc + 100, 300*nr + 120], 'Color', 'w');
        figs(end+1) = fig; %#ok<AGROW>
        tl = tiledlayout(fig, nr, nc, 'TileSpacing', 'compact', 'Padding', 'compact');
        for k = 1:nI
            ax = nexttile(tl); hold(ax, 'on');
            vals = cell(1, numel(groups)); names = cell(1, numel(groups));
            for g = 1:numel(groups)
                r = groups(g).ppr_diag.roi;
                v = r.(pulses{p, 1}).decay_sub(r.isi == isiAll(k));
                vals{g}  = v;
                names{g} = sprintf('%s (n=%d)', groups(g).name, sum(isfinite(v)));
            end
            grouped_hist(ax, vals, edges, colors, names);
            xline(ax, 0, 'k:', 'HandleVisibility', 'off');
            xlim(ax, edges([1 end])); grid(ax, 'on'); box(ax, 'on');
            title(ax, sprintf('ISI %g ms', isiAll(k)));
            legend(ax, 'Location', 'northeast', 'Interpreter', 'none');
            xlabel(ax, 'Peak amplitude per ROI (\DeltaF/F)');
            ylabel(ax, 'Relative fraction of ROIs');
        end
        title(tl, {sprintf('2AP per-ROI amplitude distributions \\bullet %s', pulses{p, 2}), ...
            sprintf(['gate-either ROIs (P1 responder OR incremental A2 \\geq %g\\sigma), pooled ' ...
            'across coverslips \\bullet x clipped to the 0.5-99.5th percentile'], ...
            config.twoAP.upper_threshold_sigma)}, ...
            'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'tex');
    end
end

function G = compute_gating_bias(gd, config)
    % PPR vs ISI under responder gates, on the average-trace estimator with the
    % SAME decay-subtracted A2 as the PPR figure. For each coverslip x ISI the
    % group average TRACE is restricted to an ROI subset, then A2/A1 of that trace:
    %   no_gate     - all QC-passed ROIs  (== PPR figure panel 3, QC reference)
    %   gate_A1     - pulse-1 responders  (selecting large A1 inflates the mean
    %                 trace's pulse-1 peak but not pulse-2 -> biases PPR DOWN)
    %   gate_A2     - pulse-2 responders  (inflates pulse-2 -> biases PPR UP)
    %   gate_either / gate_both
    % Values are the per-coverslip mean across coverslips, so no_gate matches
    % panel 3 EXACTLY. gate_A1 < no_gate < gate_A2 is the regression-to-the-mean
    % signature. Peak-max A2 is used only when the group has no decay template.
    conds   = {'no_gate', 'gate_A1', 'gate_A2', 'gate_either', 'gate_both'};
    srcName = {'all', 'gate_A1', 'gate_A2', 'gate_either', 'gate_both'};   % cs_avg_gate fields
    stim1   = config.twoAP.stim1_frame;
    T       = [];
    if isfield(gd, 'ppr_diag') && isfield(gd.ppr_diag, 'template'), T = gd.ppr_diag.template; end
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
                    pk = trace_ppr_primary(tr, stim1, ca(ii).stim2, config, ca(ii).pre, T);
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
    % shows how much a single-pulse responder gate biases the PPR.
    nG = numel(groups);
    fig = figure('Name', '2AP Gating selection-bias diagnostic', ...
        'Position', [120 130 max(820, 600*nG) 500], 'Color', 'w');
    conds = {'no_gate', 'gate_A1', 'gate_A2', 'gate_either', 'gate_both'};
    labels = {'no gate (all ROIs = PPR fig panel 3)', 'gate on A1 (P1 responders)', ...
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
    sgtitle({['2AP gating selection-bias diagnostic (per-coverslip avg traces, A2 decay-subtracted; ' ...
        'no-gate \equiv PPR fig panel 3)'], ...
        'gate-on-A1 < no-gate < gate-on-A2 \Rightarrow regression to the mean, not biology'}, ...
        'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'tex');
end

function print_gating_bias(groups, config)
    fprintf('\n=== 2AP GATING SELECTION-BIAS DIAGNOSTIC (PPR = A2/A1 of subset avg trace, A2 decay-sub) ===\n');
    fprintf('(Per-coverslip avg traces; no-gate == noGate QC reference. gate-A1 < no-gate < gate-A2 = bias.)\n');
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

%% ======================================================================
%% PPR estimator diagnostic (gate-either population, five A1/A2 definitions)
%% ======================================================================
function [names, labels] = ppr_diag_methods()
    % A1/A2 definitions compared (fields of twoap_pulse_amplitudes)
    names  = {'peak_max', 'win_mean', 'sync', 'p1lock', 'decay_sub'};
    labels = {'Peak-max', 'Window mean', 'Sync time-lock', 'P1-latency lock', 'Decay-subtracted'};
end

function [names, labels] = ppr_diag_estimators()
    % How A1/A2 become one PPR per coverslip x ISI. The first five use the
    % gate-either ROIs (P1 responder OR incremental P2 >= upper_sigma*sd); the
    % trace_all / trace_A1 / trace_A2 references show where the ungated and the
    % single-pulse-gated populations sit on the same amplitude definition.
    names  = {'trace_either', 'pooled', 'mean_ratio', 'median_ratio', 'geo_mean', ...
              'trace_all', 'trace_A1', 'trace_A2'};
    labels = {'avg trace (either)', 'pooled \Sigma A2 / \Sigma A1', 'mean of ROI ratios', ...
              'median of ROI ratios', 'geometric mean', ...
              'ref: avg trace, no gate', 'ref: avg trace, gate A1', 'ref: avg trace, gate A2'};
end

function D = compute_ppr_diagnostic(gd, frs, config)
    % PPR under five A1/A2 definitions x eight estimators, per coverslip x ISI.
    % Per-ROI amplitudes are measured on every QC-passed ROI with the shared
    % twoap_pulse_amplitudes (unclipped: negative incremental A2 is kept so no
    % estimator is pushed up by clipping). Gates: A1 = analyzer is_responder1;
    % A2 = incremental peak2_amp >= upper_sigma*sd (the residual-relative gate,
    % since the absolute P2 gate always passes on pulse 1's decay at short ISI).
    % If the per-ROI estimators agree with the gate-either average-trace ratio,
    % per-ROI PPR on that population is not dominated by divide-by-small-A1.
    methods = ppr_diag_methods();
    ests    = ppr_diag_estimators();
    nM      = numel(methods);
    [T, tisi] = pulse1_template(gd, config);
    D = struct('template', T, 'template_isi', tisi);

    % --- per-ROI amplitudes + gates (QC-passed ROIs) ---
    A1 = []; A2 = []; A2abs = []; ge = false(1, 0); g1 = false(1, 0); g2 = false(1, 0);
    cs = {}; isi = []; fidx = []; ridx = [];
    dsub = strcmp(methods, 'decay_sub');
    % Full-length per-ROI vectors in gd's ROI order (every ROI of every file,
    % NaN / false for ROIs that failed noise QC), for the per-ROI CSV.
    D.all = struct('a2_decay_sub', [], 'gate_A2', false(1, 0), 'gate_either', false(1, 0));
    for i = 1:numel(frs)
        fr = frs{i};
        if isempty(fr), continue; end
        t   = fr.twoap;
        qc  = find(isfinite(t.upper_threshold));
        fidx = [fidx, repmat(i, 1, numel(qc))];             %#ok<AGROW>
        ridx = [ridx, qc];                                  %#ok<AGROW>
        thr = t.upper_threshold - t.prestim_base;          % upper_sigma*sd per ROI
        a1f = nan(nM, numel(qc)); a2f = nan(nM, numel(qc)); a2af = nan(1, numel(qc));
        for j = 1:numel(qc)
            roi = qc(j);
            amp = twoap_pulse_amplitudes(fr.dfof_data(:, roi), t.stim1_frame, t.stim2_frame, ...
                t.prestim_base(roi), config, T);
            for m = 1:nM
                a1f(m, j) = amp.(methods{m}).a1;
                a2f(m, j) = amp.(methods{m}).a2;
            end
            a2af(j) = amp.a2_abs;          % pulse-2 peak above the pre-stim1 baseline
        end
        r1 = t.is_responder1(qc);
        r2 = isfinite(t.peak2_amp(qc)) & (t.peak2_amp(qc) >= thr(qc));
        A1  = [A1, a1f];                                    %#ok<AGROW>
        A2  = [A2, a2f];                                    %#ok<AGROW>
        A2abs = [A2abs, a2af];                              %#ok<AGROW>
        g1  = [g1, r1];                                     %#ok<AGROW>
        g2  = [g2, r2];                                     %#ok<AGROW>
        ge  = [ge, r1 | r2];                                %#ok<AGROW>
        nR  = numel(t.upper_threshold);
        fa  = nan(1, nR);   fa(qc)  = a2f(dsub, :);
        fg2 = false(1, nR); fg2(qc) = r2;
        fge = false(1, nR); fge(qc) = r1 | r2;
        D.all.a2_decay_sub = [D.all.a2_decay_sub, fa];
        D.all.gate_A2      = [D.all.gate_A2, fg2];
        D.all.gate_either  = [D.all.gate_either, fge];
        cs  = [cs, repmat({parse_coverslip_id(fr.filename)}, 1, numel(qc))]; %#ok<AGROW>
        isi = [isi, repmat(fr.isi_ms, 1, numel(qc))];       %#ok<AGROW>
    end

    % --- per coverslip x ISI values ---
    keys  = cellfun(@(c, x) sprintf('%s|%g', c, x), cs, num2cell(isi), 'UniformOutput', false);
    uk    = unique(keys(ge), 'stable');
    kISI  = cellfun(@(k) str2double(extractAfter(k, '|')), uk);
    ca    = gd.cs_avg_gate;
    caKey = arrayfun(@(e) sprintf('%s|%g', e.cs, e.isi), ca, 'UniformOutput', false);
    stim1 = config.twoAP.stim1_frame;
    trSrc = struct('trace_either', 'gate_either', 'trace_all', 'all', ...
                   'trace_A1', 'gate_A1', 'trace_A2', 'gate_A2');
    for m = 1:nM
        mm = methods{m};
        V  = struct();
        for e = 1:numel(ests), V.(ests{e}) = nan(1, numel(uk)); end
        for k = 1:numel(uk)
            sel = ge & strcmp(keys, uk{k});
            a1 = A1(m, sel); a2 = A2(m, sel);
            ok = isfinite(a1) & isfinite(a2);
            if any(ok), V.pooled(k) = sum(a2(ok)) / sum(a1(ok)); end
            pos = ok & (a1 > 0);
            if any(pos)
                r = a2(pos) ./ a1(pos);
                V.mean_ratio(k)   = mean(r);
                V.median_ratio(k) = median(r);
            end
            lg = pos & (a2 > 0);
            if any(lg), V.geo_mean(k) = exp(mean(log(a2(lg) ./ a1(lg)))); end
            ci = find(strcmp(caKey, uk{k}), 1);
            if ~isempty(ci)
                tf = fieldnames(trSrc);
                for q = 1:numel(tf)
                    if ca(ci).n.(trSrc.(tf{q})) < 1, continue; end
                    [~, ~, ~, ~, ~, amp] = trace_ppr(ca(ci).traces.(trSrc.(tf{q})), stim1, ...
                        ca(ci).stim2, config, ca(ci).pre, T);
                    V.(tf{q})(k) = amp.(mm).a2 / amp.(mm).a1;
                end
            end
        end
        for e = 1:numel(ests)
            D.(mm).(ests{e}) = summarize_by_isi(kISI, V.(ests{e}));
        end
        if strcmp(mm, 'decay_sub'), Vd = V; end
    end

    % --- per coverslip x ISI rows of the PRIMARY PPR (for the coverslip CSV) ---
    nK = numel(uk);
    D.cs_rows = struct('coverslip', {cellfun(@(k) extractBefore(k, '|'), uk, 'UniformOutput', false)}, ...
        'isi_ms', kISI, 'n_either', zeros(1, nK), 'n_A1only', zeros(1, nK), 'n_A2only', zeros(1, nK), ...
        'sum_a1', nan(1, nK), 'sum_a2_decay_sub', nan(1, nK), 'mean_a1', nan(1, nK), ...
        'mean_a2_decay_sub', nan(1, nK), 'mean_a2_abs', nan(1, nK), 'ppr_primary', nan(1, nK), ...
        'n_excluded', zeros(1, nK), 'ppr_sub_threshold', nan(1, nK), 'ppr_ref_no_gate', nan(1, nK));
    for k = 1:nK
        inK = strcmp(keys, uk{k});
        sel = ge & inK;
        a1 = A1(dsub, sel); a2 = A2(dsub, sel); ok = isfinite(a1) & isfinite(a2);
        D.cs_rows.n_either(k)  = sum(sel);
        D.cs_rows.n_A1only(k)  = sum(inK & g1 & ~g2);
        D.cs_rows.n_A2only(k)  = sum(inK & g2 & ~g1);
        D.cs_rows.sum_a1(k)    = sum(a1(ok));
        D.cs_rows.sum_a2_decay_sub(k) = sum(a2(ok));
        D.cs_rows.mean_a1(k)   = mean(a1(ok));             % mean_a2 / mean_a1 == ppr_primary
        D.cs_rows.mean_a2_decay_sub(k) = mean(a2(ok));
        a2a = A2abs(sel);
        D.cs_rows.mean_a2_abs(k) = mean(a2a(isfinite(a2a)));
        D.cs_rows.ppr_primary(k)     = Vd.pooled(k);
        D.cs_rows.ppr_ref_no_gate(k) = Vd.trace_all(k);
        D.cs_rows.n_excluded(k)      = sum(inK & ~ge);
        % Sub-threshold population: QC-passed ROIs that fail BOTH gates. Single
        % ROIs are mostly ~3 sigma (too noisy for per-ROI peaks, which compress the
        % ratio toward 1), so their PPR is taken from their AVERAGE trace.
        ci = find(strcmp(caKey, uk{k}), 1);
        if ~isempty(ci) && ca(ci).n.excluded >= 1
            [~, ~, ~, ~, ~, amp] = trace_ppr(ca(ci).traces.excluded, stim1, ca(ci).stim2, ...
                config, ca(ci).pre, T);
            D.cs_rows.ppr_sub_threshold(k) = amp.decay_sub.a2 / amp.decay_sub.a1;
        end
    end
    D.sub_threshold = summarize_by_isi(kISI, D.cs_rows.ppr_sub_threshold);
    D.amp_a1        = summarize_by_isi(kISI, D.cs_rows.mean_a1);
    D.amp_a2        = summarize_by_isi(kISI, D.cs_rows.mean_a2_decay_sub);
    D.amp_a2_abs    = summarize_by_isi(kISI, D.cs_rows.mean_a2_abs);

    % --- per-ROI ratios + counts for the distribution figure / table ---
    uISI = unique(isi(ge));
    D.isi = uISI;
    D.roi = struct('isi', isi(ge), 'cs', {cs(ge)});
    D.n_either = arrayfun(@(x) sum(ge & isi == x), uISI);
    D.n_A1only = arrayfun(@(x) sum(g1 & ~g2 & isi == x), uISI);
    D.n_A2only = arrayfun(@(x) sum(g2 & ~g1 & isi == x), uISI);
    for m = 1:nM
        a1 = A1(m, ge); a2 = A2(m, ge);
        r = a2 ./ a1; r(~(a1 > 0)) = NaN;
        D.roi.ratio.(methods{m}) = r;
        D.roi.a1.(methods{m})    = a1;
        D.roi.a2.(methods{m})    = a2;
        D.n_a1_nonpos.(methods{m}) = arrayfun(@(x) sum(~(a1 > 0) & D.roi.isi == x), uISI);
    end

    % --- representative gate-either ROIs for the decay-subtraction figure ---
    % Per short ISI (where pulse 1 still decays under pulse 2), the gate-either
    % ROIs ranked by pulse-1 amplitude, sampled at fixed quantiles so the tiles
    % span weak to strong ROIs (deterministic, not hand-picked).
    pm = strcmp(methods, 'peak_max'); dm = strcmp(methods, 'decay_sub');
    qs = [0.2 0.35 0.5 0.65 0.8 0.95];
    exISI = uISI(uISI <= 120);
    D.examples = struct('isi', {}, 'avg', {}, 'rois', {});
    for k = 1:numel(exISI)
        cand = find(ge & isi == exISI(k) & A1(pm, :) > 0 & isfinite(A2(dm, :)));
        rois = struct('trace', {}, 'base', {}, 'stim1', {}, 'stim2', {}, 'label', {}, ...
            'gate', {}, 'amp', {});
        if ~isempty(cand)
            [~, o] = sort(A1(pm, cand)); cand = cand(o);
            pick = unique(cand(max(1, round(qs * numel(cand)))), 'stable');
            for j = pick
                fr = frs{fidx(j)}; t = fr.twoap; roi = ridx(j);
                tr = double(fr.dfof_data(:, roi));
                amp = twoap_pulse_amplitudes(tr, t.stim1_frame, t.stim2_frame, ...
                    t.prestim_base(roi), config, T);
                gnames = {'A1', 'A2'};
                gtxt = strjoin(gnames([g1(j), g2(j)]), '+');
                rois(end+1) = struct('trace', tr, 'base', t.prestim_base(roi), ... %#ok<AGROW>
                    'stim1', t.stim1_frame, 'stim2', t.stim2_frame, ...
                    'label', sprintf('%s-%s ROI%d', parse_coverslip_id(fr.filename), ...
                        parse_region_id(fr.filename), roi), ...
                    'gate', gtxt, 'amp', amp);
            end
        end
        % Gate-either coverslip average traces at this ISI (one per coverslip)
        avg = struct('trace', {}, 'base', {}, 'stim1', {}, 'stim2', {}, 'label', {}, ...
            'gate', {}, 'amp', {});
        for ci = find([ca.isi] == exISI(k))
            if ca(ci).n.gate_either < 1, continue; end
            tr = double(ca(ci).traces.gate_either(:));
            pv = tr(ca(ci).pre); bs = mean(pv(isfinite(pv)));
            amp = twoap_pulse_amplitudes(tr, stim1, ca(ci).stim2, bs, config, T);
            avg(end+1) = struct('trace', tr, 'base', bs, 'stim1', stim1, ... %#ok<AGROW>
                'stim2', ca(ci).stim2, 'label', sprintf('%s avg (n=%d)', ca(ci).cs, ca(ci).n.gate_either), ...
                'gate', 'either', 'amp', amp);
        end
        D.examples(end+1) = struct('isi', exISI(k), 'avg', {avg}, 'rois', {rois});
    end
end

function [T, tisi] = pulse1_template(gd, config)
    % Pulse-1 waveform for the decay-subtracted amplitude: the mean of the group's
    % longest-ISI coverslip average traces (all QC-passed ROIs) from stim1 up to
    % stim2, baseline-subtracted and normalised to its pulse-1 peak. At the longest
    % ISI pulse 1 decays undisturbed, so this is the tail pulse 2 sits on at the
    % short ISIs. Once it first falls below decay_template_min_frac it is set to
    % 0 (pulse 1 has decayed; the rest is baseline noise/drift, not tail). Empty
    % when the longest ISI is under 250 ms (pulse 1 has not decayed before
    % pulse 2 there) or the group has no usable average trace.
    T = []; tisi = NaN;
    if ~isfield(gd, 'cs_avg') || isempty(gd.cs_avg), return; end
    ca    = gd.cs_avg;
    tisi  = max([ca.isi]);
    if tisi < 250, return; end
    sel   = find([ca.isi] == tisi);
    stim1 = config.twoAP.stim1_frame;
    len   = min(arrayfun(@(e) e.stim2 - stim1, ca(sel)));
    Y = nan(len, numel(sel));
    for j = 1:numel(sel)
        e  = ca(sel(j)); tr = double(e.trace(:));
        pv = tr(e.pre); base = mean(pv(isfinite(pv)));
        Y(:, j) = tr(stim1:stim1 + len - 1) - base;
    end
    y = mean(Y, 2, 'omitnan');
    search = max(1, round(config.twoAP.peak_search_ms / config.frame_dur_ms));
    [pk, ip] = max(y(1:min(search + 1, len)));
    if ~(pk > 0), return; end
    T = y / pk;
    kz = find(T(ip:end) < config.twoAP.decay_template_min_frac, 1);
    if ~isempty(kz), T(ip + kz - 1:end) = 0; end
end

function S = summarize_by_isi(isiE, vals)
    % Per-ISI mean +/- SEM across coverslips of one value per coverslip x ISI.
    uISI = unique(isiE);
    S = struct('isi', uISI, 'cs_vals', {cell(1, numel(uISI))}, 'mean', nan(1, numel(uISI)), ...
        'sem', nan(1, numel(uISI)), 'n_cs', zeros(1, numel(uISI)));
    for k = 1:numel(uISI)
        v = vals(isiE == uISI(k)); v = v(isfinite(v));
        S.cs_vals{k} = v; S.n_cs(k) = numel(v);
        if ~isempty(v), S.mean(k) = mean(v); S.sem(k) = std(v) / sqrt(numel(v)); end
    end
end

function figs = plot_twoap_ppr_diagnostic(groups, config)
    % One figure per group, one panel per A1/A2 definition: PPR vs ISI from the
    % five gate-either estimators (solid) and the no-gate / gate-A1 / gate-A2
    % average-trace references (dashed). Values are per-coverslip means; the
    % y-axis is shared across panels and capped so a blown-up mean-of-ratios
    % does not flatten the rest (off-scale values are in the printed table).
    [methods, mlabels] = ppr_diag_methods();
    [ests, elabels]    = ppr_diag_estimators();
    cols = [0 0 0; 0.10 0.60 0.30; 0.85 0.45 0.05; 0.55 0.25 0.75; 0.00 0.60 0.75; ...
            0.55 0.55 0.55; 0.20 0.45 0.95; 0.90 0.25 0.20];
    styl = {'-', '-', '-', '-', '-', '--', ':', ':'};
    lw   = [2.6 1.6 1.4 1.4 1.4 1.4 1.4 1.4];
    ycap = 4;
    figs = gobjects(1, 0);
    for g = 1:numel(groups)
        gd = groups(g);
        if ~isfield(gd, 'ppr_diag') || isempty(gd.ppr_diag), continue; end
        D = gd.ppr_diag;
        fig = figure('Name', sprintf('2AP PPR estimator diagnostic - %s', gd.name), ...
            'Position', [60 100 1900 520], 'Color', 'w');
        figs(end+1) = fig; %#ok<AGROW>
        tl = tiledlayout(fig, 1, numel(methods), 'TileSpacing', 'compact', 'Padding', 'compact');
        axs = gobjects(1, numel(methods)); ymax = 1;
        for m = 1:numel(methods)
            ax = nexttile(tl); hold(ax, 'on'); axs(m) = ax;
            for e = 1:numel(ests)
                S = D.(methods{m}).(ests{e}); v = isfinite(S.mean);
                if ~any(v), continue; end
                if e == 1
                    errorbar(ax, S.isi(v), S.mean(v), S.sem(v), '-o', 'Color', cols(e, :), ...
                        'MarkerFaceColor', cols(e, :), 'LineWidth', lw(e), 'MarkerSize', 5, ...
                        'CapSize', 3, 'DisplayName', elabels{e});
                else
                    plot(ax, S.isi(v), S.mean(v), styl{e}, 'Color', cols(e, :), 'LineWidth', lw(e), ...
                        'Marker', '.', 'MarkerSize', 12, 'DisplayName', elabels{e});
                end
                ymax = max(ymax, max(S.mean(v)));
            end
            set(ax, 'XScale', 'log'); grid(ax, 'on'); box(ax, 'on');
            yline(ax, 1, 'k--', 'HandleVisibility', 'off');
            xlabel(ax, 'Inter-stimulus interval (ms)');
            if m == 1, ylabel(ax, 'PPR (A2 / A1)'); end
            title(ax, mlabels{m});
            if m == numel(methods), legend(ax, 'Location', 'eastoutside', 'FontSize', 8); end
        end
        for a = axs, ylim(a, [0 min(ymax * 1.05, ycap)]); end
        title(tl, {sprintf(['2AP PPR estimator diagnostic \\bullet %s \\bullet gate-either ROIs ' ...
            '(P1 responder OR incremental A2 \\geq %g\\sigma)'], gd.name, config.twoAP.upper_threshold_sigma), ...
            sprintf(['per-coverslip mean (SEM on the avg-trace curve); dashed/dotted = avg-trace references; ' ...
            'decay template = %g ms ISI avg trace; y capped at %g'], D.template_isi, ycap)}, ...
            'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'tex');
    end
end

function figs = plot_twoap_decay_examples(groups, config)
    % Decay subtraction on representative gate-either ROIs, one figure per group:
    % rows = short ISIs (<= 120 ms), first column(s) = the gate-either coverslip
    % average trace, then ROIs sampled at fixed pulse-1 amplitude quantiles
    % (20th..95th percentile). Each tile marks the peaks the PPR uses.
    figs  = gobjects(1, 0);
    tax   = make_time_axis(config); xunit = tax.unit;
    for g = 1:numel(groups)
        gd = groups(g);
        if ~isfield(gd, 'ppr_diag') || ~isfield(gd.ppr_diag, 'examples') ...
                || isempty(gd.ppr_diag.examples) || isempty(gd.ppr_diag.template)
            continue;
        end
        D  = gd.ppr_diag; ex = D.examples;
        nA = max(arrayfun(@(e) numel(e.avg), ex));
        nC = nA + max(arrayfun(@(e) numel(e.rois), ex));
        fig = figure('Name', sprintf('2AP decay-subtraction examples - %s', gd.name), ...
            'Position', [30 30 max(1400, 270 * nC) max(600, 250 * numel(ex))], 'Color', 'w');
        figs(end+1) = fig; %#ok<AGROW>
        tl = tiledlayout(fig, numel(ex), nC, 'TileSpacing', 'compact', 'Padding', 'compact');
        for r = 1:numel(ex)
            for c = 1:numel(ex(r).avg)
                draw_decay_tile(nexttile(tl, (r - 1) * nC + c), ex(r).avg(c), D.template, config, ex(r).isi);
            end
            for c = 1:numel(ex(r).rois)
                draw_decay_tile(nexttile(tl, (r - 1) * nC + nA + c), ex(r).rois(c), D.template, config, ex(r).isi);
            end
        end
        title(tl, {sprintf(['2AP decay subtraction \\bullet %s \\bullet gate-either (left: coverslip avg trace; ' ...
            'then ROIs at the 20/35/50/65/80/95th pulse-1 amplitude percentile)'], gd.name), ...
            sprintf(['black = dF/F; orange dashed = projected pulse-1 tail (%g ms ISI template scaled to the pre-stim2 ' ...
            'level); orange solid = pulse-2 window with the tail subtracted; shaded = peak windows / pre-stim2 level window'], ...
            D.template_isi), ...
            ['o blue = A1 peak; o red = peak-max A2 (above flat pre-stim2 level); ' ...
            'diamond orange = decay-subtracted A2 peak; x in ' ...
            xunit]}, ...
            'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'tex');
    end
end

function draw_decay_tile(ax, e, T, config, isi)
    % One trace with the pulse-1 tail projection, the decay-subtracted pulse-2
    % window and the A1 / A2 peaks marked (see plot_twoap_decay_examples).
    hold(ax, 'on');
    fd = config.frame_dur_ms; s1 = e.stim1; s2 = e.stim2; tr = e.trace; nF = numel(tr);
    search = max(1, round(config.twoAP.peak_search_ms / fd));
    f0 = max(s1 - round(50 / fd), 1); f1 = min(s2 + round(150 / fd), nF);
    tax = make_time_axis(config); cv = tax.from_ms;   % same x unit as the montage
    tt = @(f) cv((f - s1) * fd);
    [r_lo, r_hi] = twoap_residual_window(s1, s2, config.twoAP.residual_window_ms, fd);
    w1e = s1 + min(search, max(s2 - s1 - 1, 1));
    w2  = s2:min(s2 + search, nF);
    orange = [0.95 0.50 0.05];
    xregion(ax, tt(s1), tt(w1e), 'FaceColor', [0.2 0.4 1.0], 'FaceAlpha', 0.07);
    xregion(ax, tt(w2(1)), tt(w2(end)), 'FaceColor', [1.0 0.3 0.3], 'FaceAlpha', 0.07);
    xregion(ax, tt(r_lo) - cv(fd / 2), tt(r_hi) + cv(fd / 2), 'FaceColor', [0.4 0.4 0.4], 'FaceAlpha', 0.15);
    plot(ax, tt(f0:f1), tr(f0:f1), 'k-', 'LineWidth', 1);
    yline(ax, e.base, ':', 'Color', [0.4 0.4 0.4]);
    xline(ax, 0, '--b'); xline(ax, tt(s2), '--r');

    a = e.amp;
    if isfinite(a.decay_scale)
        idx  = r_lo:w2(end); idx = idx(idx - s1 + 1 <= numel(T));
        plot(ax, tt(idx), e.base + a.decay_scale * T(idx - s1 + 1), '--', 'Color', orange, 'LineWidth', 1.3);
        sub = tr(w2) - a.decay_scale * T(w2 - s1 + 1);
    else
        plot(ax, tt([r_lo, w2(end)]), [a.resid a.resid], '--', 'Color', orange, 'LineWidth', 1.3);
        sub = tr(w2) - (a.resid - e.base);
    end
    plot(ax, tt(w2), sub, '-', 'Color', orange, 'LineWidth', 1.8);

    % A1: pulse-1 peak
    if isfinite(a.t1)
        f = s1 + a.t1;
        plot(ax, tt(f), tr(f), 'o', 'Color', [0.2 0.4 1.0], 'MarkerFaceColor', [0.2 0.4 1.0], 'MarkerSize', 5);
    end
    % Peak-max A2: raw pulse-2 peak
    if isfinite(a.t2)
        f = s2 + a.t2;
        plot(ax, tt(f), tr(f), 'o', 'Color', [0.85 0.15 0.15], 'MarkerFaceColor', [0.85 0.15 0.15], 'MarkerSize', 5);
    end
    % Decay-subtracted A2: peak of the subtracted window
    if isfinite(a.t2_decay)
        f  = s2 + a.t2_decay;
        yv = e.base + a.decay_sub.a2;
        plot(ax, tt(f), yv, 'd', 'Color', orange, 'MarkerFaceColor', orange, 'MarkerSize', 6);
    end
    xlim(ax, [tt(f0) tt(f1)]);
    grid(ax, 'on'); box(ax, 'on'); set(ax, 'FontSize', 7);
    ax.YAxis.Exponent = 0; ytickformat(ax, '%.3g');
    title(ax, sprintf('ISI %g | %s [%s]\nPPR peak-max %.2f | decay-sub %.2f', isi, e.label, e.gate, ...
        a.peak_max.a2 / a.peak_max.a1, a.decay_sub.a2 / a.decay_sub.a1), 'FontSize', 7, 'Interpreter', 'none');
end

function figs = plot_twoap_ppr_distributions(groups, ~)
    % Per-ROI PPR distribution of the gate-either ROIs, by ISI, one panel per
    % A1/A2 definition (ROIs pooled across coverslips: distribution shape only,
    % not a statistic). Box = median/IQR; diamond = mean of ratios; x = pooled
    % SigmaA2/SigmaA1. A mean far from the median and pooled values = tails.
    [methods, mlabels] = ppr_diag_methods();
    figs = gobjects(1, 0);
    for g = 1:numel(groups)
        gd = groups(g);
        if ~isfield(gd, 'ppr_diag') || isempty(gd.ppr_diag), continue; end
        D = gd.ppr_diag;
        fig = figure('Name', sprintf('2AP per-ROI PPR distributions - %s', gd.name), ...
            'Position', [70 110 1900 480], 'Color', 'w');
        figs(end+1) = fig; %#ok<AGROW>
        tl = tiledlayout(fig, 1, numel(methods), 'TileSpacing', 'compact', 'Padding', 'compact');
        isiLab = compose('%g', D.isi);
        for m = 1:numel(methods)
            ax = nexttile(tl); hold(ax, 'on');
            r  = D.roi.ratio.(methods{m});
            a1 = D.roi.a1.(methods{m}); a2 = D.roi.a2.(methods{m});
            x  = categorical(compose('%g', D.roi.isi), isiLab, 'Ordinal', true);
            v  = isfinite(r);
            boxchart(ax, x(v), r(v), 'BoxFaceColor', gd.color, 'MarkerStyle', 'none');
            for k = 1:numel(D.isi)
                s = D.roi.isi == D.isi(k);
                plot(ax, categorical(isiLab(k), isiLab, 'Ordinal', true), mean(r(s & v)), 'kd', ...
                    'MarkerFaceColor', 'k', 'MarkerSize', 5);
                ok = s & isfinite(a1) & isfinite(a2);
                plot(ax, categorical(isiLab(k), isiLab, 'Ordinal', true), sum(a2(ok)) / sum(a1(ok)), ...
                    'rx', 'MarkerSize', 8, 'LineWidth', 1.5);
            end
            yline(ax, 1, 'k--');
            ylim(ax, [-1 4]); grid(ax, 'on'); box(ax, 'on');
            xlabel(ax, 'ISI (ms)');
            if m == 1, ylabel(ax, 'per-ROI A2 / A1'); end
            title(ax, sprintf('%s (A1\\leq0 dropped: %d)', mlabels{m}, sum(D.n_a1_nonpos.(methods{m}))));
        end
        title(tl, {sprintf('2AP per-ROI PPR distributions \\bullet %s \\bullet gate-either ROIs (pooled across coverslips)', gd.name), ...
            'box = median / IQR (whiskers 1.5 IQR, outliers hidden); black diamond = mean of ROI ratios; x = pooled \Sigma A2 / \Sigma A1; y clipped to [-1, 4]'}, ...
            'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'tex');
    end
end

function print_ppr_diagnostic(groups, config)
    [methods, mlabels] = ppr_diag_methods();
    ests = ppr_diag_estimators();
    fprintf('\n=== 2AP PPR ESTIMATOR DIAGNOSTIC (gate-either: P1 responder OR incremental A2 >= %g sigma) ===\n', ...
        config.twoAP.upper_threshold_sigma);
    fprintf('Per-coverslip mean PPR. trEith = avg trace of gate-either ROIs; pooled = SumA2/SumA1;\n');
    fprintf('mean/median/geo = of per-ROI A2/A1 (A1<=0 dropped; geo also drops A2<=0); refs = avg trace.\n');
    for g = 1:numel(groups)
        gd = groups(g);
        if ~isfield(gd, 'ppr_diag') || isempty(gd.ppr_diag), continue; end
        D = gd.ppr_diag;
        fprintf('\n%s  (decay template: %g ms ISI)\n', gd.name, D.template_isi);
        fprintf('  ISI  n_either (A1only/A2only)\n');
        for k = 1:numel(D.isi)
            fprintf('  %5g  %5d (%d/%d)\n', D.isi(k), D.n_either(k), D.n_A1only(k), D.n_A2only(k));
        end
        for m = 1:numel(methods)
            mm = methods{m};
            fprintf('  -- %s --\n', mlabels{m});
            fprintf('    ISI  trEith pooled   mean median    geo |  noGate gateA1 gateA2 | A1<=0\n');
            for k = 1:numel(D.isi)
                x = D.isi(k);
                v = cellfun(@(e) val_at(D.(mm).(e), x), ests);
                fprintf('  %5g  %6.2f %6.2f %6.2f %6.2f %6.2f | %6.2f %6.2f %6.2f | %4d\n', ...
                    x, v, D.n_a1_nonpos.(mm)(k));
            end
        end
    end
    fprintf('==========================================================================================\n');
end

function v = val_at(S, isi)
    v = NaN;
    k = find(S.isi == isi, 1);
    if ~isempty(k), v = S.mean(k); end
end

function v = getfield_default(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function figs = plot_twoap_montage(groups, config)
    % Per group, example gate-either ROIs (the primary PPR's population: pulse-1
    % responder OR incremental A2 >= upper_sigma*sd) from the shortest-ISI
    % reference file, with each amplitude method's picks marked so they can be
    % compared by eye. Tile titles give the ROI's A1 and decay-subtracted A2,
    % the amplitudes the primary pools (markers as in the decay tiles).
    figs = [];
    orange    = [0.95 0.50 0.05];
    frame_dur = config.frame_dur_ms;
    stim1     = config.twoAP.stim1_frame;
    sync_f    = max(1, round(config.twoAP.sync_latency_ms / frame_dur));
    mc        = config.visualization.montage;
    maxROIs   = mc.twoAP_max_rois;
    tax       = make_time_axis(config);
    cv        = tax.from_ms;   % ms position -> display units (for plotting)

    for g = 1:numel(groups)
        gd  = groups(g);
        ref = gd.montage_ref;
        if isempty(ref.dfof) || isempty(ref.twoap), continue; end
        t2    = ref.twoap;
        stim2 = t2.stim2_frame;
        nCol  = size(ref.dfof, 2);
        thr   = t2.upper_threshold - t2.prestim_base;               % upper_sigma*sd
        ge    = isfinite(t2.upper_threshold) & (t2.is_responder1 | ...
            (isfinite(t2.peak2_amp) & t2.peak2_amp >= thr));
        resp  = find(ge); resp = resp(resp <= nCol);
        if isempty(resp), continue; end
        T     = gd.ppr_diag.template;
        order = resp(1:min(maxROIs, numel(resp)));

        pre_f  = round(config.visualization.montage.twoAP_pre_ms  / frame_dur);
        post_f = round(config.visualization.montage.twoAP_post_ms / frame_dur);
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
        [q_lo, q_hi] = twoap_charge_window(stim1, stim2, nF_full, config.twoAP.auc_post_ms, frame_dur);

        % With montage.rows set, ROIs beyond one rows x columns grid continue
        % on further figures (pages).
        [nr, nc, per, nP] = montage_grid(numel(order), mc);
        for pg = 1:nP
            pidx = order((pg-1)*per + 1 : min(pg*per, numel(order)));
            if nP > 1, ptxt = sprintf(' (page %d/%d)', pg, nP); else, ptxt = ''; end
            fig = figure('Name', sprintf('2AP Montage %s (ISI %g ms)%s', gd.name, ref.isi, ptxt), ...
                'Position', [50 50 1850 1000], 'Color', 'w');
            figs = [figs, fig]; %#ok<AGROW>
            tl = tiledlayout(fig, nr, nc, 'TileSpacing', 'tight', 'Padding', 'tight');

            for k = 1:numel(pidx)
                roi = pidx(k);
                nexttile(tl, k); hold on;
                % Grey band = the frames Q integrates; blue / red = peak-search windows
                xregion(cv((q_lo - stim1) * frame_dur), cv((q_hi - stim1) * frame_dur), ...
                    'FaceColor', [0.45 0.45 0.50], 'FaceAlpha', 0.12);
                xregion(0, cv(s1_end_ms), 'FaceColor', [0.2 0.4 1.0], 'FaceAlpha', 0.10);
                xregion(cv(stim2_ms), cv(stim2_ms + s2_end_ms), 'FaceColor', [1.0 0.3 0.3], 'FaceAlpha', 0.10);
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
                setappdata(gca, 'ymin', config.visualization.ymin_montage);   % see apply_ymin_zero
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
                    if f <= nF, plot(cv(stim2_ms + t1ms), ref.dfof(f, roi), 'r^', 'MarkerSize', 6, 'LineWidth', 1.2); end
                end
                % Decay-subtracted A2 (primary amplitude): orange diamond at
                % baseline + A2, the height pulse 2 adds once pulse 1's tail is removed
                amp = twoap_pulse_amplitudes(ref.dfof(:, roi), stim1, stim2, t2.prestim_base(roi), config, T);
                if isfinite(amp.t2_decay)
                    plot(cv(stim2_ms + amp.t2_decay * frame_dur), t2.prestim_base(roi) + amp.decay_sub.a2, ...
                        'd', 'Color', orange, 'MarkerFaceColor', orange, 'MarkerSize', 6);
                end
                % Title gives the amplitudes the primary PPR pools (SumA2/SumA1),
                % not a per-ROI ratio: A2-only ROIs have a tiny A1 that blows it up.
                gtxt = {'A2', 'A1', 'A1+A2'};
                gi   = 2 * t2.is_responder1(roi) + (t2.peak2_amp(roi) >= thr(roi));
                ttl = sprintf('ROI%d [%s]  A1=%.3f A2=%.3f Q=%.1f', roi, gtxt{gi}, ...
                    amp.decay_sub.a1, amp.decay_sub.a2, t2.q_pair(roi));
                title(ttl, 'FontSize', 8);
                set(gca, 'FontSize', 7); grid on; box on;
                scroll_time_axis(gca, xwin);   % show stim window; drag to scroll x
            end
            title(tl, sprintf(['2AP %s, ISI=%g ms%s \\bullet example gate-either ROIs [gate passed] (one reference file); ' ...
                'tile A1 / A2 (dF/F, A2 decay-subtracted) = what the primary PPR pools as \\SigmaA2/\\SigmaA1; ' ...
                'Q = total charge (dF/F*ms) \\bullet drag to scroll x\n' ...
                'blue / red dashed = stim1 / stim2; grey band = Q window; blue / red bands = peak-search windows; ' ...
                'o = peak-max, x = sync %gms, ' ...
                '\\Delta = p1-lock, orange diamond = decay-subtracted A2 (baseline + A2); ' ...
                'green = upper %g\\sigma / lower %g\\sigma; x in %s'], gd.name, ref.isi, ptxt, ...
                config.twoAP.sync_latency_ms, config.twoAP.upper_threshold_sigma, config.twoAP.lower_threshold_sigma, tax.unit), ...
                'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'tex');
        end
    end
end

function figs = plot_twoap_cs_average(groups, config)
    % Coverslip x ISI grand-average dF/F (QC). One tile per coverslip (rows) and
    % ISI (cols); the trace is the mean over all QC-passed ROIs (no response
    % gate). Each title gives the coverslip's PRIMARY PPR (gate-either pooled,
    % as in the PPR figure panel 1), this trace's own decay-subtracted A2/A1
    % (panel 3) and its total charge Q (panel 2).
    %
    % Split into TWO figures by ISI so the short- and long-interval panels each
    % get readable column widths: short ISI (<= 120 ms) and long ISI (> 120 ms).

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

    split = 120;   % ms: short ISIs (pulse 2 on pulse 1's decay) vs long
    bands = {isis(isis <= split), 'short ISI'; isis(isis > split), 'long ISI'};
    for b = 1:size(bands, 1)
        sub = bands{b,1};
        if isempty(sub), continue; end
        lbl = sprintf('%s (%s ms)', bands{b,2}, strjoin(compose('%g', sub), ', '));
        f = render_cs_average_fig(groups, config, rows, rowGroup, sub, lbl);
        if ~isempty(f) && isgraphics(f), figs(end+1) = f; end %#ok<AGROW>
    end
    % One y-scale for every tile in both figures (largest trace plotted)
    if ~isempty(figs)
        axs = arrayfun(@(f) getappdata(f, 'trace_axes'), figs, 'UniformOutput', false);
        share_ylim([axs{:}], config);
    end
end

function fig = render_cs_average_fig(groups, config, rows, rowGroup, isis, bandLabel)
    % Render ONE coverslip grand-average figure for the given ISI subset (cols).
    frame_dur = config.frame_dur_ms;
    stim1     = config.twoAP.stim1_frame;
    pre_f     = round(config.visualization.montage.twoAP_pre_ms  / frame_dur);
    post_f    = round(config.visualization.montage.twoAP_post_ms / frame_dur);
    search_f  = max(1, round(config.twoAP.peak_search_ms / frame_dur));
    tax       = make_time_axis(config);
    cv        = tax.from_ms;   % ms -> display units (same x as the montage)

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
            setappdata(axh, 'ymin', config.visualization.ymin_mean_trace);   % see apply_ymin_zero
            idx = find(strcmp({gd.cs_avg.cs}, csname) & [gd.cs_avg.isi] == isi, 1);
            if isempty(idx), axis off; continue; end
            rowAxes(r, cix) = axh;
            e = gd.cs_avg(idx);
            trace = e.trace; stim2 = e.stim2; nF = numel(trace);
            f0 = max(stim1 - pre_f, 1); f1 = min(stim2 + post_f, nF); fidx = f0:f1;
            t  = cv((fidx - stim1) * frame_dur); stim2_ms = (stim2 - stim1) * frame_dur;
            % Grey band = the total-charge (Q) window, the exact frames q_pair
            % integrates (twoap_charge_window, from config.twoAP.auc_post_ms).
            % Blue / red = the pulse-1 / pulse-2 peak-search windows, as in the
            % montage and decay tiles (pulse 1 capped before stim2).
            [q_lo, q_hi] = twoap_charge_window(stim1, stim2, nF, config.twoAP.auc_post_ms, frame_dur);
            xregion(cv((q_lo - stim1) * frame_dur), cv((q_hi - stim1) * frame_dur), ...
                'FaceColor', [0.45 0.45 0.50], 'FaceAlpha', 0.12);
            s1_end_ms = min(search_f, max(stim2 - stim1 - 1, 1)) * frame_dur;
            xregion(0, cv(s1_end_ms), 'FaceColor', [0.2 0.4 1.0], 'FaceAlpha', 0.10);
            xregion(cv(stim2_ms), cv(stim2_ms + config.twoAP.peak_search_ms), ...
                'FaceColor', [1.0 0.3 0.3], 'FaceAlpha', 0.10);
            % Individual region (c#) traces behind the mean show the within-
            % coverslip spread. No SEM band: ROIs are not independent samples
            % (the statistical unit is the coverslip).
            yy = trace(fidx);
            regs = e.regions; nrg = numel(regs);
            for rr = 1:nrg
                if numel(regs(rr).trace) < f1, continue; end
                plot(t, regs(rr).trace(fidx), '-', 'Color', [0.6 0.6 0.6], ...
                    'LineWidth', 0.5, 'HandleVisibility', 'off');
            end
            plot(t, yy, '-', 'Color', gd.color, 'LineWidth', 1.3, 'HandleVisibility', 'off');
            % xline spans the full axis height, so the stim markers stay full-height
            % after the y-axes are shared (share_ylim).
            xline(0, '--', 'Color', 'b', 'LineWidth', 0.8);
            xline(cv(stim2_ms), '--', 'Color', 'r', 'LineWidth', 0.8);

            [ppr_tr, q_pair] = trace_ppr_primary(trace, stim1, stim2, config, e.pre, gd.ppr_diag.template);
            cr = gd.ppr_diag.cs_rows;
            ki = find(strcmp(cr.coverslip, csname) & cr.isi_ms == isi, 1);
            ppr_p = NaN;
            if ~isempty(ki), ppr_p = cr.ppr_primary(ki); end
            title(sprintf('%s  ISI %g  PPR %.2f (trace %.2f)  Q %.1f  [%d c]', ...
                csname, isi, ppr_p, ppr_tr, q_pair, nrg), 'FontSize', 8);
            set(gca, 'FontSize', 7); grid on; box on;
            axh.YAxis.Exponent = 0;   % plain decimal y ticks (no 10^n exponent)
            ytickformat(axh, '%.3g');
            if cix == 1, ylabel(sprintf('%s / %s', gd.name, csname), 'FontSize', 8, 'FontWeight', 'bold'); end
            if r == nR, xlabel(sprintf('%s from stim1', tax.unit), 'FontSize', 7); end
        end
    end
    % y is shared across ALL tiles of both ISI figures (plot_twoap_cs_average
    % -> share_ylim), so every panel is on the largest trace's scale.
    setappdata(fig, 'trace_axes', reshape(rowAxes(isgraphics(rowAxes, 'axes')), 1, []));
    title(tl, {['2AP Coverslip grand-average dF/F \bullet ' bandLabel ' \bullet colour = mean of all ' ...
        'QC-passed ROIs (no response gate); grey = individual region (c#) means; y shared across all panels'], ...
        sprintf(['PPR = primary (gate-either \\SigmaA2/\\SigmaA1, A2 decay-subtracted; PPR fig panel 1); ' ...
        'trace = A2/A1 of this avg trace, decay-subtracted (panel 3); Q = total charge (dF/F\\cdotms) over ' ...
        '[stim1, stim2+%g ms] (panel 2)'], config.twoAP.auc_post_ms), ...
        sprintf(['blue / red dashed = stim1 / stim2; grey band = Q window (stim1 to stim2 + %g ms); ' ...
        'blue / red bands = peak-search windows; x in %s'], config.twoAP.auc_post_ms, tax.unit)}, ...
        'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'tex');
end

function [ppr_pk, q_pair, A1, A2, A2_abs, amp] = trace_ppr(trace, stim1, stim2, config, pre, template)
    % Paired-pulse metrics of a single (already-averaged) trace. Amplitudes come
    % from twoap_pulse_amplitudes, the SAME function the per-ROI analyzer uses
    % (same windows, local-max peak rule and residual). The baseline is the mean
    % over pre (the files' prestim_idx, i.e. the analyzer's baseline frames
    % with the pre-stim guard dropped); 1:stim1-1 if pre is not given.
    % template (optional) enables the decay-subtracted amplitude in amp. Returns:
    %   ppr_pk  - peak PPR = A2 / A1
    %   q_pair  - TOTAL CHARGE (dF/F*ms) over ONE window encompassing both pulses,
    %             [stim1, stim2 + auc_post_ms], above the single pre-stim1 baseline.
    %             Not a ratio; the window (ISI + post) lengthens with ISI, so
    %             baseline drift accumulates into it at long ISI.
    %   A1, A2  - pulse-1/2 peak amplitude (dF/F); A2 above the local pre-stim2 residual
    %   A2_abs  - pulse-2 peak above the ORIGINAL pre-stim1 baseline. A2_abs - A2 =
    %             how far pulse-1's decay has lifted the baseline by the time pulse 2
    %             fires; they converge once the residual settles (long ISI).
    %   amp     - full twoap_pulse_amplitudes struct (all five A1/A2 methods)
    if nargin < 5 || isempty(pre), pre = 1:max(stim1 - 1, 1); end
    if nargin < 6, template = []; end
    fd = config.frame_dur_ms; nF = numel(trace);

    pv   = trace(pre); pv = pv(isfinite(pv));
    base = mean(pv);
    amp  = twoap_pulse_amplitudes(trace, stim1, stim2, base, config, template);
    A1 = amp.peak_max.a1; A2 = amp.peak_max.a2; A2_abs = amp.a2_abs;
    ppr_pk = A2 / A1;

    % Total charge over the encompassing window (above the single pre-stim1 baseline)
    [q_lo, q_hi] = twoap_charge_window(stim1, stim2, nF, config.twoAP.auc_post_ms, fd);
    sq   = trace(q_lo:q_hi) - base; sq = sq(isfinite(sq));
    if numel(sq) >= 2, q_pair = trapz(sq) * fd; else, q_pair = NaN; end
end

function [ppr, q_pair] = trace_ppr_primary(trace, stim1, stim2, config, pre, template)
    % PPR of one average trace on the PRIMARY amplitude basis (A2 decay-subtracted
    % with the group's pulse-1 template), so every figure that annotates an
    % average-trace PPR matches the PPR figure. Falls back to peak-max A2 only
    % when the group has no template (longest ISI < 250 ms).
    [ppr, q_pair, ~, ~, ~, amp] = trace_ppr(trace, stim1, stim2, config, pre, template);
    if ~isempty(template), ppr = amp.decay_sub.a2 / amp.decay_sub.a1; end
end

function S = compute_rom_by_isi(gd, config)
    % Ungated average-trace metrics at each ISI (peak-max PPR is a legacy
    % reference; figures use q_pair, PPR comes from ppr_diag.decay_sub),
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
        [pk, q, A1, A2, A2abs] = trace_ppr(ca(k).trace, stim1, ca(k).stim2, config, ca(k).pre);
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
%% ======================================================================
%% Sync / async release (event binning)
%% ======================================================================
function SA = compute_sync_async(gd, config)
    % Sync / async release per coverslip x ISI over EVERY QC-passed ROI: the event
    % detector (twoap_sync_async_events) is its own 5-sigma gate, so no responder
    % gate is added on top (that would drop the weak, facilitating ROIs). Events
    % are binned by onset latency (<= events.sync_cutoff_ms = sync).
    %   p_sync1/2         fraction of ROIs with >= 1 sync event after pulse 1/2
    %   sync_ppf          p_sync2 / p_sync1 (paired-pulse change in sync release)
    %   async1/2_per_roi  async events per ROI minus the spontaneous rate (pre-stim1
    %                     baseline) x async-bin length
    %   async1/2_hz       the same per second of async bin
    %   async_ppf         async2 / async1 per ROI, only where both async bins have the
    %                     same length (pulse-1 window not cut short by stim2)
    %   async_frac2       async share of the events after pulse 2
    %   sync_prom / async_prom  mean event rise above its trough (dF/F)
    % Pulse-1 async events that land after stim2 are counted in the pulse-2 window.
    ev = gd.events;
    SA = struct('enabled', ev.enabled);
    if ~ev.enabled, return; end
    fd   = config.frame_dur_ms;
    keys = cellfun(@(c, x) sprintf('%s|%g', c, x), gd.coverslip_of_roi, num2cell(gd.isi_of_roi), ...
        'UniformOutput', false);
    ok   = gd.qc & isfinite(ev.n_sync1);
    uk   = unique(keys(ok), 'stable');
    nK   = numel(uk);
    z    = nan(1, nK);
    R = struct('coverslip', {cellfun(@(k) extractBefore(k, '|'), uk, 'UniformOutput', false)}, ...
        'isi_ms', cellfun(@(k) str2double(extractAfter(k, '|')), uk), 'n_roi', zeros(1, nK), ...
        'p_sync1', z, 'p_sync2', z, 'sync_ppf', z, 'async1_ms', z, 'async2_ms', z, 'spont_hz', z, ...
        'async1_per_roi', z, 'async2_per_roi', z, 'async1_hz', z, 'async2_hz', z, ...
        'async_ppf', z, 'async_frac2', z, 'sync_prom', z, 'async_prom', z);
    L = ev.list;
    for k = 1:nK
        sel = ok & strcmp(keys, uk{k});
        n = sum(sel); R.n_roi(k) = n;
        R.p_sync1(k) = mean(ev.n_sync1(sel) > 0);
        R.p_sync2(k) = mean(ev.n_sync2(sel) > 0);
        if R.p_sync1(k) > 0, R.sync_ppf(k) = R.p_sync2(k) / R.p_sync1(k); end
        sp = sum(ev.n_spont(sel)) / sum(ev.spont_ms(sel));     % spont events per ms per ROI
        R.spont_hz(k) = sp * 1000;
        w1 = median(ev.async1_ms(sel)); w2 = median(ev.async2_ms(sel));
        R.async1_ms(k) = w1; R.async2_ms(k) = w2;
        if w1 > 0
            R.async1_per_roi(k) = mean(ev.n_async1(sel)) - sp * w1;
            R.async1_hz(k)      = R.async1_per_roi(k) / w1 * 1000;
        end
        if w2 > 0
            R.async2_per_roi(k) = mean(ev.n_async2(sel)) - sp * w2;
            R.async2_hz(k)      = R.async2_per_roi(k) / w2 * 1000;
        end
        if w1 > 0 && abs(w1 - w2) < fd / 2 && R.async1_per_roi(k) > 0
            R.async_ppf(k) = R.async2_per_roi(k) / R.async1_per_roi(k);
        end
        s2 = sum(ev.n_sync2(sel)); a2 = R.async2_per_roi(k) * n;
        if isfinite(a2) && s2 + a2 > 0, R.async_frac2(k) = a2 / (s2 + a2); end
        inL = ismember(L.roi, find(sel));
        R.sync_prom(k)  = mean(L.prom(inL & L.is_sync));
        R.async_prom(k) = mean(L.prom(inL & ~L.is_sync));
    end
    SA.cs_rows = R;
    for f = {'p_sync1', 'p_sync2', 'sync_ppf', 'async1_per_roi', 'async2_per_roi', ...
             'async1_hz', 'async2_hz', 'async_ppf', 'async_frac2', 'sync_prom', 'async_prom'}
        SA.(f{1}) = summarize_by_isi(R.isi_ms, R.(f{1}));
    end

    % Onset-latency histograms (events per QC ROI per frame), pooled over the ISIs whose
    % pulse-1 window is the full async_window_ms, so both pulses have the same window.
    full  = ok & abs(ev.async1_ms - ev.async2_ms) < fd / 2;
    aw_f  = round(config.twoAP.events.async_window_ms / fd);
    edges = ((0:aw_f + 1) - 0.5) * fd;
    inF   = ismember(L.roi, find(full));
    nFull = max(sum(full), 1);
    SA.hist = struct('centers', (0:aw_f) * fd, 'n_roi', sum(full), ...
        'isi', unique(gd.isi_of_roi(full)), ...
        'p1', histcounts(L.onset_ms(inF & L.pulse == 1), edges) / nFull, ...
        'p2', histcounts(L.onset_ms(inF & L.pulse == 2), edges) / nFull, ...
        'spont', sum(ev.n_spont(full)) / max(sum(ev.spont_ms(full)), eps) * fd);
end

function fig = plot_twoap_sync_async(groups, config)
    % Sync / async event binning: latency histograms (top) and per-coverslip
    % sync / async metrics vs ISI (bottom). Groups without event binning (slow
    % frames) are skipped; returns [] when none has it.
    fig = [];
    on = arrayfun(@(g) isfield(g.sync_async, 'enabled') && g.sync_async.enabled, groups);
    if ~any(on), return; end
    groups = groups(on);
    e  = config.twoAP.events;
    fd = config.frame_dur_ms;
    fig = figure('Name', '2AP sync vs async', 'Position', [80 60 1600 820], 'Color', 'w');

    % --- Top: event onset-latency histograms, pulse 1 / pulse 2 ---
    pulses = {'p1', 'p2'};
    for pp = 1:2
        ax = subplot(2, 3, pp); hold on;
        % Side-by-side bars per frame bin (all groups share the frame grid)
        H  = groups(1).sync_async.hist;
        nb = min(arrayfun(@(g) numel(g.sync_async.hist.centers), groups));
        Y  = cell2mat(arrayfun(@(g) reshape(g.sync_async.hist.(pulses{pp})(1:nb), [], 1), ...
            groups, 'UniformOutput', false));
        sp  = max(arrayfun(@(g) g.sync_async.hist.spont, groups));
        pos = [Y(Y > 0); sp(sp > 0)];
        hb = bar(ax, H.centers(1:nb), Y, 1, 'grouped', 'EdgeColor', 'none');
        for g = 1:numel(groups)
            hb(g).FaceColor   = groups(g).color;
            hb(g).DisplayName = sprintf('%s (%d ROIs)', groups(g).name, groups(g).sync_async.hist.n_roi);
        end
        if ~isempty(pos)   % log axis: bars rise from below the smallest value / spont line
            hb(1).BaseValue = min(pos) / 2;
            hb(1).BaseLine.Visible = 'off';
        end
        if sp > 0
            yline(ax, sp, 'k:', 'spont', 'HandleVisibility', 'off', 'LabelHorizontalAlignment', 'right');
        end
        xline(ax, e.sync_cutoff_ms + fd / 2, 'k--', 'sync | async', 'HandleVisibility', 'off', ...
            'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'center', ...
            'LabelOrientation', 'horizontal');
        set(ax, 'YScale', 'log'); grid(ax, 'on'); box(ax, 'on');
        xlim(ax, [-fd, e.async_window_ms + fd]);
        xlabel(ax, sprintf('Event onset latency after stim %d (ms)', pp));
        ylabel(ax, 'Events per ROI per frame');
        legend(ax, 'Location', 'northeast');
        title(ax, {sprintf('\\bfPulse %d event onset', pp), ...
            sprintf('ISI \\geq %g ms pooled \\bullet all QC ROIs', min(H.isi))}, ...
            'Interpreter', 'tex');
    end

    % --- Top right: async share of pulse-2 events ---
    ax = subplot(2, 3, 3); hold on;
    for g = 1:numel(groups)
        S = groups(g).sync_async.async_frac2;
        draw_ppr_series(ax, S.isi, S, groups(g).color, groups(g).name, 'o');
    end
    finish_metric_axis(ax, 'Async fraction of pulse-2 events');
    title(ax, {'\bfAsync share after pulse 2', 'spont-subtracted'}, 'Interpreter', 'tex');

    % --- Bottom: P(sync), sync PPF, async rate ---
    ax = subplot(2, 3, 4); hold on;
    for g = 1:numel(groups)
        SA = groups(g).sync_async;
        draw_ppr_series(ax, SA.p_sync1.isi, SA.p_sync1, groups(g).color, [groups(g).name ' P1'], 'o', '--');
        draw_ppr_series(ax, SA.p_sync2.isi, SA.p_sync2, groups(g).color, [groups(g).name ' P2'], 's', '-');
    end
    finish_metric_axis(ax, 'P(sync event) per ROI');
    title(ax, {'\bfSync release probability', ...
        sprintf('fraction of ROIs with an event onset \\leq %g ms \\bullet dashed P1, solid P2', ...
        e.sync_cutoff_ms)}, ...
        'Interpreter', 'tex');

    ax = subplot(2, 3, 5); hold on;
    for g = 1:numel(groups)
        S = groups(g).sync_async.sync_ppf;
        draw_ppr_series(ax, S.isi, S, groups(g).color, groups(g).name, 'o');
    end
    finish_ppr_axis(ax, 'P(sync) P2 / P1');
    ylim(ax, [0 max(ax.YLim(2), 1.2)]);
    title(ax, {'\bfSync PPF', 'P(sync event) after pulse 2 / pulse 1'}, 'Interpreter', 'tex');

    ax = subplot(2, 3, 6); hold on;
    for g = 1:numel(groups)
        SA = groups(g).sync_async;
        draw_ppr_series(ax, SA.async1_hz.isi, SA.async1_hz, groups(g).color, [groups(g).name ' P1'], 'o', '--');
        draw_ppr_series(ax, SA.async2_hz.isi, SA.async2_hz, groups(g).color, [groups(g).name ' P2'], 's', '-');
    end
    finish_metric_axis(ax, 'Async events / s per ROI');
    title(ax, {'\bfAsync release rate', ...
        sprintf('onsets %.0f-%g ms after each stim (P1 stops at stim2), minus spont', ...
        async_start_ms(config), e.async_window_ms)}, 'Interpreter', 'tex');

    sgtitle({'2AP sync vs async release (binned on event onset latency)', ...
        sprintf(['events: local max \\geq %g\\sigma rising \\geq %g\\sigma above the trough before it' ...
        ' \\bullet per-coverslip mean \\pm SEM (dots = coverslips)'], ...
        config.twoAP.upper_threshold_sigma, e.min_prominence_sigma)}, ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'tex');
end

function ex = pick_sync_async_examples(frs, config)
    % Example ROIs for the sync / async trace figure: a seeded random sample of
    % ROIs with >= 1 async event and of sync-only ROIs (counts from
    % visualization.montage.twoAP_sync_async_rois), all ISIs, sorted by ISI. A
    % random sample (not the largest events) shows what a typical counted event
    % looks like. Keeps each ROI's trace and the analyzer's own events, so the
    % markers are exactly the events that were counted.
    nPick = config.visualization.montage.twoAP_sync_async_rois;
    ex = struct('trace', {}, 'base', {}, 'upper', {}, 'stim1', {}, 'stim2', {}, 'isi', {}, ...
        'label', {}, 'roi', {}, 'kind', {}, 'pulse', {}, 'lat_ms', {}, 'onset_ms', {}, 'is_sync', {});
    cand = zeros(0, 3);                              % [file, roi, has async]
    for i = 1:numel(frs)
        fr = frs{i};
        if isempty(fr) || ~fr.twoap.events_enabled, continue; end
        L = fr.twoap.events;
        a = unique(L.roi(~L.is_sync));
        s = setdiff(unique(L.roi(L.is_sync)), a);
        cand = [cand; repmat(i, numel(a), 1), a(:), ones(numel(a), 1); ...
                      repmat(i, numel(s), 1), s(:), zeros(numel(s), 1)]; %#ok<AGROW>
    end
    rs = RandStream('mt19937ar', 'Seed', 1);
    kinds = {'async', 'sync only'};
    for kk = 1:2
        c = find(cand(:, 3) == (kk == 1));
        c = c(randperm(rs, numel(c), min(nPick(kk), numel(c))));
        isi = arrayfun(@(j) frs{cand(j, 1)}.isi_ms, c);
        [~, o] = sort(isi);
        for j = c(o)'
            fr = frs{cand(j, 1)}; t = fr.twoap; roi = cand(j, 2);
            L  = t.events; m = L.roi == roi;
            ex(end+1) = struct('trace', fr.dfof_data(:, roi), 'base', t.prestim_base(roi), ...  %#ok<AGROW>
                'upper', t.upper_threshold(roi), 'stim1', t.stim1_frame, 'stim2', t.stim2_frame, ...
                'isi', fr.isi_ms, 'label', sprintf('%s-%s', parse_coverslip_id(fr.filename), ...
                parse_region_id(fr.filename)), 'roi', roi, 'kind', kinds{kk}, ...
                'pulse', L.pulse(m), 'lat_ms', L.lat_ms(m), 'onset_ms', L.onset_ms(m), ...
                'is_sync', L.is_sync(m));
        end
    end
end

function figs = plot_twoap_sync_async_examples(groups, config)
    % Example ROI traces with every counted event marked at its peak: filled
    % circle = sync (onset <= sync_cutoff_ms), open triangle = async; blue = after
    % stim1, red = after stim2; open square + faint solid line = the event's onset
    % frame (the one binned on); stim lines faint dashed. Light blue / red bands
    % = sync onset bins, grey = async onset bins. When the
    % ISI is long, the stretch between the two pulse windows is cut out (dotted
    % break line) so both windows are drawn at the same ms scale.
    figs = [];
    e    = config.twoAP.events;
    fd   = config.frame_dur_ms;
    aw   = e.async_window_ms;
    pre  = 30; post = 30;                            % ms shown before each stim / after each window
    pc   = {[0.2 0.4 1.0], [1.0 0.3 0.3]};
    grey = [0.5 0.5 0.5];
    tS   = round(e.sync_cutoff_ms / fd) * fd + fd / 2;   % sync | async onset boundary (ms after a stim)
    for g = 1:numel(groups)
        SA = groups(g).sync_async;
        if ~(isfield(SA, 'enabled') && SA.enabled) || ~isfield(SA, 'examples') || isempty(SA.examples)
            continue;
        end
        ex = SA.examples;
        nc = 4; nr = ceil(numel(ex) / nc);
        fig = figure('Name', sprintf('2AP sync/async examples %s', groups(g).name), ...
            'Position', [60 60 1800 260 * nr], 'Color', 'w');
        figs = [figs, fig]; %#ok<AGROW>
        tl = tiledlayout(fig, nr, nc, 'TileSpacing', 'compact', 'Padding', 'compact');
        for k = 1:numel(ex)
            x  = ex(k);
            ax = nexttile(tl); hold(ax, 'on');
            isi_ms = (x.stim2 - x.stim1) * fd;
            w1_end = min(isi_ms - fd, aw) + fd / 2;      % pulse-1 window end (stops before stim2)
            t  = ((1:numel(x.trace))' - x.stim1) * fd;   % ms after stim1
            gap = max(isi_ms - (aw + post + pre), 0);    % ms cut out between the windows
            xs = t;
            if gap > 0
                xs(t > aw + post & t < isi_ms - pre) = NaN;
                xs(t >= isi_ms - pre) = t(t >= isi_ms - pre) - gap;
            end
            s2x = isi_ms - gap;                          % stim2 position on the cut axis
            xregion(ax, 0, tS, 'FaceColor', pc{1}, 'FaceAlpha', 0.10);
            if w1_end > tS, xregion(ax, tS, w1_end, 'FaceColor', grey, 'FaceAlpha', 0.10); end
            xregion(ax, s2x, s2x + tS, 'FaceColor', pc{2}, 'FaceAlpha', 0.10);
            xregion(ax, s2x + tS, s2x + aw + fd / 2, 'FaceColor', grey, 'FaceAlpha', 0.10);
            plot(ax, xs, x.trace, 'k-', 'LineWidth', 1);
            yline(ax, x.base, '-', 'Color', [0.6 0.6 0.6]);
            yline(ax, x.upper, 'g:', 'LineWidth', 1);
            xline(ax, 0, '--', 'Color', pc{1}, 'LineWidth', 1, 'Alpha', 0.35);
            xline(ax, s2x, '--', 'Color', pc{2}, 'LineWidth', 1, 'Alpha', 0.35);
            if gap > 0, xline(ax, aw + post, 'k:', 'LineWidth', 1.2); end
            % Events: frame = stim + latency
            for j = 1:numel(x.lat_ms)
                stim = x.stim1; if x.pulse(j) == 2, stim = x.stim2; end
                f = stim + round(x.lat_ms(j) / fd);
                fo = stim + round(x.onset_ms(j) / fd);
                % Onset frame (the one binned on): faint vertical line + open square
                xline(ax, xs(fo), '-', 'Color', pc{x.pulse(j)}, 'LineWidth', 1, 'Alpha', 0.3);
                plot(ax, xs(fo), x.trace(fo), 's', 'Color', pc{x.pulse(j)}, 'MarkerSize', 6, 'LineWidth', 1.2);
                if x.is_sync(j)
                    plot(ax, xs(f), x.trace(f), 'o', 'Color', pc{x.pulse(j)}, ...
                        'MarkerFaceColor', pc{x.pulse(j)}, 'MarkerSize', 6);
                else
                    plot(ax, xs(f), x.trace(f), '^', 'Color', pc{x.pulse(j)}, ...
                        'MarkerSize', 7, 'LineWidth', 1.4);
                end
            end
            xl = [-pre, s2x + aw + post];
            vis = xs >= xl(1) & xs <= xl(2) & isfinite(x.trace);
            v = x.trace(vis);
            if ~isempty(v)
                pad = 0.08 * max(max(v) - min(v), eps);
                ylim(ax, [min(v) - pad, max(v) + pad]);
            end
            xlim(ax, xl);
            if gap > 0
                tk = 0:50:aw;
                set(ax, 'XTick', [tk, s2x + tk], 'XTickLabel', string([tk, tk]));
            end
            n1 = [sum(x.pulse == 1 & x.is_sync), sum(x.pulse == 1 & ~x.is_sync)];
            n2 = [sum(x.pulse == 2 & x.is_sync), sum(x.pulse == 2 & ~x.is_sync)];
            title(ax, sprintf('[%s] %s ROI%d  ISI %g ms  P1 %d/%d  P2 %d/%d', x.kind, x.label, ...
                x.roi, x.isi, n1, n2), 'FontSize', 8, 'FontWeight', 'normal');
            set(ax, 'FontSize', 7); grid(ax, 'on'); box(ax, 'on');
        end
        xlabel(tl, 'ms after each stim (blue = stim1, red = stim2; dotted = axis break)');
        ylabel(tl, 'dF/F');
        title(tl, {sprintf(['2AP %s \\bullet random example ROIs with async events, then sync-only ' ...
            'ROIs \\bullet tile titles: events sync/async after P1 and P2'], groups(g).name), ...
            sprintf(['\\bullet = sync peak (onset \\leq %g ms), \\Delta = async peak (onset to %g ms; ' ...
            'P1 stops at stim2); square + faint solid line = onset frame; faint dashed = stims'], ...
            e.sync_cutoff_ms, aw), ...
            sprintf(['blue / red bands = sync onset bins, grey = async onset bins; green dotted = upper ' ...
            '%g\\sigma; grey line = baseline'], config.twoAP.upper_threshold_sigma)}, ...
            'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'tex');
    end
end

function t = async_start_ms(config)
    % First async frame after a stim, in ms (one frame past the sync cutoff).
    fd = config.frame_dur_ms;
    t  = (round(config.twoAP.events.sync_cutoff_ms / fd) + 1) * fd;
end

function finish_metric_axis(ax, ylab)
    % Log-ISI axis for a non-ratio metric (no PPR = 1 line), y-min 0.
    set(ax, 'XScale', 'log');
    xlabel(ax, 'Inter-stimulus interval (ms)'); ylabel(ax, ylab);
    grid(ax, 'on'); box(ax, 'on'); legend(ax, 'Location', 'best');
    ylim(ax, [0 max(ax.YLim(2), eps)]);
end

function print_sync_async(groups, config)
    on = arrayfun(@(g) isfield(g.sync_async, 'enabled') && g.sync_async.enabled, groups);
    e  = config.twoAP.events;
    if ~any(on)
        fprintf('\nSync/async binning skipped: frame time %.4g ms > events.max_frame_ms %g.\n', ...
            config.frame_dur_ms, e.max_frame_ms);
        return;
    end
    fprintf('\n=== 2AP SYNC / ASYNC (event onset <= %g ms = sync; per-coverslip mean) ===\n', ...
        e.sync_cutoff_ms);
    fprintf('All QC-passed ROIs. Async = events/s per ROI in the async bin minus the pre-stim\n');
    fprintf(' spontaneous rate. The P1 async window stops at stim2 (none at ISI <= %.0f ms).\n', ...
        async_start_ms(config));
    for g = find(on)
        SA = groups(g).sync_async;
        fprintf('%s: spont %.4f ev/s per ROI\n', groups(g).name, mean(SA.cs_rows.spont_hz, 'omitnan'));
        for isi = unique(SA.cs_rows.isi_ms)
            [ps1, n] = mean_at_isi(SA.p_sync1, isi);
            fprintf(['   ISI %5g ms (N=%d cs): P(sync) P1 %.3f P2 %.3f (PPF %.2f) | async/s P1 %.3f ' ...
                'P2 %.3f | async share P2 %.2f\n'], isi, n, ps1, mean_at_isi(SA.p_sync2, isi), ...
                mean_at_isi(SA.sync_ppf, isi), mean_at_isi(SA.async1_hz, isi), ...
                mean_at_isi(SA.async2_hz, isi), mean_at_isi(SA.async_frac2, isi));
        end
    end
end

function T = build_twoap_sync_async_table(groups)
    % One row per group x coverslip x ISI of the sync / async metrics.
    T = table();
    for g = 1:numel(groups)
        SA = groups(g).sync_async;
        if ~(isfield(SA, 'enabled') && SA.enabled), continue; end
        r = SA.cs_rows;
        n = numel(r.isi_ms);
        if n == 0, continue; end
        f  = setdiff(fieldnames(r), {'coverslip'}, 'stable');
        Tg = table(repmat({groups(g).genotype}, n, 1), repmat({groups(g).drug}, n, 1), r.coverslip(:), ...
            'VariableNames', {'genotype', 'drug', 'coverslip'});
        for k = 1:numel(f), Tg.(f{k}) = r.(f{k})(:); end
        T = [T; sortrows(Tg, {'coverslip', 'isi_ms'})]; %#ok<AGROW>
    end
end

function T = build_twoap_table(groups)
    grp = {}; drug = {}; fname = {}; cs = {}; roi = []; isi = []; resp1 = [];
    pmode = {}; ichk = {};
    a1 = []; a2 = []; pmax = []; a1s = []; a2s = []; psync = []; a2p = []; pp1 = []; sd = [];
    q = []; resp2 = []; either = []; a2d = []; pdec = [];
    ns1 = []; na1 = []; ns2 = []; na2 = []; nsp = [];
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
        % Gate-either membership + decay-subtracted A2 (primary PPR inputs).
        % peak1_amp is the A1 the primary uses (decay subtraction leaves A1 as is).
        resp2  = [resp2;  gd.ppr_diag.all.gate_A2(:)];        %#ok<AGROW>
        either = [either; gd.ppr_diag.all.gate_either(:)];    %#ok<AGROW>
        a2d    = [a2d;    gd.ppr_diag.all.a2_decay_sub(:)];   %#ok<AGROW>
        pd     = gd.ppr_diag.all.a2_decay_sub(:) ./ gd.peak1_amp(:);
        pd(~(gd.peak1_amp(:) > 0) | ~gd.ppr_diag.all.gate_either(:)) = NaN;
        pdec   = [pdec;   pd];                                %#ok<AGROW>
        ns1    = [ns1;    gd.events.n_sync1(:)];              %#ok<AGROW>
        na1    = [na1;    gd.events.n_async1(:)];             %#ok<AGROW>
        ns2    = [ns2;    gd.events.n_sync2(:)];              %#ok<AGROW>
        na2    = [na2;    gd.events.n_async2(:)];             %#ok<AGROW>
        nsp    = [nsp;    gd.events.n_spont(:)];              %#ok<AGROW>
    end
    T = table(grp, drug, fname, cs, roi, isi, pmode, logical(resp1), logical(resp2), logical(either), ...
        a1, a2, a2d, pdec, pmax, a1s, a2s, psync, a2p, pp1, q, ns1, na1, ns2, na2, nsp, sd, ichk, ...
        'VariableNames', {'genotype', 'drug', 'filename', 'coverslip', 'roi', 'isi_ms', 'ppf_mode', ...
        'is_responder1', 'is_responder2_incr', 'gate_either', ...
        'peak1_amp', 'peak2_amp', 'a2_decay_sub', 'ppr_decay_sub', 'ppr_max', 'a1_sync', 'a2_sync', 'ppr_sync', ...
        'a2_p1lock', 'ppr_p1lock', 'q_pair', 'n_sync1', 'n_async1', 'n_sync2', 'n_async2', 'n_spont', ...
        'prestim_sd', 'isi_check'});
end

function T = build_twoap_cs_table(groups)
    % One row per group x coverslip x ISI: the PRIMARY PPR (pooled SumA2/SumA1 over
    % the gate-either ROIs, A2 decay-subtracted), its inputs, gate counts and the
    % ungated average-trace reference. This is the statistical-unit table.
    T = table();
    for g = 1:numel(groups)
        r = groups(g).ppr_diag.cs_rows;
        n = numel(r.isi_ms);
        if n == 0, continue; end
        Tg = table(repmat({groups(g).genotype}, n, 1), repmat({groups(g).drug}, n, 1), ...
            r.coverslip(:), r.isi_ms(:), r.n_either(:), r.n_A1only(:), r.n_A2only(:), ...
            r.sum_a1(:), r.sum_a2_decay_sub(:), r.mean_a1(:), r.mean_a2_decay_sub(:), r.mean_a2_abs(:), r.ppr_primary(:), ...
            r.n_excluded(:), r.ppr_sub_threshold(:), r.ppr_ref_no_gate(:), ...
            'VariableNames', {'genotype', 'drug', 'coverslip', 'isi_ms', 'n_gate_either', ...
            'n_A1_only', 'n_A2_only', 'sum_a1', 'sum_a2_decay_sub', 'mean_a1', 'mean_a2_decay_sub', 'mean_a2_abs', ...
            'ppr_primary', 'n_sub_threshold', 'ppr_sub_threshold', 'ppr_ref_no_gate'});
        T = [T; sortrows(Tg, {'coverslip', 'isi_ms'})]; %#ok<AGROW>
    end
end

function print_summary_2ap(groups, config)
    fprintf('\n=== 2AP SUMMARY (PPR = A2/A1; per-coverslip mean at each ISI) ===\n');
    fprintf('(Statistical unit = coverslip. Regions/files within a coverslip are\n');
    fprintf(' pseudoreplicates. With 1 animal/genotype these are technical replicates.)\n');
    fprintf('PRIMARY = pooled SumA2/SumA1 over gate-either ROIs (P1 responder OR\n');
    fprintf(' incremental A2 >= %g sigma), A2 decay-subtracted. noGate = ungated avg-trace\n', ...
        config.twoAP.upper_threshold_sigma);
    fprintf(' ratio on the same amplitudes (QC reference). subThr = avg trace of ROIs failing\n');
    fprintf(' both gates (decay-sub). A1/A2 = mean per gate-either ROI (A2 decay-sub).\n');
    for g = 1:numel(groups)
        gd   = groups(g);
        rom  = compute_rom_by_isi(gd, config);
        D    = gd.ppr_diag;
        fprintf('%s: %d files, %d coverslips, %d pulse-1 responders (decay template %g ms ISI)\n', ...
            gd.name, gd.num_files, gd.num_coverslips, gd.n_responders1, D.template_isi);
        for k = 1:numel(gd.isi_list)
            isi = gd.isi_list(k);
            [mP, nP]   = mean_at_isi(D.decay_sub.pooled, isi);
            mRef       = mean_at_isi(D.decay_sub.trace_all, isi);
            mQ         = mean_at_rom(rom.q_pair, rom.isi, isi);
            mA1        = mean_at_isi(D.amp_a1, isi);
            mA2        = mean_at_isi(D.amp_a2, isi);
            mSub       = mean_at_isi(D.sub_threshold, isi);
            ki = find(D.isi == isi, 1); nE = 0;
            if ~isempty(ki), nE = D.n_either(ki); end
            fprintf(['   ISI %5g ms: PPR %.2f (N=%d cs, %d ROIs) | subThr %.2f | noGate %.2f | Q %.1f | ' ...
                'A1=%.4f A2=%.4f\n'], isi, mP, nP, nE, mSub, mRef, mQ, mA1, mA2);
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
