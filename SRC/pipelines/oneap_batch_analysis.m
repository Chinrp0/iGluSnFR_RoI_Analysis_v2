function results = oneap_batch_analysis(dataset, options)
    % ONEAP_BATCH_ANALYSIS - Batch analysis of 1AP evoked responses
    %
    % Takes a dataset loaded by load_dataset_folder (data + Recordings.csv), so
    % all timing - frame interval, frame count, stim frames - comes from the
    % metadata. Every file is processed by oneap_pipeline, which analyses each
    % trial separately; the resulting (ROI, trial) "units" are pooled per group
    % (genotype x drug). Unit-level metrics (amplitude, kinetics, ...) are
    % therefore over SUCCESSFUL trials, and per-ROI release probability across
    % trials is reported alongside. Statistics compare genotypes within each
    % drug condition, on per-coverslip aggregates.
    %
    % Inputs:
    %   dataset - struct from load_dataset_folder
    %   options - struct: verbose, createPlots
    %
    % Output:
    %   results - struct with .groups, .comparisons, .plot_handles,
    %             .oneap_metrics_table (one row per ROI x trial),
    %             .oneap_roi_table (one row per ROI),
    %             .oneap_coverslip_table (one row per coverslip), .processing_info

    if nargin < 2, options = struct(); end
    if ~isfield(options, 'verbose'),     options.verbose = true; end
    if ~isfield(options, 'createPlots'), options.createPlots = true; end
    if ~isstruct(dataset) || ~isfield(dataset, 'files')
        error('oneap_batch_analysis:input', ...
            'oneap_batch_analysis takes a dataset from load_dataset_folder.');
    end

    config = local_config();

    fprintf('=== 1AP Evoked Response Batch Analysis ===\n');
    fprintf('Folder: %s\n', dataset.folder);

    %% Step 1: files + one shared time base
    files = dataset.files;
    if isempty(files)
        error('oneap_batch_analysis:noFiles', 'The dataset has no loadable files.');
    end
    % Groups are pooled onto one time axis, so all files must share the frame
    % interval. Differing intervals are a mixed folder, not something to average.
    fi = arrayfun(@(f) f.meta.frame_interval_ms, files);
    if (max(fi) - min(fi)) > config.metadata.interval_tol_frac * median(fi)
        error('oneap_batch_analysis:mixedInterval', ...
            'Files have different Frame_Interval_ms (%s); analyse them as separate folders.', ...
            strjoin(unique(cellstr(num2str(fi(:), '%.4g'))), ', '));
    end
    pre_f  = round(config.oneAP.epoch_pre_ms  / median(fi));
    post_f = round(config.oneAP.epoch_post_ms / median(fi));
    config = apply_recording_timing(config, median(fi), pre_f + post_f + 1);
    config.oneAP.stim_frame = pre_f + 1;            % stim index within each epoch
    fprintf('  Frame interval %.4g ms (metadata) | epoch %d frames, stim at %d\n', ...
        config.frame_dur_ms, config.expected_frames, config.oneAP.stim_frame);

    %% Step 2: process every file
    fprintf('\nProcessing files...\n');
    nFiles = numel(files);
    fileResults = cell(nFiles, 1);
    tic;
    for i = 1:nFiles
        try
            fileResults{i} = oneap_pipeline(files(i), struct('verbose', false), config);
            fr = fileResults{i};
            qc = fr.oneap.qc_counts;
            rs = fr.roi_summary;
            fprintf('  [%d/%d] %s (%s, %s): %d ROIs x %d trials, %d/%d successful trials (%.0f%%), %d ROIs responded at least once\n', ...
                i, nFiles, fr.filename, fr.genotype, fr.drug, fr.numROIs, fr.numTrials, ...
                fr.oneap.n_responders, fr.numUnits, 100 * fr.oneap.responder_fraction, ...
                sum(rs.n_success >= 1));
            fprintf('        gates: %d sync-peak resp, %d async-peak resp, %d subthreshold, %d noisy, %d excluded\n', ...
                qc.responder_sync, qc.responder_async, qc.subthreshold, qc.noise_excluded, qc.excluded);
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

    %% Step 3: aggregate per group (genotype x drug)
    fprintf('\nAggregating by genotype x drug...\n');
    ok = ~cellfun(@isempty, fileResults(:)');
    genotypes = {files.genotype};
    drugs     = {files.drug};
    [geno_list, drug_list] = genotype_drug_order(genotypes(ok), drugs(ok), config);
    groups = [];
    for d = 1:numel(drug_list)
        for g = 1:numel(geno_list)
            mask = ok & strcmp(genotypes, geno_list{g}) & strcmp(drugs, drug_list{d});
            if ~any(mask), continue; end
            name = sprintf('%s %s', geno_list{g}, drug_list{d});
            gd = aggregate_group(fileResults(mask), name, ...
                group_color(geno_list{g}, d, numel(drug_list), config));
            gd.genotype = geno_list{g};
            gd.drug     = drug_list{d};
            gd.avg_fit  = fit_group_double(gd, config);   % averaged double-exp decay
            gd.cs_avg   = coverslip_average_traces(gd, config);
            if isempty(groups), groups = gd; else, groups(end+1) = gd; end %#ok<AGROW>
            fprintf('  %s: %d files, %d coverslips, %d ROIs, %d/%d successful trials (%.0f%%)\n', ...
                gd.name, gd.num_files, gd.num_coverslips, gd.n_rois, gd.n_responders, ...
                gd.total_rois, 100 * gd.n_responders / max(gd.total_rois, 1));
        end
    end

    %% Step 4: statistics - genotype pairs within each drug condition
    comparisons = struct('group_a', {}, 'group_b', {}, 'drug', {}, 'stats', {});
    for d = 1:numel(drug_list)
        idx = find(strcmp({groups.drug}, drug_list{d}));
        for a = 1:numel(idx)
            for b = a+1:numel(idx)
                ga = groups(idx(a)); gb = groups(idx(b));
                fprintf('\n%s vs %s (Mann-Whitney U, per-coverslip means)...\n', ga.name, gb.name);
                comparisons(end+1) = struct('group_a', ga.name, 'group_b', gb.name, ...
                    'drug', drug_list{d}, 'stats', compare_groups(ga, gb)); %#ok<AGROW>
            end
        end
    end
    % The figures draw one p-value bar between groups 1 and 2, so they are only
    % annotated when exactly those two groups exist; otherwise see comparisons.
    if numel(groups) == 2 && numel(comparisons) == 1
        comparison = comparisons(1).stats;
    else
        comparison = struct();
    end

    %% Step 5: plots
    plot_handles = struct();
    if options.createPlots
        fprintf('\nCreating 1AP figures...\n');
        plot_handles.average_traces = plot_oneap_average_traces(groups, config);
        plot_handles.cs_average     = plot_oneap_cs_average(groups, config);
        plot_handles.amplitude      = plot_oneap_amplitude(groups, comparison);
        plot_handles.kinetics       = plot_oneap_kinetics(groups, comparison);
        plot_handles.auc_split      = plot_oneap_auc_split(groups, comparison);
        plot_handles.sync_async     = plot_oneap_sync_async(groups, config);
        plot_handles.responders     = plot_oneap_responders(groups, comparison);
        plot_handles.release_prob   = plot_oneap_release_prob(groups, comparison);
        plot_handles.async_fraction = plot_oneap_async_fraction(groups, comparison, config);
        plot_handles.qc             = plot_oneap_qc(groups, config);
        plot_handles.montage        = plot_oneap_montage(groups, config);
        plot_handles.marginal       = plot_oneap_marginal(groups, config);
        fprintf('  Figures created.\n');
    end

    %% Step 6: metrics tables
    oneap_metrics_table = build_metrics_table(groups);
    oneap_roi_table     = build_roi_table(groups);
    oneap_coverslip_table = build_coverslip_table(groups);

    %% Assemble results
    results = struct();
    results.dataType            = '1AP';
    results.groups              = groups;
    results.comparisons         = comparisons;
    results.comparison          = comparison;
    results.plot_handles        = plot_handles;
    results.oneap_metrics_table = oneap_metrics_table;
    results.oneap_roi_table     = oneap_roi_table;
    results.oneap_coverslip_table = oneap_coverslip_table;
    results.processing_info     = struct('folder_path', dataset.folder, ...
        'metadata_file', dataset.metadata_file, ...
        'process_time', process_time, 'timestamp', datetime('now'), ...
        'frame_interval_ms', config.frame_dur_ms, 'epoch_frames', config.expected_frames, ...
        'epoch_stim_frame', config.oneAP.stim_frame);

    print_summary(groups, comparisons);
end

%% ======================================================================
%% Aggregation / stats helpers
%% ======================================================================
function gd = aggregate_group(fileResults, name, color)
    % Columns of every per-"ROI" array here are (ROI, trial) units from
    % oneap_pipeline; unit_roi / unit_trial map them back. The roi_* arrays are
    % one entry per physical ROI (release statistics across its trials).
    gd = struct();
    gd.name = name;
    gd.color = color;
    gd.num_files = numel(fileResults);
    unit_roi = []; unit_trial = []; unit_flags = {};
    roi_p_response = []; roi_potency = []; roi_efficacy = [];
    roi_n_success = []; roi_n_valid = []; roi_coverslip = {}; roi_file = []; roi_index = [];

    dfof_all = [];
    peak_amplitude = []; time_to_peak_ms = []; auc = []; tau_ms = []; tau_r2 = [];
    auc_sync = []; auc_async = [];
    is_responder = []; prestim_sd = []; sync_count = []; async_count = []; n_spont = [];
    release_events = struct('lat_ms', {}, 'onset_ms', {}, 'is_sync', {});
    events_ms = struct('sync', NaN, 'async', NaN, 'spont', NaN); events_enabled = false;
    total_events = []; peak_frame = []; upper_threshold = []; lower_threshold = [];
    qc_category = []; crossed_upper = [];
    tau_A = []; tau_offset = []; tau_peak_frame = [];
    event_peak_times_ms = {};
    file_of_roi = []; filenames = {}; coverslip_of_roi = {}; region_of_roi = {};
    temporal_hist_sum = []; bin_centers_ms = [];

    ref_frames = [];

    for i = 1:numel(fileResults)
        fr = fileResults{i};
        if isempty(fr), continue; end
        o = fr.oneap;
        filenames{end+1} = fr.filename; %#ok<AGROW>
        fi = numel(filenames);   % index into filenames (skipped files leave no gap)

        % Experimental unit = coverslip (Cs#). Multiple regions (c#) on one
        % coverslip are pseudoreplicates and must collapse to it before stats.
        csid = parse_coverslip_id(fr.filename);
        coverslip_of_roi = [coverslip_of_roi, repmat({csid}, 1, fr.numUnits)]; %#ok<AGROW>
        region_of_roi    = [region_of_roi, repmat({parse_region_id(fr.filename)}, 1, fr.numUnits)]; %#ok<AGROW>
        unit_roi   = [unit_roi,   o.unit_roi];    %#ok<AGROW>
        unit_trial = [unit_trial, o.unit_trial];  %#ok<AGROW>
        unit_flags = [unit_flags, o.unit_flags];  %#ok<AGROW>

        rs = fr.roi_summary;
        roi_p_response = [roi_p_response, rs.p_response]; %#ok<AGROW>
        roi_potency    = [roi_potency,    rs.potency];    %#ok<AGROW>
        roi_efficacy   = [roi_efficacy,   rs.efficacy];   %#ok<AGROW>
        roi_n_success  = [roi_n_success,  rs.n_success];  %#ok<AGROW>
        roi_n_valid    = [roi_n_valid,    rs.n_valid];    %#ok<AGROW>
        roi_coverslip  = [roi_coverslip,  repmat({csid}, 1, fr.numROIs)]; %#ok<AGROW>
        roi_file       = [roi_file,       repmat(fi, 1, fr.numROIs)];     %#ok<AGROW>
        roi_index      = [roi_index,      1:fr.numROIs];                  %#ok<AGROW>

        peak_amplitude  = [peak_amplitude,  o.peak_amplitude];   %#ok<AGROW>
        time_to_peak_ms = [time_to_peak_ms, o.time_to_peak_ms];  %#ok<AGROW>
        auc             = [auc,             o.auc];              %#ok<AGROW>
        auc_sync        = [auc_sync,        o.auc_sync];         %#ok<AGROW>
        auc_async       = [auc_async,       o.auc_async];        %#ok<AGROW>
        tau_ms          = [tau_ms,          o.tau_ms];           %#ok<AGROW>
        tau_r2          = [tau_r2,          o.tau_r2];           %#ok<AGROW>
        tau_A           = [tau_A,           o.tau_A];            %#ok<AGROW>
        tau_offset      = [tau_offset,      o.tau_offset];       %#ok<AGROW>
        tau_peak_frame  = [tau_peak_frame,  o.tau_peak_frame];   %#ok<AGROW>
        is_responder    = [is_responder,    o.is_responder];     %#ok<AGROW>
        prestim_sd      = [prestim_sd,      o.prestim_sd];       %#ok<AGROW>
        sync_count      = [sync_count,      o.sync_count];       %#ok<AGROW>
        async_count     = [async_count,     o.async_count];      %#ok<AGROW>
        n_spont         = [n_spont,         o.n_spont];          %#ok<AGROW>
        release_events  = [release_events,  o.release_events];   %#ok<AGROW>
        if o.events_enabled && isfield(o, 'events_ms')
            events_enabled = true; events_ms = o.events_ms;
        end
        total_events    = [total_events,    o.total_events];     %#ok<AGROW>
        peak_frame      = [peak_frame,      o.peak_frame];       %#ok<AGROW>
        upper_threshold = [upper_threshold, o.upper_threshold];  %#ok<AGROW>
        lower_threshold = [lower_threshold, o.lower_threshold];  %#ok<AGROW>
        qc_category     = [qc_category,     o.qc_category];       %#ok<AGROW>
        crossed_upper   = [crossed_upper,   o.crossed_upper];     %#ok<AGROW>
        event_peak_times_ms = [event_peak_times_ms, o.event_peak_times_ms]; %#ok<AGROW>
        file_of_roi     = [file_of_roi,     repmat(fi, 1, fr.numUnits)]; %#ok<AGROW>

        if isempty(temporal_hist_sum)
            temporal_hist_sum = sum(o.temporal_hist, 2);
            bin_centers_ms = o.bin_centers_ms;
        else
            temporal_hist_sum = temporal_hist_sum + sum(o.temporal_hist, 2);
        end

        % dF/F matrix (only stack frame-aligned files)
        if isempty(ref_frames)
            ref_frames = fr.numFrames;
        end
        if fr.numFrames == ref_frames
            dfof_all = [dfof_all, fr.dfof_data]; %#ok<AGROW>
        else
            dfof_all = [dfof_all, nan(ref_frames, fr.numUnits)]; %#ok<AGROW>
        end
    end

    gd.dfof_all          = dfof_all;
    gd.peak_amplitude    = peak_amplitude;
    gd.time_to_peak_ms   = time_to_peak_ms;
    gd.auc               = auc;
    gd.auc_sync          = auc_sync;
    gd.auc_async         = auc_async;
    gd.tau_ms            = tau_ms;
    gd.tau_r2            = tau_r2;
    gd.tau_A             = tau_A;
    gd.tau_offset        = tau_offset;
    gd.tau_peak_frame    = tau_peak_frame;
    gd.is_responder      = logical(is_responder);
    gd.prestim_sd        = prestim_sd;
    gd.sync_count        = sync_count;
    gd.async_count       = async_count;
    gd.n_spont           = n_spont;
    gd.release_events    = release_events;
    gd.events_enabled    = events_enabled;   % sync/async binning ran (fast frames)
    gd.events_ms         = events_ms;        % sync / async / spont bin lengths (ms)
    gd.total_events      = total_events;
    gd.peak_frame        = peak_frame;
    gd.upper_threshold   = upper_threshold;
    gd.lower_threshold   = lower_threshold;
    gd.qc_category       = qc_category;
    gd.crossed_upper     = crossed_upper;
    gd.event_peak_times_ms = event_peak_times_ms;
    gd.file_of_roi       = file_of_roi;
    gd.filenames         = filenames;
    gd.file_ids          = file_ids(filenames);
    gd.coverslip_of_roi  = coverslip_of_roi;
    gd.region_of_roi     = region_of_roi;
    gd.temporal_hist_sum = temporal_hist_sum;
    gd.bin_centers_ms    = bin_centers_ms;
    gd.unit_roi          = unit_roi;
    gd.unit_trial        = unit_trial;
    gd.unit_flags        = unit_flags;
    gd.total_rois        = numel(peak_amplitude);   % (ROI, trial) units
    gd.n_responders      = sum(gd.is_responder);    % successful trials
    gd.num_coverslips    = numel(unique(coverslip_of_roi));

    gd.roi_p_response    = roi_p_response;
    gd.roi_potency       = roi_potency;
    gd.roi_efficacy      = roi_efficacy;
    gd.roi_n_success     = roi_n_success;
    gd.roi_n_valid       = roi_n_valid;
    gd.roi_coverslip     = roi_coverslip;
    gd.roi_file          = roi_file;
    gd.roi_index         = roi_index;
    gd.roi_active        = roi_n_success >= 1;       % responded on >= 1 trial
    gd.n_rois            = numel(roi_p_response);
end

function vals = coverslip_means_roi(gd, field)
    % Per-coverslip mean of a per-ROI release statistic over ACTIVE ROIs (>= 1
    % successful trial). ROIs that never responded are not counted as release
    % sites with P = 0: an ROI with no response on any trial may not be a
    % functional synapse at all, and including them would mix "how many ROIs
    % are synapses" into "how reliably does a synapse release".
    cs   = gd.roi_coverslip;
    base = gd.(field);
    keep = isfinite(base) & gd.roi_active;
    uniq = unique(cs, 'stable');
    vals = nan(1, numel(uniq));
    for i = 1:numel(uniq)
        m = strcmp(cs, uniq{i}) & keep;
        if any(m), vals(i) = mean(base(m)); end
    end
    vals = vals(isfinite(vals));
end

function vals = coverslip_means(gd, field, useResponders)
    % Collapse ROIs to ONE value per coverslip (mean of finite, responder-
    % gated values). This is the per-replicate aggregation that removes the
    % within-coverslip pseudoreplication before any statistical test.
    cs   = gd.coverslip_of_roi;
    base = gd.(field);
    keep = isfinite(base);
    if useResponders, keep = keep & gd.is_responder; end
    uniq = unique(cs, 'stable');
    vals = nan(1, numel(uniq));
    for i = 1:numel(uniq)
        m = strcmp(cs, uniq{i}) & keep;
        if any(m), vals(i) = mean(base(m)); end
    end
    vals = vals(isfinite(vals));
end

function fr = coverslip_resp_frac(gd)
    % Per-coverslip responder fraction (responders / total ROIs on that
    % coverslip). One value per coverslip, for a per-replicate comparison.
    cs   = gd.coverslip_of_roi;
    uniq = unique(cs, 'stable');
    fr   = nan(1, numel(uniq));
    for i = 1:numel(uniq)
        m = strcmp(cs, uniq{i});
        if any(m), fr(i) = mean(gd.is_responder(m)); end
    end
    fr = fr(isfinite(fr));
end

function fr = coverslip_async_frac(gd)
    % Per-coverslip ASYNCHRONOUS-release fraction: pooled over the coverslip's
    % ROIs, async events / (sync + async) events. One value per coverslip, for a
    % per-replicate comparison of the slow / asynchronous tail that is visible in
    % the averaged trace but not captured by the per-ROI mono tau (which is
    % truncated at the lower threshold) or the peak-dominated AUC.
    %
    % Event-pooled (not ROI-averaged) so it matches the post-stim PSTH directly
    % and avoids per-ROI 0/0 ratios; async is spontaneous-subtracted as in 2AP
    % (event_totals). Coverslips with no classified events are dropped (NaN).
    cs   = gd.coverslip_of_roi;
    uniq = unique(cs, 'stable');
    fr   = nan(1, numel(uniq));
    for i = 1:numel(uniq)
        [s, a] = event_totals(gd, strcmp(cs, uniq{i}));
        if s + a > 0, fr(i) = a / (s + a); end
    end
    fr = fr(isfinite(fr));
end

function [s, a, a_raw, sp] = event_totals(gd, mask)
    % Sync and async release-event totals over the (ROI, trial) units in MASK.
    % The async count is spontaneous-subtracted like 2AP's async share: the
    % pre-stim event rate of those units x the async bin length is removed
    % (floored at 0). s / a are NaN when event binning did not run.
    s = NaN; a = NaN; a_raw = NaN; sp = NaN;
    if ~gd.events_enabled, return; end
    if nargin < 2, mask = true(size(gd.sync_count)); end
    s     = sum(gd.sync_count(mask), 'omitnan');
    a_raw = sum(gd.async_count(mask), 'omitnan');
    sp    = sum(gd.n_spont(mask), 'omitnan') * gd.events_ms.async / gd.events_ms.spont;
    a     = max(a_raw - sp, 0);
end

function cs_avg = coverslip_average_traces(gd, config)
    % Per-coverslip average dF/F trace over EVERY scorable (ROI, trial) unit
    % (qc 1-3: successes and failures alike), pooled across the coverslip's
    % regions. Not gated on response, so it is selection-free: its peak is the
    % coverslip's mean evoked response per trial (efficacy-like). Alongside it:
    % the successful-trials-only mean (potency-like) and each region's (c#) mean,
    % to show the within-coverslip spread the average collapses.
    cs_avg = struct('cs', {}, 'trace', {}, 'sem', {}, 'n_units', {}, 'n_rois', {}, ...
        'trace_resp', {}, 'sem_resp', {}, 'n_resp', {}, 'regions', {}, 'peak', {}, 'ttp_ms', {}, 'auc', {}, ...
        'peak_resp', {}, 'auc_resp', {});
    if isempty(gd.dfof_all), return; end
    nU    = size(gd.dfof_all, 2);
    valid = gd.qc_category(1:nU) >= 1 & gd.qc_category(1:nU) <= 3;
    cs    = gd.coverslip_of_roi(1:nU);
    rg    = gd.region_of_roi(1:nU);
    uniq  = unique(cs, 'stable');
    for i = 1:numel(uniq)
        m = strcmp(cs, uniq{i}) & valid;
        if ~any(m), continue; end
        M   = gd.dfof_all(:, m);
        n   = sum(isfinite(M), 2);
        mu  = mean(M, 2, 'omitnan')';
        sem = (std(M, 0, 2, 'omitnan') ./ sqrt(max(n, 1)))';
        mr  = m & gd.is_responder(1:nU);
        tr  = nan(size(mu)); ser = nan(size(mu));
        if any(mr)
            R   = gd.dfof_all(:, mr);
            tr  = mean(R, 2, 'omitnan')';
            ser = (std(R, 0, 2, 'omitnan') ./ sqrt(max(sum(isfinite(R), 2), 1)))';
        end
        regions = struct('label', {}, 'trace', {}, 'n', {});
        ur = unique(rg(m), 'stable');
        for r = 1:numel(ur)
            mm = m & strcmp(rg, ur{r});
            regions(end+1) = struct('label', ur{r}, ...                      %#ok<AGROW>
                'trace', mean(gd.dfof_all(:, mm), 2, 'omitnan')', 'n', nnz(mm));
        end
        [pk, ttp, a]   = avg_trace_metrics(mu, config);
        [pkr, ~, ar]   = avg_trace_metrics(tr, config);
        n_rois = size(unique([gd.file_of_roi(m); gd.unit_roi(m)]', 'rows'), 1);
        cs_avg(end+1) = struct('cs', uniq{i}, 'trace', mu, 'sem', sem, ...    %#ok<AGROW>
            'n_units', nnz(m), 'n_rois', n_rois, 'trace_resp', tr, 'sem_resp', ser, 'n_resp', nnz(mr), ...
            'regions', regions, 'peak', pk, 'ttp_ms', ttp, 'auc', a, ...
            'peak_resp', pkr, 'auc_resp', ar);
    end
end

function [pk, ttp_ms, a] = avg_trace_metrics(trace, config)
    % Peak (dF/F above the pre-stim mean), time to peak (ms) and AUC
    % (dF/F*ms) of an already-averaged trace, both over config.oneAP.auc_window_ms.
    % The peak window is the whole AUC window (not the 50 ms tau_peak_window)
    % because some coverslip averages peak late (~100 ms) and would be missed.
    pk = NaN; ttp_ms = NaN; a = NaN;
    stim = config.oneAP.stim_frame;
    fd   = 1000 / config.frame_rate;
    nF   = numel(trace);
    if stim < 2 || stim > nF, return; end
    pre  = trace(1:stim-1); base = mean(pre(isfinite(pre)));
    if ~isfinite(base), return; end
    aw = config.oneAP.auc_window_ms;
    f0 = max(stim + round(aw(1) / fd), 1); f1 = min(stim + round(aw(2) / fd), nF);
    w = trace(f0:f1);
    w(~isfinite(w)) = -Inf;
    [mx, rel] = max(w);
    if isfinite(mx), pk = mx - base; ttp_ms = (f0 + rel - 1 - stim) * fd; end
    seg = trace(f0:f1) - base;
    seg = seg(isfinite(seg));
    if numel(seg) >= 2, a = trapz(seg) * fd; end
end

function comparison = compare_groups(ga, gb)
    % Statistics are run on PER-COVERSLIP aggregates (one value per coverslip),
    % NOT on pooled ROIs. Pooling ROIs is pseudoreplication: ROIs within a
    % coverslip are not independent, so a pooled test reports an inflated n and
    % an artificially small p. Each metric is collapsed to a coverslip mean
    % first, then compared with Mann-Whitney U across coverslips. Note: with a
    % single animal per genotype the coverslip is a technical (within-animal)
    % replicate, so p-values describe this preparation, not the genotype
    % population - that needs >=3 animals/genotype with animal as the unit.
    comparison = struct('test_type', ...
        'Mann-Whitney U on per-coverslip means (n = coverslips)');
    na = ga.name; nb = gb.name;

    metrics = {'peak_amplitude', 'time_to_peak_ms', 'auc', 'auc_sync', 'auc_async', 'tau_ms'};
    for m = 1:numel(metrics)
        f = metrics{m};
        a = coverslip_means(ga, f, true);   % one value per coverslip (successful trials)
        b = coverslip_means(gb, f, true);
        p = ranksum_or_nan(a, b);
        comparison.(f) = p;
        fprintf('  %-16s %s=%.4g (N=%d cs), %s=%.4g (N=%d cs), p=%.4g\n', ...
            f, na, safemean(a), numel(a), nb, safemean(b), numel(b), p);
    end

    % Per-ROI release statistics across trials (active ROIs, per coverslip):
    % P(response), potency (amplitude of successes), efficacy (failures = 0).
    roi_metrics = {'roi_p_response', 'roi_potency', 'roi_efficacy'};
    for m = 1:numel(roi_metrics)
        f = roi_metrics{m};
        a = coverslip_means_roi(ga, f);
        b = coverslip_means_roi(gb, f);
        p = ranksum_or_nan(a, b);
        comparison.(f) = p;
        fprintf('  %-16s %s=%.4g (N=%d cs), %s=%.4g (N=%d cs), p=%.4g\n', ...
            f, na, safemean(a), numel(a), nb, safemean(b), numel(b), p);
    end

    % Trial success rate: per-coverslip fraction of (ROI, trial) units that
    % responded, compared across coverslips (not ROI-pooled: pseudoreplication).
    a = coverslip_resp_frac(ga);
    b = coverslip_resp_frac(gb);
    p = ranksum_or_nan(a, b);
    comparison.responder_fraction_p = p;
    comparison.responder_fraction   = safemean(a) - safemean(b);
    fprintf('  %-16s %s=%.1f%% (N=%d cs), %s=%.1f%% (N=%d cs), p=%.4g\n', ...
        'success_rate', na, 100*safemean(a), numel(a), nb, 100*safemean(b), numel(b), p);

    % Asynchronous-release fraction: per-coverslip async/(sync+async) events,
    % compared across coverslips. This is the direct test of the slow tail seen
    % in the averaged trace (and the post-stim PSTH), which the per-ROI tau / AUC
    % metrics do not capture.
    a = coverslip_async_frac(ga);
    b = coverslip_async_frac(gb);
    p = ranksum_or_nan(a, b);
    comparison.async_fraction_p = p;
    comparison.async_fraction   = safemean(a) - safemean(b);
    fprintf('  %-16s %s=%.1f%% (N=%d cs), %s=%.1f%% (N=%d cs), p=%.4g\n', ...
        'async_frac', na, 100*safemean(a), numel(a), nb, 100*safemean(b), numel(b), p);
end

function p = ranksum_or_nan(a, b)
    p = NaN;
    if numel(a) >= 1 && numel(b) >= 1 && (numel(a) + numel(b)) >= 3
        p = ranksum(a, b);
    end
end

function v = responder_vals(gd, field)
    v = gd.(field)(gd.is_responder);
    v = v(isfinite(v));
end

function m = safemean(x)
    if isempty(x), m = NaN; else, m = mean(x, 'omitnan'); end
end

%% ======================================================================
%% Plots
%% ======================================================================
function fig = plot_oneap_average_traces(groups, config)
    fig = figure('Name', 'OneAP Average Traces', 'Position', [100 100 1000 650]);
    hold on;
    stim_frame = config.oneAP.stim_frame;
    frame_dur  = 1000 / config.frame_rate;
    tax        = make_time_axis(config);

    legend_entries = {};
    for g = 1:numel(groups)
        gd = groups(g);
        if isempty(gd.dfof_all), continue; end
        resp = gd.is_responder;
        if ~any(resp), resp = true(1, size(gd.dfof_all, 2)); end
        M = gd.dfof_all(:, resp);
        n = sum(isfinite(M), 2);
        mu  = mean(M, 2, 'omitnan');
        sem = std(M, 0, 2, 'omitnan') ./ sqrt(max(n, 1));
        t = tax.from_ms(((1:size(M,1))' - stim_frame) * frame_dur);

        fill([t; flipud(t)], [mu+sem; flipud(mu-sem)], gd.color, ...
            'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(t, mu, 'Color', gd.color, 'LineWidth', 2);
        legend_entries{end+1} = sprintf('%s (n=%d responders)', gd.name, sum(resp)); %#ok<AGROW>

        % Overlay the averaged double-exponential decay fit and annotate
        % tau_fast / tau_slow / slow-amplitude fraction in the legend.
        af = gd.avg_fit;
        if isfield(af, 'tau_fast') && isfinite(af.tau_fast)
            plot(tax.from_ms(af.fit_t), af.fit_y, '--', 'Color', gd.color * 0.6, 'LineWidth', 1.5);
            legend_entries{end+1} = sprintf(['%s fit: \\tau_{fast}=%.0f ms, ' ...
                '\\tau_{slow}=%.0f ms, slow=%.0f%%'], gd.name, af.tau_fast, ...
                af.tau_slow, 100*af.slow_frac); %#ok<AGROW>
        end
    end

    yl = ylim;
    plot([0 0], yl, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    text(0, yl(2), ' stim', 'VerticalAlignment', 'top', 'FontSize', 9);
    xlim(tax.from_ms([-200, max(config.oneAP.async_window_ms(2), config.oneAP.auc_window_ms(2)) + 100]));   % default view
    ylabel('\DeltaF/F (mean \pm SEM)');
    title('1AP Average Evoked Trace (double-exponential decay fit)');
    legend(legend_entries, 'Location', 'best');
    grid on; box on;
    setappdata(gca, 'ymin', config.visualization.ymin_mean_trace);   % see apply_ymin_zero
    % Label in the chosen unit; override xlim with config.x_axis.range_frames if set.
    format_time_axis(gca, tax, 'Time from stimulus', stim_frame);
end

function fig = plot_oneap_cs_average(groups, config)
    % Coverslip-average dF/F. One row per group: the first tile is the group
    % average OF the coverslip averages (mean +/- SEM across coverslips, n =
    % coverslips - the statistical unit), the rest are one tile per coverslip
    % (mean +/- SEM over successful trials, dashed = all scorable trials,
    % grey = each region c#). y is shared across all tiles so groups compare.
    fig = [];
    nG = numel(groups);
    nCs = arrayfun(@(g) numel(g.cs_avg), groups);
    if nG == 0 || all(nCs == 0), return; end
    stim_frame = config.oneAP.stim_frame;
    frame_dur  = 1000 / config.frame_rate;
    tax        = make_time_axis(config);
    xwin       = tax.from_ms([-200, max(config.oneAP.async_window_ms(2), config.oneAP.auc_window_ms(2)) + 100]);
    aw         = tax.from_ms(config.oneAP.auc_window_ms);
    nC = max(nCs) + 1;

    fig = figure('Name', 'OneAP Coverslip average dF/F', ...
        'Position', [40 40 max(1200, 260*nC) max(450, 280*nG)], 'Color', 'w');
    tl = tiledlayout(fig, nG, nC, 'TileSpacing', 'tight', 'Padding', 'tight');
    axs = gobjects(1, 0);
    for g = 1:nG
        gd = groups(g);
        if nCs(g) == 0, continue; end
        t = tax.from_ms(((1:numel(gd.cs_avg(1).trace)) - stim_frame) * frame_dur);

        % Group tile: each coverslip's successful-trial average (thin) + their
        % mean +/- SEM; dashed = mean of the all-scorable coverslip averages
        ax = nexttile(tl, (g-1)*nC + 1); hold on; axs(end+1) = ax; %#ok<AGROW>
        C  = vertcat(gd.cs_avg.trace_resp);            % [coverslips x frames]
        mu = mean(C, 1, 'omitnan');
        se = std(C, 0, 1, 'omitnan') ./ sqrt(max(sum(isfinite(C), 1), 1));
        xregion(aw(1), aw(2), 'FaceColor', [0.45 0.45 0.50], 'FaceAlpha', 0.10);
        fill([t, fliplr(t)], [mu + se, fliplr(mu - se)], gd.color, ...
            'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(t, C', '-', 'Color', [gd.color 0.45], 'LineWidth', 0.7);
        plot(t, mean(vertcat(gd.cs_avg.trace), 1, 'omitnan'), '--', 'Color', gd.color, 'LineWidth', 1);
        plot(t, mu, '-', 'Color', gd.color, 'LineWidth', 2);
        xline(0, 'k--');
        [pk, ~, a] = avg_trace_metrics(mu, config);
        title(sprintf('%s: mean of %d coverslips\npeak=%.3f AUC=%.2f', ...
            gd.name, size(C, 1), pk, a), 'FontSize', 8, 'Interpreter', 'none');
        ylabel(sprintf('%s\n\\DeltaF/F', gd.name), 'FontSize', 8, 'FontWeight', 'bold');
        style_cs_tile(ax, tax, xwin, g == nG, config);

        % One tile per coverslip
        for c = 1:nCs(g)
            e  = gd.cs_avg(c);
            ax = nexttile(tl, (g-1)*nC + 1 + c); hold on; axs(end+1) = ax; %#ok<AGROW>
            xregion(aw(1), aw(2), 'FaceColor', [0.45 0.45 0.50], 'FaceAlpha', 0.10);
            fill([t, fliplr(t)], [e.trace_resp + e.sem_resp, fliplr(e.trace_resp - e.sem_resp)], gd.color, ...
                'FaceAlpha', 0.2, 'EdgeColor', 'none');
            hr = gobjects(0);
            for r = 1:numel(e.regions)
                hr = plot(t, e.regions(r).trace, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
            end
            hs = plot(t, e.trace, '--', 'Color', gd.color, 'LineWidth', 1);
            hm = plot(t, e.trace_resp, '-', 'Color', gd.color, 'LineWidth', 1.5);
            xline(0, 'k--');
            if c == 1   % spell out the three line types once per row
                h = [hm, hs, hr];
                names = {'coverslip mean, successful trials only', ...
                    'coverslip mean, all scorable trials', ...
                    'each region (c#), all scorable trials'};
                legend(h, names(1:numel(h)), 'Location', 'northeast', 'FontSize', 7, 'Box', 'off');
            end
            title(sprintf('%s  peak=%.3f AUC=%.2f\n%d scorable ROIs, %d/%d successful trials, %d c', ...
                e.cs, e.peak_resp, e.auc_resp, e.n_rois, e.n_resp, e.n_units, numel(e.regions)), ...
                'FontSize', 8, 'Interpreter', 'none');
            style_cs_tile(ax, tax, xwin, g == nG, config);
        end
    end
    share_ylim(axs, config);
    title(tl, sprintf(['1AP Coverslip average \\DeltaF/F\n' ...
        'solid = mean \\pm SEM over successful trials only (peak/AUC from this); dashed = all scorable trials (successes + failures); ' ...
        'grey = each region c# (all trials)\n' ...
        'First column = mean of coverslip averages (n = coverslips); shaded = AUC / peak window ' ...
        '%g-%g ms; x in %s'], config.oneAP.auc_window_ms, tax.unit), ...
        'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'tex');
end

function style_cs_tile(ax, tax, xwin, isBottom, config)
    set(ax, 'FontSize', 7); grid(ax, 'on'); box(ax, 'on');
    ax.YAxis.Exponent = 0;
    scroll_time_axis(ax, xwin);   % stim window shown; drag to scroll x
    setappdata(ax, 'ymin', config.visualization.ymin_mean_trace);   % see apply_ymin_zero
    if isBottom, xlabel(ax, sprintf('Time from stimulus (%s)', tax.unit), 'FontSize', 7); end
end

function af = fit_group_double(gd, config)
    % Group-averaged responder decay, fit with a DOUBLE exponential from the
    % averaged peak to the END of the trace (no fixed window). Returns the time
    % constants, the slow-amplitude fraction (asynchronous-release index), the
    % fit R^2 and the fitted curve (fit_t in ms from stim, fit_y in dF/F).
    af = struct('tau_fast', NaN, 'tau_slow', NaN, 'slow_frac', NaN, ...
        'r2', NaN, 'n_roi', 0, 'fit_t', [], 'fit_y', []);
    if isempty(gd.dfof_all) || ~any(gd.is_responder), return; end

    stim_frame = config.oneAP.stim_frame;
    frame_dur  = 1000 / config.frame_rate;
    nFrames    = size(gd.dfof_all, 1);

    mu = mean(gd.dfof_all(:, gd.is_responder), 2, 'omitnan');   % [frames x 1]
    af.n_roi = sum(gd.is_responder);

    % Decay onset = peak of the averaged trace within the tau-peak window
    peak_end   = min(stim_frame + round(config.oneAP.tau_peak_window_ms / frame_dur), nFrames);
    seg        = mu(stim_frame:peak_end);
    seg(~isfinite(seg)) = -Inf;
    [~, rel]   = max(seg);
    onset      = stim_frame + rel - 1;
    if nFrames <= onset + 4, return; end

    dec_idx = onset:nFrames;
    y = mu(dec_idx);
    tt = (0:numel(dec_idx)-1)' * frame_dur;          % ms from peak
    [tf, ts, Af, As, Cfit, r2] = fit_decay_double(tt, y);
    af.r2 = r2;
    if isfinite(tf) && isfinite(ts)
        af.tau_fast  = tf;
        af.tau_slow  = ts;
        af.slow_frac = As / max(Af + As, eps);
        af.fit_t = (dec_idx(:) - stim_frame) * frame_dur;               % ms from stim (col)
        af.fit_y = Af*exp(-tt/tf) + As*exp(-tt/ts) + Cfit;              % col
    end
end

function fig = plot_oneap_amplitude(groups, comparison)
    % Evoked peak amplitude: per-coverslip SuperPlot (the statistic) plus the
    % distributions behind it (ROIs pooled across coverslips: shape only).
    % Left: SuperPlot; middle: every successful trial; right: per-ROI potency.
    fig = figure('Name', 'OneAP Evoked Peak Amplitude', 'Position', [150 120 1700 560]);
    tl  = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, '1AP Evoked Peak Amplitude', 'FontSize', 14, 'FontWeight', 'bold');
    subtitle(tl, 'left: per-coverslip SuperPlot (statistic)  |  histograms: pooled across coverslips (shape only)', ...
        'FontSize', 9);
    cols = vertcat(groups.color);
    succ = arrayfun(@(gd) gd.peak_amplitude(gd.is_responder), groups, 'UniformOutput', false);
    pot  = arrayfun(@(gd) gd.roi_potency(gd.roi_active), groups, 'UniformOutput', false);

    plot_group_box(nexttile(tl), groups, 'peak_amplitude', true, 'Evoked peak amplitude (\DeltaF/F)', ...
        'Per-coverslip means (successful trials)', comparison, 'peak_amplitude');

    ax = nexttile(tl); hold(ax, 'on');
    grouped_hist(ax, succ, 40, cols, n_labels(groups, succ));
    xlabel(ax, 'Evoked peak amplitude (\DeltaF/F)'); ylabel(ax, 'Relative fraction of trials');
    title(ax, 'All successful trials'); legend(ax, 'Location', 'northeast', 'Interpreter', 'none');
    grid(ax, 'on'); box(ax, 'on');

    ax = nexttile(tl); hold(ax, 'on');
    grouped_hist(ax, pot, 30, cols, n_labels(groups, pot));
    xlabel(ax, 'Potency per ROI (\DeltaF/F)'); ylabel(ax, 'Relative fraction of active ROIs');
    title(ax, 'Per-ROI potency (mean of successes)');
    legend(ax, 'Location', 'northeast', 'Interpreter', 'none'); grid(ax, 'on'); box(ax, 'on');
end

function labels = n_labels(groups, vals)
    % Legend names 'group (n=#finite values)'.
    labels = cellfun(@(n, x) sprintf('%s (n=%d)', n, sum(isfinite(x))), ...
        {groups.name}, vals, 'UniformOutput', false);
end

function plot_metric_panels(fig, groups, metrics, comparison, figTitle)
    % One column per metric: per-coverslip SuperPlot on top, side-by-side
    % relative-frequency histogram of the successful-trial values below
    % (pooled across coverslips: distribution shape only, not the statistic).
    % metrics: rows of {field, axis label, title}.
    nM   = size(metrics, 1);
    cols = vertcat(groups.color);
    tl   = tiledlayout(fig, 2, nM, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, figTitle, 'FontSize', 13, 'FontWeight', 'bold');
    subtitle(tl, 'top: per-coverslip SuperPlot (statistic)  |  bottom: successful trials pooled across coverslips (shape only; x clipped to 0.5-99.5th percentile)', ...
        'FontSize', 9);
    for k = 1:nM
        [f, ylab, ttl] = metrics{k, :};
        ax = nexttile(tl, k);
        plot_group_box(ax, groups, f, true, ylab, ttl, comparison, f);
        ax = nexttile(tl, nM + k); hold(ax, 'on'); axes(ax);
        v = arrayfun(@(gd) gd.(f)(gd.is_responder), groups, 'UniformOutput', false);
        allv = [v{:}]; allv = allv(isfinite(allv));
        if ~isempty(allv)
            lim = prctile(allv, [0.5 99.5]);   % a few outliers don't set the range
            if lim(2) <= lim(1), lim = lim(1) + [-0.5 0.5]; end
            grouped_hist(gca, v, linspace(lim(1), lim(2), 31), cols, n_labels(groups, v));
            xlim(lim);
        end
        xlabel(ylab); ylabel('Relative fraction of trials');
        title([ttl ' distribution']);
        legend('Location', 'northeast', 'Interpreter', 'none'); grid on; box on;
    end
end

function fig = plot_oneap_kinetics(groups, comparison)
    fig = figure('Name', 'OneAP Kinetics', 'Position', [200 60 1500 850]);
    plot_metric_panels(fig, groups, { ...
        'time_to_peak_ms', 'Time to peak (ms)',      'Time to Peak'; ...
        'auc',             'AUC (\DeltaF/F\cdotms)', 'Area Under Curve'; ...
        'tau_ms',          '\tau decay (ms)',        'Decay \tau'}, comparison, ...
        '1AP Response Kinetics (responders)');
end

function plot_group_box(ax, groups, field, useResponders, ylab, ttl, comparison, pfield)
    % SuperPlot: faint ROI cloud (shows spread) PLUS the per-coverslip means
    % (large black-edged dots) - the latter are the actual statistical units.
    % The box and the annotated p-value both describe the coverslip level, so
    % the figure can't imply more independence than the test assumes.
    axes(ax); hold on;
    nG = numel(groups);
    cs_vals = []; labels = cell(1, nG);
    for g = 1:nG
        gd = groups(g);

        % Per-coverslip box/whisker (native, drawn behind the points)
        cm = coverslip_means(gd, field, useResponders);
        draw_box(ax, g, cm, gd.color, 0.5);

        % Faint ROI-level points (descriptive only)
        v = gd.(field);
        if useResponders, v = v(gd.is_responder); end
        v = v(isfinite(v));
        if ~isempty(v)
            jit = (rand(1, numel(v)) - 0.5) * 0.3;
            scatter(g + jit, v, 10, gd.color, 'filled', 'MarkerFaceAlpha', 0.18);
        end

        % Per-coverslip means: the box and stats operate on THESE
        if ~isempty(cm)
            cjit = (rand(1, numel(cm)) - 0.5) * 0.18;
            scatter(g + cjit, cm, 55, gd.color, 'filled', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 1);
        end
        cs_vals = [cs_vals, cm]; %#ok<AGROW>
        labels{g} = sprintf('N=%d cs, %d ROIs', numel(cm), numel(v));
    end
    set_group_ticks(ax, {groups.name}, labels);
    xlim([0.5, nG + 0.5]);
    ylabel(ylab); title(ttl); grid on; box on;
    annotate_pvalue(comparison, pfield, cs_vals);
end

function annotate_pvalue(comparison, pfield, all_vals)
    all_vals = all_vals(isfinite(all_vals));   % drop NaN padding from empty groups
    if isempty(comparison) || ~isfield(comparison, pfield) || isempty(all_vals), return; end
    p = comparison.(pfield);
    if isnan(p), return; end
    y = max(all_vals) + 0.12 * (max(all_vals) - min(all_vals) + eps);
    if p < 0.001, ptxt = 'p < 0.001 ***';
    elseif p < 0.01, ptxt = sprintf('p = %.3f **', p);
    elseif p < 0.05, ptxt = sprintf('p = %.3f *', p);
    else, ptxt = sprintf('p = %.3f (ns)', p);
    end
    line([1 2], [y y], 'Color', 'k');
    text(1.5, y, ptxt, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
    % Y starts at 0; leave headroom for the p-value text. The span uses the
    % data's distance from 0 so identical values (e.g. every coverslip at
    % P(response) = 1) still give an increasing range.
    span = max(abs([all_vals, y]));
    if span == 0, span = 1; end
    ylim([min(0, min(all_vals)), y + 0.1 * span]);
end

function fig = plot_oneap_auc_split(groups, comparison)
    % AUC split by release timing: total / synchronous / asynchronous, all
    % responder-gated and aggregated per coverslip (same units as the total AUC).
    % The asynchronous panel is the magnitude of the slow tail seen in the
    % averaged trace - the quantity the per-ROI mono tau cannot capture.
    fig = figure('Name', 'OneAP AUC Split', 'Position', [180 60 1500 850]);
    plot_metric_panels(fig, groups, { ...
        'auc',       'AUC (\DeltaF/F\cdotms)', 'Total AUC'; ...
        'auc_sync',  'AUC (\DeltaF/F\cdotms)', 'Synchronous AUC'; ...
        'auc_async', 'AUC (\DeltaF/F\cdotms)', 'Asynchronous AUC (slow tail)'}, comparison, ...
        '1AP Area Under Curve by release timing (responders)');
end

function fig = plot_oneap_sync_async(groups, config)
    % Sync / async release events from the shared 2AP detector (onset-binned,
    % config.twoAP.events). Returns [] when no group has event binning (frames
    % slower than events.max_frame_ms), as the 2AP sync/async figure does.
    fig = [];
    if ~any([groups.events_enabled])
        fprintf('  1AP sync/async figures skipped: frames slower than %g ms (events.max_frame_ms).\n', ...
            config.twoAP.events.max_frame_ms);
        return;
    end
    groups = groups([groups.events_enabled]);
    e   = config.twoAP.events;
    fig = figure('Name', 'OneAP Sync vs Async', 'Position', [150 100 1700 580]);
    tl  = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, '1AP Synchronous / Asynchronous Release', 'FontSize', 14, 'FontWeight', 'bold');
    subtitle(tl, sprintf(['release events as in 2AP: onset \\leq %g ms = sync, later = async ' ...
        '(to %g ms); async spont-subtracted'], e.sync_cutoff_ms, e.async_window_ms), 'FontSize', 9);
    tax = make_time_axis(config);

    % Left: group-pooled sync / async share of all events
    ax = nexttile(tl); hold(ax, 'on');
    nG = numel(groups);
    fracs = zeros(nG, 2);
    for g = 1:nG
        [s, a] = event_totals(groups(g));
        fracs(g, :) = [s, a] / max(s + a, 1);
    end
    hb = bar(ax, fracs, 'stacked');
    hb(1).FaceColor = [0.2 0.6 0.2]; hb(2).FaceColor = [0.8 0.5 0.1];
    set(ax, 'XTick', 1:nG, 'XTickLabel', {groups.name}, 'TickLabelInterpreter', 'none');
    ylabel(ax, 'Fraction of events'); ylim(ax, [0 1]);
    legend(ax, {'Synchronous', 'Asynchronous'}, 'Location', 'southoutside', 'Orientation', 'horizontal');
    title(ax, 'Event share (pooled)'); grid(ax, 'on'); box(ax, 'on');

    % Middle: distribution of the per-trial async share async/(sync+async) over
    % (ROI, trial) units with >= 1 event (raw counts: a per-trial spont
    % correction would be mostly noise); legend = group-pooled share.
    ax = nexttile(tl); hold(ax, 'on');
    vals = cell(1, nG); labels = cell(1, nG);
    for g = 1:nG
        gd = groups(g);
        tot = gd.sync_count + gd.async_count;
        ok  = tot > 0;                                  % NaN (unscored) compares false
        vals{g} = gd.async_count(ok) ./ tot(ok);
        [s, a] = event_totals(gd);
        labels{g} = sprintf('%s (%d trials, pooled %.0f%% async)', gd.name, numel(vals{g}), ...
            100 * a / max(s + a, 1));
    end
    grouped_hist(ax, vals, 0:0.1:1, vertcat(groups.color), labels);
    xlim(ax, [0 1]); xlabel(ax, 'Asynchronous share of events per trial (raw)');
    ylabel(ax, 'Fraction of trials with events');
    legend(ax, 'Location', 'north', 'Interpreter', 'none');
    title(ax, 'Per-trial asynchronous share'); grid(ax, 'on'); box(ax, 'on');

    % Right: event ONSET latency histogram. Time 0 = the refined stim frame
    % (Stim_Frames_Refined), the same zero every 1AP metric uses.
    % Side-by-side bars per bin, each group normalised to its own event total
    % (relative frequency) so groups with different ROI counts compare directly.
    ax = nexttile(tl); hold(ax, 'on'); axes(ax);
    has = arrayfun(@(g) ~isempty(g.temporal_hist_sum), groups);
    if any(has)
        ctr = groups(find(has, 1)).bin_centers_ms(:);
        Y = zeros(numel(ctr), nG);
        for g = find(has)
            h = groups(g).temporal_hist_sum(:);
            Y(:, g) = h / max(sum(h), 1);
        end
        hb = bar(tax.from_ms(ctr), Y(:, has), 1, 'grouped', 'EdgeColor', 'none');
        gi = find(has);
        for k = 1:numel(hb)
            hb(k).FaceColor   = groups(gi(k)).color;
            hb(k).DisplayName = groups(gi(k)).name;
        end
    end
    yl = ylim;
    % Sync region = onsets up to the cutoff (frame-rounded), to the bin's half edge
    sync_hi = tax.from_ms(round(e.sync_cutoff_ms / config.frame_dur_ms) * config.frame_dur_ms ...
        + config.oneAP.bin_size_ms / 2);
    hp = patch([tax.from_ms(-config.oneAP.bin_size_ms / 2) sync_hi sync_hi tax.from_ms(-config.oneAP.bin_size_ms / 2)], ...
        [0 0 yl(2) yl(2)], [0.2 0.6 0.2], 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    uistack(hp, 'bottom');
    ylabel(sprintf('Relative frequency (%d ms bins)', config.oneAP.bin_size_ms));
    title('Event onset latency (green = sync)');
    legend('Location', 'best'); grid on; box on;
    format_time_axis(gca, tax, 'Event onset after stimulus');
end

function fig = plot_oneap_responders(groups, comparison)
    fig = figure('Name', 'OneAP Responder Fraction', 'Position', [300 150 1300 580]);
    tl  = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, '1AP trial success rate', 'FontSize', 13, 'FontWeight', 'bold');
    subtitle(tl, 'ROI x trial units with a validated evoked response', 'FontSize', 9);
    nG = numel(groups);
    vals = cell(1, nG); labels = cell(1, nG); fr = zeros(1, nG);
    for g = 1:nG
        gd = groups(g);
        vals{g} = 100 * coverslip_resp_frac(gd);
        fr(g) = 100 * gd.n_responders / max(gd.total_rois, 1);
        labels{g} = sprintf('%s (N=%d cs; pooled %d/%d = %.0f%%)', gd.name, numel(vals{g}), ...
            gd.n_responders, gd.total_rois, fr(g));
    end

    % Left: group-pooled rate (bar) + per-coverslip rates (dots)
    draw_fraction_bars(nexttile(tl), groups, fr, vals, 'Successful trials (%)', ...
        'Pooled bar, coverslip dots', comparison, 'responder_fraction_p');

    % Right: histogram of per-coverslip success rates (the statistical units)
    ax = nexttile(tl); hold(ax, 'on');
    grouped_hist(ax, vals, 0:10:100, vertcat(groups.color), labels);
    xlim(ax, [0 100]); xlabel(ax, 'Successful trials per coverslip (%)');
    ylabel(ax, 'Fraction of coverslips');
    legend(ax, 'Location', 'northeast', 'Interpreter', 'none');
    grid(ax, 'on'); box(ax, 'on');
    title(ax, 'Distribution of per-coverslip success rates');
end

function draw_fraction_bars(ax, groups, pooled, cs_vals, ylab, ttl, comparison, pfield)
    % Bar = group-pooled percentage; dots = per-coverslip percentages (the
    % statistical units); bracket = Mann-Whitney p between groups 1 and 2.
    axes(ax); hold(ax, 'on');
    nG = numel(groups);
    hb = bar(pooled, 0.6, 'FaceColor', 'flat');
    for g = 1:nG, hb.CData(g, :) = groups(g).color; end
    labels = cell(1, nG);
    for g = 1:nG
        cf = cs_vals{g};
        if ~isempty(cf)
            scatter(g + (rand(1, numel(cf)) - 0.5) * 0.25, cf, 45, groups(g).color, ...
                'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
        end
        labels{g} = sprintf('pooled %.0f%%, N=%d cs', pooled(g), numel(cf));
    end
    set_group_ticks(ax, {groups.name}, labels);
    xlim([0.4, nG + 0.6]);
    ylabel(ylab); title(ttl); grid on; box on;
    top = max([pooled, cs_vals{:}, 1]);
    ylim([0, top * 1.25]);   % room for the p-value bracket above 100%
    pt = pvalue_text(comparison, pfield, nG);
    if ~isempty(pt)
        y = top * 1.1;
        line([1 2], [y y], 'Color', 'k');
        text(1.5, y, pt, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontWeight', 'bold');
    end
end

function set_group_ticks(ax, names, details)
    % x ticks = group names only; the per-group detail (n, pooled value) goes in
    % the x label, one line per group. Axes layout reserves room for a
    % multi-line x label, but not for multi-line tick labels (those overlapped
    % the panel below and were clipped).
    set(ax, 'XTick', 1:numel(names), 'XTickLabel', names, 'TickLabelInterpreter', 'none', ...
        'XTickLabelRotation', 0);
    xlabel(ax, cellfun(@(n, d) sprintf('%s: %s', n, d), names, details, 'UniformOutput', false), ...
        'Interpreter', 'none', 'FontSize', 8);
end

function t = pvalue_text(comparison, pfield, nG)
    % One-line Mann-Whitney p-value for a figure title ('' when not available).
    t = '';
    if nG >= 2 && ~isempty(comparison) && isfield(comparison, pfield) ...
            && ~isnan(comparison.(pfield))
        p = comparison.(pfield);
        if p < 0.001, t = 'p < 0.001 *** (per coverslip)';
        else, t = sprintf('p = %.3f (per coverslip)', p); end
    end
end

function fig = plot_oneap_release_prob(groups, comparison)
    % Per-ROI release statistics across trials, over ACTIVE ROIs (>= 1 success):
    %   (1) distribution of P(response) = successes / valid trials
    %   (2) P(response), (3) potency, (4) efficacy - per-coverslip SuperPlots
    fig = figure('Name', 'OneAP Release Probability', 'Position', [140 60 1600 860]);
    tl  = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, '1AP release reliability across trials', 'FontSize', 14, 'FontWeight', 'bold');
    nG = numel(groups);

    ax = nexttile(tl); hold(ax, 'on'); axes(ax);
    labels = cell(1, nG); vals = cell(1, nG);
    maxT = max([1, arrayfun(@(g) max([g.roi_n_valid, 1]), groups)]);
    for g = 1:nG
        gd = groups(g);
        v = gd.roi_p_response(gd.roi_active & isfinite(gd.roi_p_response));
        labels{g} = sprintf('%s (%d active / %d ROIs)', gd.name, numel(v), gd.n_rois);
        vals{g} = v;
    end
    grouped_hist(gca, vals, ((0:maxT+1) - 0.5) / maxT, ...   % centred on k/maxT
        vertcat(groups.color), labels);
    xlim([0 1.05]); xlabel('P(response) per ROI'); ylabel('Fraction of active ROIs');
    title('Per-ROI response probability across trials');
    legend('Location', 'northwest', 'Interpreter', 'none'); grid on; box on;

    plot_roi_box(nexttile(tl), groups, 'roi_p_response', 'P(response)', ...
        'Response probability (active ROIs)', comparison);
    ylim([0 1.1]);
    plot_roi_box(nexttile(tl), groups, 'roi_potency', '\DeltaF/F', ...
        'Potency (mean amplitude of successes)', comparison);
    plot_roi_box(nexttile(tl), groups, 'roi_efficacy', '\DeltaF/F', ...
        'Efficacy (mean amplitude, failures = 0)', comparison);
end

function plot_roi_box(ax, groups, field, ylab, ttl, comparison)
    % SuperPlot of a per-ROI statistic: faint active-ROI cloud + per-coverslip
    % means (the statistical units), as plot_group_box does for trial units.
    axes(ax); hold on;
    nG = numel(groups);
    cs_vals = []; labels = cell(1, nG);
    for g = 1:nG
        gd = groups(g);
        cm = coverslip_means_roi(gd, field);
        draw_box(ax, g, cm, gd.color, 0.5);
        v = gd.(field)(gd.roi_active);
        v = v(isfinite(v));
        if ~isempty(v)
            scatter(g + (rand(1, numel(v)) - 0.5) * 0.3, v, 10, gd.color, ...
                'filled', 'MarkerFaceAlpha', 0.18);
        end
        if ~isempty(cm)
            scatter(g + (rand(1, numel(cm)) - 0.5) * 0.18, cm, 55, gd.color, ...
                'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
        end
        cs_vals = [cs_vals, cm]; %#ok<AGROW>
        labels{g} = sprintf('N=%d cs, %d ROIs', numel(cm), numel(v));
    end
    set_group_ticks(ax, {groups.name}, labels);
    xlim([0.5, nG + 0.5]);
    ylabel(ylab); title(ttl); grid on; box on;
    annotate_pvalue(comparison, field, cs_vals);
end

function fig = plot_oneap_async_fraction(groups, comparison, config)
    % Per-coverslip ASYNCHRONOUS-release fraction: async/(sync+async) release
    % events (shared 2AP detector, async spont-subtracted - event_totals).
    % Left: group-pooled fraction (bar, matches the onset histogram) + per-
    % coverslip fractions (dots, the statistical units); p-value = Mann-Whitney
    % U across coverslips. Right: histogram of the per-coverslip fractions.
    % Returns [] when event binning did not run (slow frames).
    fig = [];
    if ~any([groups.events_enabled]), return; end
    e   = config.twoAP.events;
    fig = figure('Name', 'OneAP Asynchronous Fraction', 'Position', [320 160 1300 580]);
    tl  = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, '1AP asynchronous release fraction', 'FontSize', 13, 'FontWeight', 'bold');
    subtitle(tl, sprintf(['async / all release events; onset > %g ms (to %g ms) = async, ' ...
        'as in 2AP; async spont-subtracted'], e.sync_cutoff_ms, e.async_window_ms), 'FontSize', 9);
    nG = numel(groups);
    vals = cell(1, nG); labels = cell(1, nG); fr = zeros(1, nG);
    for g = 1:nG
        gd = groups(g);
        [s, a] = event_totals(gd);
        vals{g} = 100 * coverslip_async_frac(gd);
        fr(g) = 100 * a / max(s + a, 1);
        labels{g} = sprintf('%s (N=%d cs; pooled %.0f%%)', gd.name, numel(vals{g}), fr(g));
    end

    draw_fraction_bars(nexttile(tl), groups, fr, vals, 'Asynchronous fraction (%)', ...
        'Pooled bar, coverslip dots', comparison, 'async_fraction_p');

    ax = nexttile(tl); hold(ax, 'on');
    grouped_hist(ax, vals, 0:5:100, vertcat(groups.color), labels);
    top = max([vals{:}, 10]);
    xlim(ax, [0, min(100, 5 * ceil(top / 5) + 5)]);
    xlabel(ax, 'Asynchronous fraction per coverslip (%)');
    ylabel(ax, 'Fraction of coverslips');
    legend(ax, 'Location', 'northeast', 'Interpreter', 'none');
    grid(ax, 'on'); box(ax, 'on');
    title(ax, 'Distribution of per-coverslip async fractions');
end

function fig = plot_oneap_qc(groups, config)
    fig = figure('Name', 'OneAP QC', 'Position', [120 80 1300 850]);
    tax = make_time_axis(config);

    % Distributions are side-by-side relative-frequency histograms.
    cols  = vertcat(groups.color);
    names = {groups.name};

    % (1) pre-stim noise SD distribution
    subplot(2,2,1); hold on;
    grouped_hist(gca, {groups.prestim_sd}, 30, cols, names);
    xlabel('Pre-stim noise SD (\DeltaF/F)'); ylabel('Fraction of trials');
    title('Pre-stimulus noise (threshold basis)'); legend('Location','best'); grid on; box on;

    % (2) peak-time histogram (verify peak ~ expected), relative to stim
    subplot(2,2,2); hold on;
    binw = tax.from_frames_rel(1);   % one-frame bins in display units
    vals = arrayfun(@(gd) tax.from_frames_rel(gd.peak_frame(gd.is_responder) ...
        - config.oneAP.stim_frame), groups, 'UniformOutput', false);   % time from stim
    allv = [vals{:}]; allv = allv(isfinite(allv));
    if ~isempty(allv)
        edges = ((floor(min(allv) / binw) - 0.5) : (ceil(max(allv) / binw) + 0.5)) * binw;
        grouped_hist(gca, vals, edges, cols, names);
    end
    yl = ylim;
    plot([0 0], yl, 'k--', 'HandleVisibility', 'off');
    text(0, yl(2), ' stim', 'VerticalAlignment', 'top', 'FontSize', 8);
    ylabel('Fraction of responders');
    title('Evoked peak time (responders)'); legend('Location','best'); grid on; box on;
    format_time_axis(gca, tax, 'Peak time from stimulus');

    % (3) tau fit R^2 distribution
    subplot(2,2,3); hold on;
    grouped_hist(gca, {groups.tau_r2}, 0:0.05:1, cols, names);
    xlabel('Decay fit R^2'); ylabel('Fraction of fits');
    title('\tau fit quality'); legend('Location','best'); grid on; box on;

    % (4) tau vs amplitude scatter (responders)
    subplot(2,2,4); hold on;
    for g = 1:numel(groups)
        gd = groups(g);
        m = gd.is_responder & isfinite(gd.tau_ms) & isfinite(gd.peak_amplitude);
        if ~any(m), continue; end
        scatter(gd.peak_amplitude(m), gd.tau_ms(m), 16, gd.color, 'filled', ...
            'MarkerFaceAlpha', 0.4, 'DisplayName', gd.name);
    end
    xlabel('Peak amplitude (\DeltaF/F)'); ylabel('\tau decay (ms)');
    title('Decay \tau vs peak amplitude'); legend('Location','best'); grid on; box on;

    sgtitle('1AP Quality Control', 'FontSize', 14, 'FontWeight', 'bold');
end

function figs = plot_oneap_montage(groups, config)
    % Per-group montage of example ROI traces (responders first) with the
    % pre-stim threshold, stim marker and detected peak.
    figs = [];
    stim_frame = config.oneAP.stim_frame;
    frame_dur  = 1000 / config.frame_rate;
    tax        = make_time_axis(config);
    mc         = config.visualization.montage;
    maxROIs    = mc.oneAP_max_rois;
    ylim_fixed = mc.oneAP_ylim;   % shared y-axis for all panels (dF/F; [] = auto)

    % Show a generous pre-stim segment so the resting noise is visible against
    % the evoked response (the response itself is fast: ~0-30 ms post-stim).
    win_pre  = round(mc.oneAP_pre_ms / frame_dur);
    win_post = round((config.oneAP.async_window_ms(2) + mc.oneAP_post_extra_ms) / frame_dur);

    for g = 1:numel(groups)
        gd = groups(g);
        if isempty(gd.dfof_all), continue; end

        % Prioritise responders, then fill with non-responders
        order = [find(gd.is_responder), find(~gd.is_responder)];
        order = order(order <= size(gd.dfof_all, 2));
        order = order(1:min(maxROIs, numel(order)));
        if isempty(order), continue; end

        nF_full = size(gd.dfof_all, 1);
        [f0, f1] = montage_window(tax, stim_frame, win_pre, win_post, nF_full);
        fidx  = f0:f1;
        t     = tax.from_ms((fidx - stim_frame) * frame_dur);
        tfull = tax.from_ms(((1:nF_full) - stim_frame) * frame_dur);  % FULL-trace axis
        xwin  = [t(1) t(end)];                                        % initial x-window

        % tiledlayout with tight spacing/padding gives each panel more
        % horizontal room than subplot's default margins. With montage.rows set,
        % ROIs beyond one rows x columns grid continue on further figures.
        [nr, nc, per, nP] = montage_grid(numel(order), mc);
        for pg = 1:nP
            pidx = order((pg-1)*per + 1 : min(pg*per, numel(order)));
            ptxt = page_text(pg, nP);
            fig = figure('Name', sprintf('OneAP Montage %s%s', gd.name, ptxt), ...
                'Position', [50 40 1850 1000], 'Color', 'w');
            figs = [figs, fig]; %#ok<AGROW>
            tl = tiledlayout(fig, nr, nc, 'TileSpacing', 'tight', 'Padding', 'tight');

            for k = 1:numel(pidx)
                roi = pidx(k);
                nexttile(tl, k); hold on;
                % Plot the WHOLE trace; the view starts on the stim window and the
                % rest is a horizontal drag away (see scroll_time_axis).
                plot(tfull, gd.dfof_all(:, roi), 'k-', 'LineWidth', 1);
                if isfinite(gd.upper_threshold(roi))
                    yline(gd.upper_threshold(roi), 'g:', 'LineWidth', 1);
                end
                if ~isempty(ylim_fixed), ylim(montage_ylim(ylim_fixed, gd.dfof_all(fidx, roi))); end
                xline(0, '--', 'Color', 'k', 'LineWidth', 1);
                % All counted peaks (whole trace): blue = synchronous, red = asynchronous
                plot_peak_markers(gd, roi, stim_frame, frame_dur, 1, nF_full, ...
                    config.oneAP.sync_window_ms, config.oneAP.async_window_ms, tax);
                % Fitted mono-exponential decay (magenta) where a valid tau exists
                plot_tau_fit(gd, roi, stim_frame, frame_dur, nF_full, tax);
                lbl = unit_label(gd, roi);
                if gd.is_responder(roi) && isfinite(gd.tau_ms(roi))
                    ttl = sprintf('%s R \\tau=%.0fms', lbl, gd.tau_ms(roi));
                elseif gd.is_responder(roi)
                    ttl = sprintf('%s R', lbl);
                else
                    ttl = sprintf('%s -', lbl);
                end
                title(ttl, 'FontSize', 8, 'Interpreter', 'tex');
                set(gca, 'FontSize', 7); grid on; box on;
                scroll_time_axis(gca, xwin);   % show stim window; drag to scroll x
                setappdata(gca, 'ymin', config.visualization.ymin_montage);   % see apply_ymin_zero
            end
            title(tl, sprintf(['1AP ROI montage: %s%s (responders first; drag to scroll x, whole trace plotted; ' ...
                'events at peak: blue=sync onset, red=async onset, magenta=\\tau fit; x in %s)'], gd.name, ptxt, tax.unit), ...
                'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'tex');
        end
    end
end

function figs = plot_oneap_marginal(groups, config)
    % Montage of the MARGINAL synchronous responders - the weakest evoked calls,
    % the ones sitting just above the detection threshold. These are where a
    % responder/non-responder mis-call is most likely, so inspecting them by eye
    % is the check on whether the filter is keeping noise. Same panel grid
    % format and styling as plot_oneap_montage (extended baseline, green upper /
    % cyan lower threshold, stim marker, peak markers, tau fit).
    figs = [];
    stim_frame = config.oneAP.stim_frame;
    frame_dur  = 1000 / config.frame_rate;
    tax        = make_time_axis(config);
    mc         = config.visualization.montage;
    maxROIs    = mc.oneAP_marginal_max_rois;
    ylim_fixed = mc.oneAP_ylim;      % shared y-axis for all panels (dF/F; [] = auto)

    win_pre  = round(mc.oneAP_pre_ms / frame_dur);   % generous pre-stim baseline (see noise)
    win_post = round((config.oneAP.async_window_ms(2) + mc.oneAP_post_extra_ms) / frame_dur);

    for g = 1:numel(groups)
        gd = groups(g);
        if isempty(gd.dfof_all), continue; end
        nCol = size(gd.dfof_all, 2);

        % Marginal = responders ranked by ascending SNR of their LARGEST evoked
        % peak (peak amplitude / pre-stim noise); show the weakest 20 - the
        % borderline calls where a mis-classification is most likely.
        resp = find(gd.is_responder);
        resp = resp(resp <= nCol);
        snr = gd.peak_amplitude(resp) ./ gd.prestim_sd(resp);
        [~, ord] = sort(snr, 'ascend');
        order = resp(ord);
        order = order(1:min(maxROIs, numel(order)));
        if isempty(order), continue; end

        nF_full = size(gd.dfof_all, 1);
        [f0, f1] = montage_window(tax, stim_frame, win_pre, win_post, nF_full);
        fidx  = f0:f1;
        t     = tax.from_ms((fidx - stim_frame) * frame_dur);
        tfull = tax.from_ms(((1:nF_full) - stim_frame) * frame_dur);  % FULL-trace axis
        xwin  = [t(1) t(end)];                                        % initial x-window

        [nr, nc, per, nP] = montage_grid(numel(order), mc);
        for pg = 1:nP
            pidx = order((pg-1)*per + 1 : min(pg*per, numel(order)));
            ptxt = page_text(pg, nP);
            fig = figure('Name', sprintf('OneAP Marginal responders %s%s', gd.name, ptxt), ...
                'Position', [50 40 1850 1000], 'Color', 'w');
            figs = [figs, fig]; %#ok<AGROW>
            tl = tiledlayout(fig, nr, nc, 'TileSpacing', 'tight', 'Padding', 'tight');

            for k = 1:numel(pidx)
                roi = pidx(k);
                nexttile(tl, k); hold on;
                % Whole trace plotted; view starts on the stim window, drag to scroll.
                plot(tfull, gd.dfof_all(:, roi), 'k-', 'LineWidth', 1);
                if isfinite(gd.upper_threshold(roi))
                    yline(gd.upper_threshold(roi), 'g:', 'LineWidth', 1);
                end
                if isfinite(gd.lower_threshold(roi))
                    yline(gd.lower_threshold(roi), 'c:', 'LineWidth', 1);
                end
                if ~isempty(ylim_fixed), ylim(montage_ylim(ylim_fixed, gd.dfof_all(fidx, roi))); end
                xline(0, '--', 'Color', 'k', 'LineWidth', 1);
                plot_peak_markers(gd, roi, stim_frame, frame_dur, 1, nF_full, ...
                    config.oneAP.sync_window_ms, config.oneAP.async_window_ms, tax);
                plot_tau_fit(gd, roi, stim_frame, frame_dur, nF_full, tax);
                snr_roi = gd.peak_amplitude(roi) / gd.prestim_sd(roi);
                title(sprintf('%s SNR=%.1f', unit_label(gd, roi), snr_roi), ...
                    'FontSize', 8, 'Interpreter', 'tex');
                set(gca, 'FontSize', 7); grid on; box on;
                scroll_time_axis(gca, xwin);   % show stim window; drag to scroll x
                setappdata(gca, 'ymin', config.visualization.ymin_montage);   % see apply_ymin_zero
            end
            title(tl, sprintf(['1AP marginal responders: %s%s (weakest evoked calls, ' ...
                'low SNR first; drag to scroll x; green=%g\\sigma, cyan=%g\\sigma, ' ...
                'events at peak: blue=sync onset, red=async onset, magenta=\\tau fit; x in %s)'], gd.name, ptxt, ...
                config.oneAP.upper_threshold_sigma, config.oneAP.lower_threshold_sigma, tax.unit), ...
                'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'tex');
        end
    end
end

function yl = montage_ylim(ylim_fixed, seg)
    % Shared montage y-range, with the top raised (10% headroom) when the trace
    % in the initial view window goes above it, so large responses aren't cut.
    seg = seg(isfinite(seg));
    top = ylim_fixed(2);
    if ~isempty(seg), top = max(top, 1.1 * max(seg)); end
    yl = [ylim_fixed(1), top];
end

function s = page_text(pg, nP)
    % ' (page 2/3)' suffix for multi-page montages; '' for a single page.
    if nP > 1, s = sprintf(' (page %d/%d)', pg, nP); else, s = ''; end
end

function s = unit_label(gd, u)
    % "<Cs#-c#> ROI<n> T<trial>" for a (ROI, trial) unit column of the group:
    % recording (see file_ids), ROI = data column in that recording's CSV,
    % T = k-th stim of the recording (Recordings.csv order). Same string as the
    % unit_id column of the metrics table, so montage -> Excel -> raw traces.
    s = sprintf('%s ROI%d T%d', gd.file_ids{gd.file_of_roi(u)}, gd.unit_roi(u), gd.unit_trial(u));
end

function ids = file_ids(filenames)
    % recording_id per file (e.g. 'Cs1-c6-02'); the file stem instead when two
    % files in the group would share an ID.
    ids = cellfun(@recording_id, filenames, 'UniformOutput', false);
    [~, ~, j] = unique(ids);
    dup = accumarray(j(:), 1) > 1;
    for k = find(dup(j(:)))'
        [~, ids{k}] = fileparts(filenames{k});
    end
end

function plot_tau_fit(gd, roi, stim_frame, frame_dur, nFrames, tax)
    % Overlay the fitted mono-exponential decay for a responder ROI, in dF/F
    % units: yhat(t) = tau_offset + tau_A * exp(-t/tau), drawn from the evoked
    % peak over ~5 time constants (enough to show the decay) as a magenta line.
    if roi > numel(gd.tau_ms) || ~isfinite(gd.tau_ms(roi)) ...
            || ~isfinite(gd.tau_peak_frame(roi)) || ~isfinite(gd.tau_A(roi))
        return;
    end
    pf   = gd.tau_peak_frame(roi);
    fend = min(pf + round(5 * gd.tau_ms(roi) / frame_dur), nFrames);
    if fend <= pf, return; end
    fr   = pf:fend;
    tt   = (fr - pf) * frame_dur;                       % ms from peak
    yhat = gd.tau_offset(roi) + gd.tau_A(roi) * exp(-tt / gd.tau_ms(roi));
    x    = tax.from_ms((fr - stim_frame) * frame_dur);  % display units from stim
    plot(x, yhat, 'm-', 'LineWidth', 1.5);
end

function plot_peak_markers(gd, roi, stim_frame, frame_dur, f0, f1, sync_win, async_win, tax)
    % Overlay every release event on a per-ROI trace, at its PEAK, coloured by
    % its ONSET bin (the shared 2AP detector, as counted in sync/async_count):
    %   blue filled circle  = synchronous event (onset <= events.sync_cutoff_ms)
    %   red  filled circle  = asynchronous event
    % Without event binning (slow frames) the responder peaks are shown in
    % black instead. The y value is read back from the aligned dF/F trace.
    % (sync_win / async_win are unused; kept so the callers stay unchanged.)
    if gd.events_enabled && roi <= numel(gd.release_events)
        ev = gd.release_events(roi);
        times = ev.lat_ms; issync = ev.is_sync;
    elseif roi <= numel(gd.event_peak_times_ms)
        times = gd.event_peak_times_ms{roi}; issync = [];
    else
        return;
    end
    for i = 1:numel(times)
        tt = times(i);
        fr = stim_frame + round(tt / frame_dur);
        if fr < f0 || fr > f1, continue; end
        yv = gd.dfof_all(fr, roi);
        if ~isfinite(yv), continue; end
        xx = tax.from_ms(tt);   % display units from stim
        if isempty(issync), c = 'k'; elseif issync(i), c = 'b'; else, c = 'r'; end
        plot(xx, yv, 'o', 'MarkerEdgeColor', c, 'MarkerFaceColor', c, 'MarkerSize', 5);
    end
end

%% ======================================================================
%% Output / utility
%% ======================================================================
function T = build_metrics_table(groups)
    % One row per (ROI, trial) unit. step1_flags = that trial's Step1_Flags.
    genotype = {}; drug = {}; filename = {}; coverslip = {}; roi = []; trial = [];
    flags = {}; responder = []; unit_id = {};
    peak_amp = []; ttp = []; auc = []; auc_s = []; auc_a = []; tau = []; tau_r2 = [];
    sync_c = []; async_c = []; spont_c = []; tot_ev = []; pre_sd = [];

    for g = 1:numel(groups)
        gd = groups(g);
        n = gd.total_rois;
        drug     = [drug;     repmat({gd.drug}, n, 1)]; %#ok<AGROW>
        trial    = [trial;    gd.unit_trial(:)];        %#ok<AGROW>
        flags    = [flags;    gd.unit_flags(:)];        %#ok<AGROW>
        for j = 1:n
            genotype{end+1,1} = gd.genotype; %#ok<AGROW>
            fidx = gd.file_of_roi(j);
            if fidx >= 1 && fidx <= numel(gd.filenames)
                filename{end+1,1} = gd.filenames{fidx}; %#ok<AGROW>
            else
                filename{end+1,1} = ''; %#ok<AGROW>
            end
            unit_id{end+1,1} = unit_label(gd, j); %#ok<AGROW>
            if j <= numel(gd.coverslip_of_roi)
                coverslip{end+1,1} = gd.coverslip_of_roi{j}; %#ok<AGROW>
            else
                coverslip{end+1,1} = ''; %#ok<AGROW>
            end
        end
        roi       = [roi;       gd.unit_roi(:)];          %#ok<AGROW>
        responder = [responder; gd.is_responder(:)];      %#ok<AGROW>
        peak_amp  = [peak_amp;  gd.peak_amplitude(:)];     %#ok<AGROW>
        ttp       = [ttp;       gd.time_to_peak_ms(:)];    %#ok<AGROW>
        auc       = [auc;       gd.auc(:)];                %#ok<AGROW>
        auc_s     = [auc_s;     gd.auc_sync(:)];           %#ok<AGROW>
        auc_a     = [auc_a;     gd.auc_async(:)];          %#ok<AGROW>
        tau       = [tau;       gd.tau_ms(:)];             %#ok<AGROW>
        tau_r2    = [tau_r2;    gd.tau_r2(:)];             %#ok<AGROW>
        sync_c    = [sync_c;    gd.sync_count(:)];         %#ok<AGROW>
        async_c   = [async_c;   gd.async_count(:)];        %#ok<AGROW>
        spont_c   = [spont_c;   gd.n_spont(:)];            %#ok<AGROW>
        tot_ev    = [tot_ev;    gd.total_events(:)];       %#ok<AGROW>
        pre_sd    = [pre_sd;    gd.prestim_sd(:)];         %#ok<AGROW>
    end

    T = table(unit_id, genotype, drug, filename, coverslip, roi, trial, logical(responder), peak_amp, ttp, auc, ...
        auc_s, auc_a, tau, tau_r2, sync_c, async_c, spont_c, tot_ev, pre_sd, flags, ...
        'VariableNames', {'unit_id', 'genotype', 'drug', 'filename', 'coverslip', 'roi', 'trial', 'is_responder', ...
        'peak_amplitude', 'time_to_peak_ms', 'auc', 'auc_sync', 'auc_async', 'tau_ms', 'tau_r2', ...
        'sync_count', 'async_count', 'spont_count', 'total_events', 'prestim_sd', 'step1_flags'});
end

function T = build_roi_table(groups)
    % One row per physical ROI: release statistics across its trials.
    parts = cell(numel(groups), 1);
    for g = 1:numel(groups)
        gd = groups(g);
        n = gd.n_rois;
        roi_id = arrayfun(@(f, r) sprintf('%s ROI%d', gd.file_ids{f}, r), ...
            gd.roi_file(:), gd.roi_index(:), 'UniformOutput', false);
        parts{g} = table(roi_id, repmat({gd.genotype}, n, 1), repmat({gd.drug}, n, 1), ...
            gd.filenames(gd.roi_file)', gd.roi_coverslip(:), gd.roi_index(:), ...
            gd.roi_n_valid(:), gd.roi_n_success(:), gd.roi_p_response(:), ...
            gd.roi_potency(:), gd.roi_efficacy(:), ...
            'VariableNames', {'roi_id', 'genotype', 'drug', 'filename', 'coverslip', 'roi', ...
            'n_valid_trials', 'n_success', 'p_response', 'potency', 'efficacy'});
    end
    T = vertcat(parts{:});
end

function T = build_coverslip_table(groups)
    % One row per coverslip: the per-coverslip values the statistics compare
    % (same definitions as coverslip_means / coverslip_means_roi /
    % coverslip_resp_frac / coverslip_async_frac), plus peak / time-to-peak /
    % AUC of the coverslip's average over ALL scorable trials, failures included
    % (gd.cs_avg.trace; the allscorable_* columns).
    trialF = {'peak_amplitude', 'time_to_peak_ms', 'auc', 'auc_sync', 'auc_async', 'tau_ms'};
    roiF   = {'roi_p_response', 'roi_potency', 'roi_efficacy'};
    rows = {};
    for g = 1:numel(groups)
        gd = groups(g);
        uniq = unique(gd.coverslip_of_roi, 'stable');
        for i = 1:numel(uniq)
            mu = strcmp(gd.coverslip_of_roi, uniq{i});     % (ROI, trial) units
            mr = strcmp(gd.roi_coverslip, uniq{i});        % physical ROIs
            act = mr & gd.roi_active;
            [s, a] = event_totals(gd, mu);   % async spont-subtracted (NaN at slow frames)
            row = {gd.genotype, gd.drug, uniq{i}, numel(unique(gd.region_of_roi(mu))), ...
                nnz(mr), nnz(act), nnz(mu), nnz(mu & gd.is_responder), mean(gd.is_responder(mu))};
            for f = 1:numel(trialF)
                v = gd.(trialF{f})(mu & gd.is_responder);
                row{end+1} = safemean(v(isfinite(v))); %#ok<AGROW>
            end
            for f = 1:numel(roiF)
                v = gd.(roiF{f})(act);
                row{end+1} = safemean(v(isfinite(v))); %#ok<AGROW>
            end
            if s + a > 0, row{end+1} = a / (s + a); %#ok<AGROW>
            else, row{end+1} = NaN; end %#ok<AGROW>
            k = find(strcmp({gd.cs_avg.cs}, uniq{i}), 1);
            if isempty(k), row = [row, {NaN, NaN, NaN}]; %#ok<AGROW>
            else, row = [row, {gd.cs_avg(k).peak, gd.cs_avg(k).ttp_ms, gd.cs_avg(k).auc}]; %#ok<AGROW>
            end
            rows(end+1, :) = row; %#ok<AGROW>
        end
    end
    names = [{'genotype', 'drug', 'coverslip', 'n_regions', 'n_rois', 'n_active_rois', ...
        'n_trials', 'n_success', 'success_rate'}, trialF, ...
        {'p_response', 'potency', 'efficacy', 'async_fraction', ...
        'allscorable_peak', 'allscorable_ttp_ms', 'allscorable_auc'}];
    if isempty(rows)
        T = cell2table(cell(0, numel(names)), 'VariableNames', names);
    else
        T = cell2table(rows, 'VariableNames', names);
    end
end

function print_summary(groups, comparisons)
    fprintf('\n=== 1AP SUMMARY ===\n');
    fprintf('(Statistical unit = coverslip. ROIs and regions/files within a\n');
    fprintf(' coverslip are pseudoreplicates. With 1 animal/genotype these are\n');
    fprintf(' technical replicates - report N animals for a population claim.)\n');
    for g = 1:numel(groups)
        gd = groups(g);
        ra = responder_vals(gd, 'peak_amplitude');
        rt = responder_vals(gd, 'tau_ms');
        act = gd.roi_active;
        fprintf('%s: %d files, %d coverslips, %d ROIs (%d active), %d/%d successful trials (%.0f%%)\n', ...
            gd.name, gd.num_files, gd.num_coverslips, gd.n_rois, sum(act), gd.n_responders, ...
            gd.total_rois, 100*gd.n_responders/max(gd.total_rois,1));
        fprintf('   active ROIs: P(response) %.2f | potency %.4f dF/F | efficacy %.4f dF/F\n', ...
            safemean(gd.roi_p_response(act)), safemean(gd.roi_potency(act)), ...
            safemean(gd.roi_efficacy(act)));
        fprintf('   successful trials: peak %.4f dF/F | time-to-peak %.1f ms | tau %.0f ms\n', ...
            safemean(ra), safemean(responder_vals(gd,'time_to_peak_ms')), safemean(rt));
        af = gd.avg_fit;
        if isfield(af, 'tau_fast') && isfinite(af.tau_fast)
            fprintf('   avg decay (double-exp, n=%d): tau_fast %.0f ms | tau_slow %.0f ms | slow-frac %.0f%% | R^2 %.3f\n', ...
                af.n_roi, af.tau_fast, af.tau_slow, 100*af.slow_frac, af.r2);
        end
        c = gd.qc_category;
        fprintf('   QC gates: %d sync-peak resp | %d async-peak resp | %d subthreshold | %d noisy | %d excluded\n', ...
            sum(c==1), sum(c==2), sum(c==3), sum(c==4), sum(c==0));
    end
    for k = 1:numel(comparisons)
        c = comparisons(k).stats;
        fprintf(['%s vs %s p-values: P(resp)=%.3g, potency=%.3g, efficacy=%.3g, ' ...
            'peak=%.3g, ttp=%.3g, auc=%.3g, tau=%.3g, success=%.3g, async_frac=%.3g\n'], ...
            comparisons(k).group_a, comparisons(k).group_b, c.roi_p_response, ...
            c.roi_potency, c.roi_efficacy, c.peak_amplitude, c.time_to_peak_ms, ...
            c.auc, c.tau_ms, c.responder_fraction_p, c.async_fraction_p);
    end
    fprintf('===================\n');
end

function config = local_config()
    evalc('config = tracenorm_config();');
    config.dataType = '1AP';
end
