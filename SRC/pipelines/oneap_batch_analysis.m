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
    %             .oneap_roi_table (one row per ROI), .processing_info

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
        plot_handles.amplitude      = plot_oneap_metric(groups, 'peak_amplitude', true, ...
            'Evoked Peak Amplitude (\DeltaF/F)', 'Evoked Peak Amplitude (successful trials)', comparison, 'peak_amplitude');
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

    %% Assemble results
    results = struct();
    results.dataType            = '1AP';
    results.groups              = groups;
    results.comparisons         = comparisons;
    results.comparison          = comparison;
    results.plot_handles        = plot_handles;
    results.oneap_metrics_table = oneap_metrics_table;
    results.oneap_roi_table     = oneap_roi_table;
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
    is_responder = []; prestim_sd = []; sync_count = []; async_count = [];
    total_events = []; peak_frame = []; upper_threshold = []; lower_threshold = [];
    qc_category = []; crossed_upper = [];
    tau_A = []; tau_offset = []; tau_peak_frame = [];
    event_peak_times_ms = {};
    file_of_roi = []; filenames = {}; coverslip_of_roi = {};
    temporal_hist_sum = []; bin_centers_ms = [];

    ref_frames = [];

    for i = 1:numel(fileResults)
        fr = fileResults{i};
        if isempty(fr), continue; end
        o = fr.oneap;
        filenames{end+1} = fr.filename; %#ok<AGROW>

        % Experimental unit = coverslip (Cs#). Multiple regions (c#) on one
        % coverslip are pseudoreplicates and must collapse to it before stats.
        csid = parse_coverslip_id(fr.filename);
        coverslip_of_roi = [coverslip_of_roi, repmat({csid}, 1, fr.numUnits)]; %#ok<AGROW>
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
        roi_file       = [roi_file,       repmat(i, 1, fr.numROIs)];      %#ok<AGROW>
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
        total_events    = [total_events,    o.total_events];     %#ok<AGROW>
        peak_frame      = [peak_frame,      o.peak_frame];       %#ok<AGROW>
        upper_threshold = [upper_threshold, o.upper_threshold];  %#ok<AGROW>
        lower_threshold = [lower_threshold, o.lower_threshold];  %#ok<AGROW>
        qc_category     = [qc_category,     o.qc_category];       %#ok<AGROW>
        crossed_upper   = [crossed_upper,   o.crossed_upper];     %#ok<AGROW>
        event_peak_times_ms = [event_peak_times_ms, o.event_peak_times_ms]; %#ok<AGROW>
        file_of_roi     = [file_of_roi,     repmat(i, 1, fr.numUnits)]; %#ok<AGROW>

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
    gd.total_events      = total_events;
    gd.peak_frame        = peak_frame;
    gd.upper_threshold   = upper_threshold;
    gd.lower_threshold   = lower_threshold;
    gd.qc_category       = qc_category;
    gd.crossed_upper     = crossed_upper;
    gd.event_peak_times_ms = event_peak_times_ms;
    gd.file_of_roi       = file_of_roi;
    gd.filenames         = filenames;
    gd.coverslip_of_roi  = coverslip_of_roi;
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
    % and avoids per-ROI 0/0 ratios. Coverslips with no classified events are
    % dropped (NaN).
    cs   = gd.coverslip_of_roi;
    uniq = unique(cs, 'stable');
    fr   = nan(1, numel(uniq));
    for i = 1:numel(uniq)
        m   = strcmp(cs, uniq{i});
        s   = sum(gd.sync_count(m));
        a   = sum(gd.async_count(m));
        tot = s + a;
        if tot > 0, fr(i) = a / tot; end
    end
    fr = fr(isfinite(fr));
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
    xlim(tax.from_ms([-200, config.oneAP.async_window_ms(2) + 100]));   % default view
    ylabel('\DeltaF/F (mean \pm SEM)');
    title('1AP Average Evoked Trace (double-exponential decay fit)');
    legend(legend_entries, 'Location', 'best');
    grid on; box on;
    % Label in the chosen unit; override xlim with config.x_axis.range_frames if set.
    format_time_axis(gca, tax, 'Time from stimulus', stim_frame);
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

function fig = plot_oneap_metric(groups, field, useResponders, ylab, ttl, comparison, pfield)
    fig = figure('Name', ['OneAP ' ttl], 'Position', [150 120 700 600]);
    plot_group_box(gca, groups, field, useResponders, ylab, ttl, comparison, pfield);
end

function fig = plot_oneap_kinetics(groups, comparison)
    fig = figure('Name', 'OneAP Kinetics', 'Position', [200 120 1300 500]);
    subplot(1,3,1);
    plot_group_box(gca, groups, 'time_to_peak_ms', true, 'Time to peak (ms)', ...
        'Time to Peak', comparison, 'time_to_peak_ms');
    subplot(1,3,2);
    plot_group_box(gca, groups, 'auc', true, 'AUC (\DeltaF/F\cdotms)', ...
        'Area Under Curve', comparison, 'auc');
    subplot(1,3,3);
    plot_group_box(gca, groups, 'tau_ms', true, '\tau decay (ms)', ...
        'Decay \tau', comparison, 'tau_ms');
    sgtitle('1AP Response Kinetics (responders)', 'FontSize', 14, 'FontWeight', 'bold');
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
        labels{g} = sprintf('%s (N=%d cs, %d ROIs)', gd.name, numel(cm), numel(v));
    end
    set(ax, 'XTick', 1:nG, 'XTickLabel', labels);
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
    ylim([min(all_vals) - 0.05*(max(all_vals)-min(all_vals)+eps), y + 0.1*(max(all_vals)-min(all_vals)+eps)]);
end

function fig = plot_oneap_auc_split(groups, comparison)
    % AUC split by release timing: total / synchronous / asynchronous, all
    % responder-gated and aggregated per coverslip (same units as the total AUC).
    % The asynchronous panel is the magnitude of the slow tail seen in the
    % averaged trace - the quantity the per-ROI mono tau cannot capture.
    fig = figure('Name', 'OneAP AUC Split', 'Position', [180 110 1300 500]);
    subplot(1,3,1);
    plot_group_box(gca, groups, 'auc', true, 'AUC (\DeltaF/F\cdotms)', ...
        'Total AUC', comparison, 'auc');
    subplot(1,3,2);
    plot_group_box(gca, groups, 'auc_sync', true, 'AUC (\DeltaF/F\cdotms)', ...
        'Synchronous AUC', comparison, 'auc_sync');
    subplot(1,3,3);
    plot_group_box(gca, groups, 'auc_async', true, 'AUC (\DeltaF/F\cdotms)', ...
        'Asynchronous AUC (slow tail)', comparison, 'auc_async');
    sgtitle('1AP Area Under Curve by release timing (responders)', ...
        'FontSize', 14, 'FontWeight', 'bold');
end

function fig = plot_oneap_sync_async(groups, config)
    fig = figure('Name', 'OneAP Sync vs Async', 'Position', [250 100 1200 550]);
    tax = make_time_axis(config);

    % Left: sync/async fraction of events per group
    subplot(1,2,1); hold on;
    nG = numel(groups);
    fracs = zeros(nG, 2);
    labels = cell(1, nG);
    for g = 1:nG
        gd = groups(g);
        s = sum(gd.sync_count);
        a = sum(gd.async_count);
        tot = s + a;
        if tot == 0, tot = 1; end
        fracs(g, :) = [s, a] / tot;
        labels{g} = gd.name;
    end
    hb = bar(fracs, 'stacked');
    hb(1).FaceColor = [0.2 0.6 0.2]; hb(2).FaceColor = [0.8 0.5 0.1];
    set(gca, 'XTick', 1:nG, 'XTickLabel', labels);
    ylabel('Fraction of events'); ylim([0 1]);
    legend({'Synchronous', 'Asynchronous'}, 'Location', 'bestoutside');
    title('Synchronous vs Asynchronous event share'); grid on; box on;

    % Right: post-stim temporal histogram (PSTH)
    subplot(1,2,2); hold on;
    for g = 1:nG
        gd = groups(g);
        if isempty(gd.temporal_hist_sum), continue; end
        plot(tax.from_ms(gd.bin_centers_ms), gd.temporal_hist_sum, '-o', 'Color', gd.color, ...
            'LineWidth', 1.8, 'MarkerFaceColor', gd.color, 'MarkerSize', 4, ...
            'DisplayName', gd.name);
    end
    yl = ylim;
    sync_hi = tax.from_ms(config.oneAP.sync_window_ms(2));
    patch([0 sync_hi sync_hi 0], [0 0 yl(2) yl(2)], [0.2 0.6 0.2], ...
        'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    ylabel(sprintf('Event count (%d ms bins)', config.oneAP.bin_size_ms));
    title('Post-stimulus event time histogram');
    legend('Location', 'best'); grid on; box on;
    format_time_axis(gca, tax, 'Time from stimulus');
    sgtitle('1AP Synchronous / Asynchronous Release', 'FontSize', 14, 'FontWeight', 'bold');
end

function fig = plot_oneap_responders(groups, comparison)
    fig = figure('Name', 'OneAP Responder Fraction', 'Position', [300 150 700 550]);
    hold on;
    nG = numel(groups);
    fr = zeros(1, nG); labels = cell(1, nG);
    for g = 1:nG
        gd = groups(g);
        fr(g) = 100 * gd.n_responders / max(gd.total_rois, 1);
        labels{g} = sprintf('%s\n(%d/%d)', gd.name, gd.n_responders, gd.total_rois);
    end
    hb = bar(fr, 0.6, 'FaceColor', 'flat');
    for g = 1:nG, hb.CData(g, :) = groups(g).color; end

    % Overlay per-coverslip responder fractions (the statistical units)
    for g = 1:nG
        cf = 100 * coverslip_resp_frac(groups(g));
        if isempty(cf), continue; end
        cjit = (rand(1, numel(cf)) - 0.5) * 0.25;
        scatter(g + cjit, cf, 45, groups(g).color, 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1);
        labels{g} = sprintf('%s\n(%d/%d trials, N=%d cs)', groups(g).name, ...
            groups(g).n_responders, groups(g).total_rois, numel(cf));
    end
    set(gca, 'XTick', 1:nG, 'XTickLabel', labels);
    ylabel('Successful trials (%)');
    title('Trial success rate (ROI x trial units with a validated evoked response)');
    grid on; box on; ylim([0 100]);

    if ~isempty(comparison) && isfield(comparison, 'responder_fraction_p') ...
            && ~isnan(comparison.responder_fraction_p) && nG >= 2
        p = comparison.responder_fraction_p;
        y = max(fr) + 6;
        line([1 2], [y y], 'Color', 'k');
        if p < 0.001, ptxt = 'p < 0.001 ***'; else, ptxt = sprintf('p = %.3f', p); end
        text(1.5, y + 1, ptxt, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
end

function fig = plot_oneap_release_prob(groups, comparison)
    % Per-ROI release statistics across trials, over ACTIVE ROIs (>= 1 success):
    %   (1) distribution of P(response) = successes / valid trials
    %   (2) P(response), (3) potency, (4) efficacy - per-coverslip SuperPlots
    fig = figure('Name', 'OneAP Release Probability', 'Position', [140 90 1600 800]);
    nG = numel(groups);

    subplot(2,2,1); hold on;
    labels = cell(1, nG);
    maxT = max([1, arrayfun(@(g) max([g.roi_n_valid, 1]), groups)]);
    for g = 1:nG
        gd = groups(g);
        v = gd.roi_p_response(gd.roi_active & isfinite(gd.roi_p_response));
        labels{g} = sprintf('%s (%d active / %d ROIs)', gd.name, numel(v), gd.n_rois);
        if isempty(v), continue; end
        histogram(v, 'BinEdges', ((0:maxT+1) - 0.5) / maxT, ...   % centred on k/maxT
            'Normalization', 'probability', 'DisplayStyle', 'stairs', ...
            'EdgeColor', gd.color, 'LineWidth', 2);
    end
    xlim([0 1.05]); xlabel('P(response) per ROI'); ylabel('Fraction of active ROIs');
    title('Per-ROI response probability across trials');
    legend(labels, 'Location', 'northwest', 'Interpreter', 'none'); grid on; box on;

    subplot(2,2,2);
    plot_roi_box(gca, groups, 'roi_p_response', 'P(response)', ...
        'Response probability (active ROIs)', comparison);
    ylim([0 1.1]);
    subplot(2,2,3);
    plot_roi_box(gca, groups, 'roi_potency', '\DeltaF/F', ...
        'Potency (mean amplitude of successes)', comparison);
    subplot(2,2,4);
    plot_roi_box(gca, groups, 'roi_efficacy', '\DeltaF/F', ...
        'Efficacy (mean amplitude, failures = 0)', comparison);
    sgtitle('1AP release reliability across trials', 'FontSize', 14, 'FontWeight', 'bold');
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
        labels{g} = sprintf('%s (N=%d cs, %d ROIs)', gd.name, numel(cm), numel(v));
    end
    set(ax, 'XTick', 1:nG, 'XTickLabel', labels, 'TickLabelInterpreter', 'none');
    xlim([0.5, nG + 0.5]);
    ylabel(ylab); title(ttl); grid on; box on;
    annotate_pvalue(comparison, field, cs_vals);
end

function fig = plot_oneap_async_fraction(groups, comparison, config)
    % Per-coverslip ASYNCHRONOUS-release fraction: async/(sync+async) events.
    % Bars = group-pooled fraction (matches the PSTH); dots = per-coverslip
    % fractions (the statistical units); p-value = Mann-Whitney U across
    % coverslips. This is the per-replicate test of the slow tail that is obvious
    % in the averaged evoked trace but is NOT captured by the per-ROI mono tau
    % (truncated at the lower threshold) or the peak-dominated AUC.
    fig = figure('Name', 'OneAP Asynchronous Fraction', 'Position', [320 160 700 550]);
    hold on;
    nG = numel(groups);
    fr = zeros(1, nG); labels = cell(1, nG); cf_all = [];
    for g = 1:nG
        gd = groups(g);
        s = sum(gd.sync_count); a = sum(gd.async_count);
        tot = s + a; if tot == 0, tot = 1; end
        fr(g) = 100 * a / tot;
        labels{g} = gd.name;
    end
    hb = bar(fr, 0.6, 'FaceColor', 'flat');
    for g = 1:nG, hb.CData(g, :) = groups(g).color; end

    % Overlay per-coverslip async fractions (the statistical units)
    for g = 1:nG
        cf = 100 * coverslip_async_frac(groups(g));
        cf_all = [cf_all, cf]; %#ok<AGROW>
        if isempty(cf), continue; end
        cjit = (rand(1, numel(cf)) - 0.5) * 0.25;
        scatter(g + cjit, cf, 45, groups(g).color, 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1);
        labels{g} = sprintf('%s\n(N=%d cs)', groups(g).name, numel(cf));
    end
    set(gca, 'XTick', 1:nG, 'XTickLabel', labels);
    ylabel('Asynchronous fraction (%)');
    title(sprintf('Asynchronous release fraction (events >%g ms / all post-stim events)', ...
        config.oneAP.sync_window_ms(2)));
    grid on; box on;
    top = max([fr, cf_all, 1]);
    ylim([0, min(100, top * 1.25)]);

    if ~isempty(comparison) && isfield(comparison, 'async_fraction_p') ...
            && ~isnan(comparison.async_fraction_p) && nG >= 2
        p = comparison.async_fraction_p;
        y = top * 1.12;
        line([1 2], [y y], 'Color', 'k');
        if p < 0.001, ptxt = 'p < 0.001 ***'; else, ptxt = sprintf('p = %.3f', p); end
        text(1.5, y + top*0.02, ptxt, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
end

function fig = plot_oneap_qc(groups, config)
    fig = figure('Name', 'OneAP QC', 'Position', [120 80 1300 850]);
    tax = make_time_axis(config);

    % (1) pre-stim noise SD distribution
    subplot(2,2,1); hold on;
    for g = 1:numel(groups)
        gd = groups(g);
        v = gd.prestim_sd(isfinite(gd.prestim_sd));
        if isempty(v), continue; end
        histogram(v, 30, 'FaceColor', gd.color, 'FaceAlpha', 0.5, ...
            'EdgeColor', 'none', 'DisplayName', gd.name);
    end
    xlabel('Pre-stim noise SD (\DeltaF/F)'); ylabel('ROI count');
    title('Pre-stimulus noise (threshold basis)'); legend('Location','best'); grid on; box on;

    % (2) peak-time histogram (verify peak ~ expected), relative to stim
    subplot(2,2,2); hold on;
    binw = tax.from_frames_rel(1);   % one-frame bins in display units
    for g = 1:numel(groups)
        gd = groups(g);
        v = gd.peak_frame(gd.is_responder);
        v = v(isfinite(v));
        if isempty(v), continue; end
        v = tax.from_frames_rel(v - config.oneAP.stim_frame);   % time from stim
        histogram(v, 'BinWidth', binw, 'FaceColor', gd.color, ...
            'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', gd.name);
    end
    yl = ylim;
    plot([0 0], yl, 'k--', 'HandleVisibility', 'off');
    text(0, yl(2), ' stim', 'VerticalAlignment', 'top', 'FontSize', 8);
    ylabel('Responder count');
    title('Evoked peak time (responders)'); legend('Location','best'); grid on; box on;
    format_time_axis(gca, tax, 'Peak time from stimulus');

    % (3) tau fit R^2 distribution
    subplot(2,2,3); hold on;
    for g = 1:numel(groups)
        gd = groups(g);
        v = gd.tau_r2(isfinite(gd.tau_r2));
        if isempty(v), continue; end
        histogram(v, 0:0.05:1, 'FaceColor', gd.color, 'FaceAlpha', 0.5, ...
            'EdgeColor', 'none', 'DisplayName', gd.name);
    end
    xlabel('Decay fit R^2'); ylabel('ROI count');
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
    maxROIs = 40;
    ylim_fixed = [-0.02, 0.06];   % shared y-axis for all panels (dF/F)

    % Show a generous pre-stim segment so the resting noise is visible against
    % the evoked response (the response itself is fast: ~0-30 ms post-stim).
    win_pre  = round(150 / frame_dur);
    win_post = round((config.oneAP.async_window_ms(2) + 80) / frame_dur);

    for g = 1:numel(groups)
        gd = groups(g);
        if isempty(gd.dfof_all), continue; end

        % Prioritise responders, then fill with non-responders
        order = [find(gd.is_responder), find(~gd.is_responder)];
        order = order(order <= size(gd.dfof_all, 2));
        order = order(1:min(maxROIs, numel(order)));
        if isempty(order), continue; end

        fig = figure('Name', sprintf('OneAP Montage %s', gd.name), ...
            'Position', [50 40 1850 1000], 'Color', 'w');
        figs = [figs, fig]; %#ok<AGROW>
        nR = numel(order);
        nc = 5; nr = ceil(nR / nc);
        % tiledlayout with tight spacing/padding gives each panel more
        % horizontal room than subplot's default margins.
        tl = tiledlayout(fig, nr, nc, 'TileSpacing', 'tight', 'Padding', 'tight');

        nF_full = size(gd.dfof_all, 1);
        [f0, f1] = montage_window(tax, stim_frame, win_pre, win_post, nF_full);
        fidx  = f0:f1;
        t     = tax.from_ms((fidx - stim_frame) * frame_dur);
        tfull = tax.from_ms(((1:nF_full) - stim_frame) * frame_dur);  % FULL-trace axis
        xwin  = [t(1) t(end)];                                        % initial x-window

        for k = 1:nR
            roi = order(k);
            nexttile(tl, k); hold on;
            % Plot the WHOLE trace; the view starts on the stim window and the
            % rest is a horizontal drag away (see scroll_time_axis).
            plot(tfull, gd.dfof_all(:, roi), 'k-', 'LineWidth', 1);
            if isfinite(gd.upper_threshold(roi))
                yline(gd.upper_threshold(roi), 'g:', 'LineWidth', 1);
            end
            ylim(ylim_fixed);
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
        end
        title(tl, sprintf(['1AP ROI montage: %s (responders first; drag to scroll x, whole trace plotted; ' ...
            'blue=sync peak, red=async peak, magenta=\\tau fit; x in %s)'], gd.name, tax.unit), ...
            'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'tex');
    end
end

function figs = plot_oneap_marginal(groups, config)
    % Montage of the MARGINAL synchronous responders - the weakest evoked calls,
    % the ones sitting just above the detection threshold. These are where a
    % responder/non-responder mis-call is most likely, so inspecting them by eye
    % is the check on whether the filter is keeping noise. Same 5x4 = 20-panel
    % format and styling as plot_oneap_montage (extended baseline, green upper /
    % cyan lower threshold, stim marker, peak markers, tau fit).
    figs = [];
    stim_frame = config.oneAP.stim_frame;
    frame_dur  = 1000 / config.frame_rate;
    tax        = make_time_axis(config);
    maxROIs    = 20;                 % 5 col x 4 rows, matching the montage
    ylim_fixed = [-0.02, 0.06];      % shared y-axis for all panels (dF/F)

    win_pre  = round(150 / frame_dur);   % generous pre-stim baseline (see noise)
    win_post = round((config.oneAP.async_window_ms(2) + 80) / frame_dur);

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

        fig = figure('Name', sprintf('OneAP Marginal responders %s', gd.name), ...
            'Position', [50 40 1850 1000], 'Color', 'w');
        figs = [figs, fig]; %#ok<AGROW>
        nR = numel(order);
        nc = 5; nr = ceil(nR / nc);
        tl = tiledlayout(fig, nr, nc, 'TileSpacing', 'tight', 'Padding', 'tight');

        nF_full = size(gd.dfof_all, 1);
        [f0, f1] = montage_window(tax, stim_frame, win_pre, win_post, nF_full);
        fidx  = f0:f1;
        t     = tax.from_ms((fidx - stim_frame) * frame_dur);
        tfull = tax.from_ms(((1:nF_full) - stim_frame) * frame_dur);  % FULL-trace axis
        xwin  = [t(1) t(end)];                                        % initial x-window

        for k = 1:nR
            roi = order(k);
            nexttile(tl, k); hold on;
            % Whole trace plotted; view starts on the stim window, drag to scroll.
            plot(tfull, gd.dfof_all(:, roi), 'k-', 'LineWidth', 1);
            if isfinite(gd.upper_threshold(roi))
                yline(gd.upper_threshold(roi), 'g:', 'LineWidth', 1);
            end
            if isfinite(gd.lower_threshold(roi))
                yline(gd.lower_threshold(roi), 'c:', 'LineWidth', 1);
            end
            ylim(ylim_fixed);
            xline(0, '--', 'Color', 'k', 'LineWidth', 1);
            plot_peak_markers(gd, roi, stim_frame, frame_dur, 1, nF_full, ...
                config.oneAP.sync_window_ms, config.oneAP.async_window_ms, tax);
            plot_tau_fit(gd, roi, stim_frame, frame_dur, nF_full, tax);
            snr_roi = gd.peak_amplitude(roi) / gd.prestim_sd(roi);
            title(sprintf('%s SNR=%.1f', unit_label(gd, roi), snr_roi), ...
                'FontSize', 8, 'Interpreter', 'tex');
            set(gca, 'FontSize', 7); grid on; box on;
            scroll_time_axis(gca, xwin);   % show stim window; drag to scroll x
        end
        title(tl, sprintf(['1AP marginal responders: %s (weakest evoked calls, ' ...
            'low SNR first; drag to scroll x; green=%g\\sigma, cyan=%g\\sigma, ' ...
            'blue=sync peak, red=async peak, magenta=\\tau fit; x in %s)'], gd.name, ...
            config.oneAP.upper_threshold_sigma, config.oneAP.lower_threshold_sigma, tax.unit), ...
            'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'tex');
    end
end

function s = unit_label(gd, u)
    % "F<file> ROI<n> T<trial>" for a (ROI, trial) unit column of the group.
    s = sprintf('F%d ROI%d T%d', gd.file_of_roi(u), gd.unit_roi(u), gd.unit_trial(u));
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
    % Overlay every counted peak on a per-ROI trace:
    %   blue filled circle  = synchronous peak (peak time in sync_win)
    %   red  filled circle  = asynchronous peak (peak time in async_win)
    % Peak times (ms from stim) come straight from the analyzer; the y value is
    % read back from the aligned dF/F trace at the corresponding frame.
    if roi > numel(gd.event_peak_times_ms), return; end
    times = gd.event_peak_times_ms{roi};
    for i = 1:numel(times)
        tt = times(i);
        fr = stim_frame + round(tt / frame_dur);
        if fr < f0 || fr > f1, continue; end
        yv = gd.dfof_all(fr, roi);
        if ~isfinite(yv), continue; end
        xx = tax.from_ms(tt);   % display units from stim
        if tt >= sync_win(1) && tt <= sync_win(2)
            plot(xx, yv, 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
        elseif tt > async_win(1) && tt <= async_win(2)
            plot(xx, yv, 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
        end
    end
end

%% ======================================================================
%% Output / utility
%% ======================================================================
function T = build_metrics_table(groups)
    % One row per (ROI, trial) unit. step1_flags = that trial's Step1_Flags.
    genotype = {}; drug = {}; filename = {}; coverslip = {}; roi = []; trial = [];
    flags = {}; responder = [];
    peak_amp = []; ttp = []; auc = []; auc_s = []; auc_a = []; tau = []; tau_r2 = [];
    sync_c = []; async_c = []; tot_ev = []; pre_sd = [];

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
        tot_ev    = [tot_ev;    gd.total_events(:)];       %#ok<AGROW>
        pre_sd    = [pre_sd;    gd.prestim_sd(:)];         %#ok<AGROW>
    end

    T = table(genotype, drug, filename, coverslip, roi, trial, logical(responder), peak_amp, ttp, auc, ...
        auc_s, auc_a, tau, tau_r2, sync_c, async_c, tot_ev, pre_sd, flags, ...
        'VariableNames', {'genotype', 'drug', 'filename', 'coverslip', 'roi', 'trial', 'is_responder', ...
        'peak_amplitude', 'time_to_peak_ms', 'auc', 'auc_sync', 'auc_async', 'tau_ms', 'tau_r2', ...
        'sync_count', 'async_count', 'total_events', 'prestim_sd', 'step1_flags'});
end

function T = build_roi_table(groups)
    % One row per physical ROI: release statistics across its trials.
    parts = cell(numel(groups), 1);
    for g = 1:numel(groups)
        gd = groups(g);
        n = gd.n_rois;
        parts{g} = table(repmat({gd.genotype}, n, 1), repmat({gd.drug}, n, 1), ...
            gd.filenames(gd.roi_file)', gd.roi_coverslip(:), gd.roi_index(:), ...
            gd.roi_n_valid(:), gd.roi_n_success(:), gd.roi_p_response(:), ...
            gd.roi_potency(:), gd.roi_efficacy(:), ...
            'VariableNames', {'genotype', 'drug', 'filename', 'coverslip', 'roi', ...
            'n_valid_trials', 'n_success', 'p_response', 'potency', 'efficacy'});
    end
    T = vertcat(parts{:});
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
