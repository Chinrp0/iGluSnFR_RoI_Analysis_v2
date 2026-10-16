function results = spont_batch_analysis(dataset, options)
    % SPONT_BATCH_ANALYSIS - Batch analysis of spontaneous glutamate release
    %
    % Counterpart to oneap_batch_analysis / twoap_batch_analysis for spontaneous
    % recordings, taking a dataset from load_dataset_folder (one frame interval;
    % run_pipeline splits mixed folders first). Frame interval and frame count
    % come from Recordings.csv per file. Every file is processed once by
    % spont_pipeline (rolling-median
    % baseline -> dF/F -> corrected Schmitt event detection). Files are grouped
    % by genotype x drug: genotype from config.groups (parse_group, fails loudly
    % on an unrecognised or ambiguous name) and drug from the "+drug" token
    % (parse_drug). The plots compare TWO groups at a time, so the comparison
    % figures are made once per genotype pair within each drug condition.
    %
    % Inputs:
    %   dataset - struct from load_dataset_folder (spontaneous recordings)
    %   options - struct: verbose, createPlots, create_*_plot switches
    %             (frame_rate / recording_duration_s are set from the metadata)
    %
    % Output:
    %   results - .groups (one aggregate per genotype x drug group),
    %             .comparisons (one per pair: .wt_data = group A,
    %             .mut_data = group B, .comparison, .plot_handles, ...),
    %             .processing_info

    if nargin < 2
        options = struct();
    end

    % Set defaults
    if ~isfield(options, 'verbose'), options.verbose = true; end
    if ~isfield(options, 'createPlots'), options.createPlots = true; end
    if ~isfield(options, 'useParallel'), options.useParallel = true; end
    if ~isfield(options, 'continueOnError'), options.continueOnError = true; end
    if ~isfield(options, 'savePlots'), options.savePlots = false; end

    if ~isstruct(dataset) || ~isfield(dataset, 'files')
        error('spont_batch_analysis:input', ...
            'spont_batch_analysis takes a dataset from load_dataset_folder.');
    end
    folder_path = dataset.folder;
    evalc('config = tracenorm_config();');
    config.dataType = 'Spont';

    fprintf('=== Spontaneous Batch Analysis (genotype x drug) ===\n');
    fprintf('Folder: %s\n', folder_path);

    %% Step 1: files, groups and the shared time base (from Recordings.csv)
    entries = dataset.files;
    if isempty(entries)
        error('spont_batch_analysis:noFiles', 'The dataset has no loadable files.');
    end
    fi = arrayfun(@(f) f.meta.frame_interval_ms, entries);
    if (max(fi) - min(fi)) > config.metadata.interval_tol_frac * median(fi)
        error('spont_batch_analysis:mixedInterval', ...
            'Files have different Frame_Interval_ms; split them first (run_pipeline does).');
    end
    nF = arrayfun(@(f) f.meta.n_frames, entries);
    config = apply_recording_timing(config, median(fi), max(nF));
    % The plot modules take the time base from options (rates, raster length).
    options.frame_rate           = config.frame_rate;
    options.recording_duration_s = config.recording_duration_s;
    options.max_iei_s            = config.visualization.spont.iei_max_s;   % IEI plot range
    options.iei_bin_size_ms      = config.visualization.spont.iei_bin_ms;  % IEI histogram bin
    options.viz                  = config.visualization.spont;    % bins, subsampling, grids
    options.ymin_trace           = config.visualization.ymin_montage;  % example-trace y floor
    fprintf('  Frame interval %.4g ms (%.4g Hz), %d-%d frames (metadata)\n', ...
        config.frame_dur_ms, config.frame_rate, min(nF), max(nF));

    files     = {entries.path};
    genotypes = {entries.genotype};
    drugs     = {entries.drug};
    if options.verbose
        fprintf('\nStep 1: Files by genotype x drug...\n');
        for i = 1:numel(entries)
            fprintf('  %-12s %-10s %s\n', genotypes{i}, drugs{i}, entries(i).name);
        end
    end

    %% Step 2: process every file once
    if options.verbose
        fprintf('\nStep 2: Processing %d files...\n', numel(files));
    end
    tic;
    file_results = process_files(entries, config);
    process_time = toc;
    ok = ~cellfun(@isempty, file_results);
    fprintf('  Processed %d/%d files successfully in %.3f s\n', sum(ok), numel(files), process_time);

    %% Step 3: aggregate per genotype x drug group
    if options.verbose
        fprintf('\nStep 3: Aggregating by genotype x drug...\n');
    end
    [geno_list, drug_list] = genotype_drug_order(genotypes(ok), drugs(ok), config);
    groups = struct('name', {}, 'genotype', {}, 'drug', {}, 'data', {}, ...
        'file_results', {}, 'files', {});
    for d = 1:numel(drug_list)
        for g = 1:numel(geno_list)
            mask = ok & strcmp(genotypes, geno_list{g}) & strcmp(drugs, drug_list{d});
            if ~any(mask), continue; end
            name = sprintf('%s %s', geno_list{g}, drug_list{d});
            gd = aggregate_condition_data_integrated(file_results(mask), files(mask), name);
            gd.color = group_color(geno_list{g}, d, numel(drug_list), config);
            groups(end+1) = struct('name', name, 'genotype', geno_list{g}, ...
                'drug', drug_list{d}, 'data', gd, 'file_results', {file_results(mask)}, ...
                'files', {files(mask)}); %#ok<AGROW>
            fprintf('  %s: %d files, %d total ROIs, %d total events\n', ...
                name, gd.num_files, gd.total_rois, gd.total_events);
        end
    end

    %% Steps 4-6: statistics + plots for each genotype pair within a drug
    comparisons = {};
    for d = 1:numel(drug_list)
        idx = find(strcmp({groups.drug}, drug_list{d}));
        for a = 1:numel(idx)
            for b = a+1:numel(idx)
                comparisons{end+1} = compare_pair(groups(idx(a)), groups(idx(b)), ...
                    folder_path, process_time, options); %#ok<AGROW>
            end
        end
    end
    comparisons = [comparisons{:}];
    if isempty(comparisons)
        warning('spont_batch_analysis:noPairs', ...
            'No drug condition has two genotypes; groups aggregated but nothing compared.');
    end

    results = struct();
    results.groups          = groups;
    results.comparisons     = comparisons;
    results.processing_info = struct('folder_path', folder_path, 'files', {files}, ...
        'genotypes', {genotypes}, 'drugs', {drugs}, 'process_time', process_time, ...
        'timestamp', datetime('now'), ...
        'metadata_file', dataset.metadata_file, 'frame_interval_ms', config.frame_dur_ms, ...
        'frame_rate', config.frame_rate, 'recording_duration_s', config.recording_duration_s);
end

function results = compare_pair(ga, gb, folder_path, process_time, options)
    % Two-group comparison (A = ga, B = gb) with the existing spont plots. The
    % result is shaped like the old WT-vs-R213W results (wt_data = A slot,
    % mut_data = B slot) because every plot module reads those fields; the
    % real group names/colours travel in .condition / .color.
    fprintf('\n--- %s vs %s ---\n', ga.name, gb.name);
    wt_data  = ga.data;
    mut_data = gb.data;

    try
        %% Step 4: Statistical comparison
        if options.verbose
            fprintf('Step 4: Statistical comparison...\n');
        end
        comparison = compare_conditions_statistical_integrated(wt_data, mut_data);

        %% Step 5: Create core comparison visualizations
        plot_handles = [];
        if options.createPlots
            if options.verbose
                fprintf('Step 5: Creating core comparison plots...\n');
            end
            plot_handles = create_plots_integrated_final(wt_data, mut_data, comparison, options.viz);
            if options.savePlots
                save_condition_plots_integrated(plot_handles, folder_path);
            end
        end

        %% Step 5.5: assemble the pair's results (the plot modules read these)
        results = struct();
        results.group_a = ga.name;
        results.group_b = gb.name;
        results.drug    = ga.drug;
        results.wt_data = wt_data;
        results.mut_data = mut_data;
        results.comparison = comparison;
        results.file_results = struct('wt_results', {ga.file_results}, 'mut_results', {gb.file_results});
        results.plot_handles = plot_handles;
        results.processing_info = struct();
        results.processing_info.folder_path = folder_path;
        results.processing_info.wt_files = ga.files;
        results.processing_info.mut_files = gb.files;
        results.processing_info.process_time = process_time;
        results.processing_info.timestamp = datetime('now');
        results.processing_info.frame_rate = options.frame_rate;
        results.processing_info.recording_duration_s = options.recording_duration_s;

        %% Step 6: Create ADDITIONAL modular visualizations
        additional_plots = struct();
        
        if options.createPlots && options.verbose
            fprintf('\nStep 6: Creating additional analysis plots...\n');
        end
        
        % Inter-Event Interval analysis
        if options.createPlots && getfield_safe(options, 'create_iei_plot', true)
            try
                plot_options = struct('frame_rate', options.frame_rate, ...
                    'max_iei_s', options.max_iei_s, ...
                    'bin_size_ms', options.iei_bin_size_ms);
                additional_plots.iei_comparison = plot_iei_comparison(results, plot_options);
            catch ME
                fprintf('  WARNING: Failed to create IEI plot: %s\n', ME.message);
            end
        end
        
        % Amplitude-Frequency correlation
        if options.createPlots && getfield_safe(options, 'create_amp_freq_correlation', true)
            try
                additional_plots.amp_freq_correlation = plot_amplitude_frequency_correlation(results);
            catch ME
                fprintf('  WARNING: Failed to create amp-freq correlation plot: %s\n', ME.message);
            end
        end
        
        % Event timing raster
        if options.createPlots && getfield_safe(options, 'create_event_raster', true)
            try
                plot_options = struct('frame_rate', options.frame_rate, ...
                    'recording_duration_s', options.recording_duration_s, ...
                    'max_rois_per_plot', options.viz.raster_max_rois, ...
                    'rate_bin_s', options.viz.rate_bin_s);
                additional_plots.event_raster = plot_event_raster(results, plot_options);
            catch ME
                fprintf('  WARNING: Failed to create event raster plot: %s\n', ME.message);
            end
        end
        
        % Biological variability
        if options.createPlots && getfield_safe(options, 'create_biological_variability', true)
            try
                additional_plots.biological_variability = plot_biological_variability(results);
            catch ME
                fprintf('  WARNING: Failed to create biological variability plot: %s\n', ME.message);
            end
        end
        
        % ROI recruitment curves
        if options.createPlots && getfield_safe(options, 'create_recruitment_curves', true)
            try
                additional_plots.recruitment_curves = plot_roi_recruitment_curves(results, ...
                    struct('thresholds_hz', options.viz.recruitment_thresholds_hz, ...
                           'key_thresholds_hz', options.viz.recruitment_key_thresholds_hz));
            catch ME
                fprintf('  WARNING: Failed to create recruitment curves: %s\n', ME.message);
            end
        end
        
        % Cumulative event count
        if options.createPlots && getfield_safe(options, 'create_cumulative_events', true)
            try
                plot_options = struct('frame_rate', options.frame_rate, ...
                    'recording_duration_s', options.recording_duration_s, ...
                    'smooth_points', options.viz.cumulative_smooth_points);
                additional_plots.cumulative_events = plot_cumulative_events_over_time(results, plot_options);
            catch ME
                fprintf('  WARNING: Failed to create cumulative events plot: %s\n', ME.message);
            end
        end
        
        % ROI trace comparison
        if options.createPlots && getfield_safe(options, 'create_roi_traces', true)
            try
                plot_options = struct('frame_rate', options.frame_rate, 'ymin', options.ymin_trace);
                additional_plots.roi_traces = plot_condition_roi_traces(results, plot_options);
            catch ME
                fprintf('  WARNING: Failed to create ROI traces: %s\n', ME.message);
            end
        end
        
        % Store additional plots in results
        results.additional_plots = additional_plots;
        
        % Print summary
        if options.verbose
            print_condition_summary_integrated(results);
        end
        
    catch ME
        fprintf('ERROR in integrated batch analysis: %s\n', ME.message);
        if ~isempty(ME.stack)
            fprintf('  In: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
        end
        rethrow(ME);
    end
end

%% === INTEGRATED PLOTTING FUNCTIONS ===

function plot_handles = create_plots_integrated_final(wt_data, mut_data, comparison, viz)
    % Core condition-comparison plots (active-ROI focus); bins etc. from viz
    % (config.visualization.spont)
    
    plot_handles = struct();
    fprintf('  Creating improved condition comparison plots...\n');
    
    % Frequency comparison with mean labels on box plot
    plot_handles.frequency_comparison = create_frequency_plot_final(wt_data, mut_data, comparison, viz);
    
    % Amplitude comparison
    plot_handles.amplitude_comparison = create_amplitude_plot_final(wt_data, mut_data, comparison, viz);
    
    % Active ROI comparison (no 0-event ROIs)
    plot_handles.activity_comparison = create_activity_plot_final(wt_data, mut_data, comparison, viz);
    
    % File summary
    plot_handles.file_summary = create_file_summary_final(wt_data, mut_data, comparison);
    
    fprintf('    Created 4 improved condition comparison plots\n');
end

function fig = create_frequency_plot_final(wt_data, mut_data, comparison, viz)
    [lblA, lblB, colA, colB] = spont_group_labels(wt_data, mut_data);   % group labels/colours
    % FREQUENCY PLOT: viz.freq_bin_hz bins, mean labels on box plot
    
    fig = figure('Name', ['Event Frequency Comparison: ' lblA ' vs ' lblB ' (Active ROIs Only)'], ...
        'Position', [100, 100, 1200, 800]);
    
    wt_freqs = wt_data.all_roi_frequencies(wt_data.all_roi_frequencies > 0);
    mut_freqs = mut_data.all_roi_frequencies(mut_data.all_roi_frequencies > 0);
    
    % Subplot 1: histogram with viz.freq_bin_hz bins
    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile(tl, 1);
    max_freq = max([wt_freqs, mut_freqs]);
    edges = 0:viz.freq_bin_hz:(max_freq + viz.freq_bin_hz);

    hold on;
    histogram(wt_freqs, edges, 'FaceColor', colA, 'FaceAlpha', 0.7, ...
        'EdgeColor', 'none', 'DisplayName', sprintf([lblA ' (n=%d active ROIs)'], length(wt_freqs)));
    histogram(mut_freqs, edges, 'FaceColor', colB, 'FaceAlpha', 0.7, ...
        'EdgeColor', 'none', 'DisplayName', sprintf([lblB ' (n=%d active ROIs)'], length(mut_freqs)));
    
    xlabel('Event Frequency (Hz)');
    ylabel('Active ROI Count');
    title('Event Frequency Distribution (Active ROIs Only)');
    legend('Location', 'best');
    grid on;
    
    % Store means for subplot 2
    wt_mean = mean(wt_freqs);
    mut_mean = mean(mut_freqs);
    
    % Subplot 2: Box plot with mean labels
    nexttile(tl, 2);
    all_freqs = [wt_freqs, mut_freqs];
    groups = [ones(1, length(wt_freqs)), 2*ones(1, length(mut_freqs))];
    
    boxplot(all_freqs, groups, 'Labels', {lblA, lblB}, 'Colors', [colA; colB], 'Symbol', '+');
    ylabel('Event Frequency (Hz)');
    title('Frequency Box Plot Comparison');
    grid on;
    
    % Add mean frequency labels
    hold on;
    y_range = max(all_freqs) - min(all_freqs);
    label_box_means([wt_mean, mut_mean], {lblA, lblB}, [colA; colB], y_range, '%.3f Hz');
    
    % Statistical annotation
    if isfield(comparison, 'frequency_pvalue') && ~isnan(comparison.frequency_pvalue)
        y_max = max(all_freqs);
        y_text = y_max + y_range * 0.20;
        
        if comparison.frequency_pvalue < 0.001
            p_text = 'p < 0.001 ***';
        else
            p_text = sprintf('p = %.4f', comparison.frequency_pvalue);
        end
        text(1.5, y_text, p_text, 'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
        
        % Significance bracket
        line([1, 2], [y_text*0.95, y_text*0.95], 'Color', 'k', 'LineWidth', 1);
        line([1, 1], [y_text*0.95, y_text*0.9], 'Color', 'k', 'LineWidth', 1);
        line([2, 2], [y_text*0.95, y_text*0.9], 'Color', 'k', 'LineWidth', 1);
        
        ylim([0, y_text + y_range*0.05]);
    end
    
    % Subplot 3: Cumulative distributions
    nexttile(tl, 3);
    [wt_f, wt_x] = ecdf(wt_freqs);
    [mut_f, mut_x] = ecdf(mut_freqs);
    
    plot(wt_x, wt_f, '-', 'Color', colA, 'LineWidth', 2, 'DisplayName', lblA);
    hold on;
    plot(mut_x, mut_f, '-', 'Color', colB, 'LineWidth', 2, 'DisplayName', lblB);
    
    xlabel('Event Frequency (Hz)');
    ylabel('Cumulative Probability');
    title('Cumulative Distribution Functions');
    legend('Location', 'best');
    grid on;
    
    % Subplot 4: Summary
    nexttile(tl, 4); axis off;
    summary_text = {
        'FREQUENCY COMPARISON SUMMARY';
        '';
        sprintf([lblA ' (n = %d active ROIs):'], length(wt_freqs));
        sprintf('  Mean: %.4f ± %.4f Hz', wt_mean, std(wt_freqs));
        sprintf('  Median: %.4f Hz', median(wt_freqs));
        '';
        sprintf([lblB ' (n = %d active ROIs):'], length(mut_freqs));
        sprintf('  Mean: %.4f ± %.4f Hz', mut_mean, std(mut_freqs));
        sprintf('  Median: %.4f Hz', median(mut_freqs));
        '';
        'STATISTICAL TEST:';
        'Mann-Whitney U test';
        p_label(comparison, 'frequency_pvalue');
    };
    
    text(0.05, 0.95, summary_text, 'FontSize', 10, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', 'Units', 'normalized', 'FontName', 'FixedWidth');
    
    title(tl, ['Event Frequency Comparison: ' lblA ' vs ' lblB ' (Active ROIs Only)'], 'FontSize', 16, 'FontWeight', 'bold');
end

function fig = create_activity_plot_final(wt_data, mut_data, comparison, viz)
    [lblA, lblB, colA, colB] = spont_group_labels(wt_data, mut_data);   % group labels/colours
    % ACTIVITY PLOT: Active ROIs only (no 0-event ROIs)
    
    fig = figure('Name', ['Active ROI Comparison: ' lblA ' vs ' lblB], 'Position', [300, 100, 1200, 800]);
    
    wt_events = wt_data.all_roi_event_counts(wt_data.all_roi_event_counts > 0);
    mut_events = mut_data.all_roi_event_counts(mut_data.all_roi_event_counts > 0);
    
    % Events per active ROI
    tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile(tl, 1);
    if ~isempty(wt_events) && ~isempty(mut_events)
        max_events = max([wt_events, mut_events]);
        edges = 0.5:1:(max_events + 0.5);
        
        hold on;
        histogram(wt_events, edges, 'FaceColor', colA, 'FaceAlpha', 0.7, ...
            'EdgeColor', 'none', 'DisplayName', sprintf([lblA ' (n=%d active ROIs)'], length(wt_events)));
        histogram(mut_events, edges, 'FaceColor', colB, 'FaceAlpha', 0.7, ...
            'EdgeColor', 'none', 'DisplayName', sprintf([lblB ' (n=%d active ROIs)'], length(mut_events)));
        
        xlabel('Events per ROI'); ylabel('Active ROI Count');
        title('Events per Active ROI Distribution'); legend('Location', 'best'); grid on;
    end
    
    % Activity categories (active ROIs only)
    nexttile(tl, 2);
    if ~isempty(wt_events) && ~isempty(mut_events)
        activity_bins = activity_categories(viz.activity_bin_edges);

        wt_counts = zeros(1, size(activity_bins, 1));
        mut_counts = zeros(1, size(activity_bins, 1));
        
        for i = 1:size(activity_bins, 1)
            wt_counts(i) = sum(activity_bins{i, 2}(wt_events));
            mut_counts(i) = sum(activity_bins{i, 2}(mut_events));
        end
        
        hb = bar([wt_counts; mut_counts]', 'grouped');
        hb(1).FaceColor = colA; hb(2).FaceColor = colB;
        xlabel('Activity Category'); ylabel('Active ROI Count');
        title('Active ROI Categories (No Inactive ROIs)');
        set(gca, 'XTickLabel', activity_bins(:, 1));
        legend(lblA, lblB, 'Location', 'best'); grid on; xtickangle(45);
    end
    
    title(tl, ['Active ROI Comparison: ' lblA ' vs ' lblB ' (Inactive ROIs Excluded)'], 'FontSize', 16, 'FontWeight', 'bold');
end

%% === OTHER REQUIRED FUNCTIONS ===

function fig = create_amplitude_plot_final(wt_data, mut_data, comparison, viz)
    [lblA, lblB, colA, colB] = spont_group_labels(wt_data, mut_data);   % group labels/colours
    % AMPLITUDE PLOT: Full comparison with histograms and box plots
    
    fig = figure('Name', ['Event Amplitude Comparison: ' lblA ' vs ' lblB], 'Position', [200, 100, 1200, 800]);
    
    wt_amps = wt_data.all_individual_events;
    mut_amps = mut_data.all_individual_events;
    
    if isempty(wt_amps) || isempty(mut_amps)
        text(0.5, 0.5, 'Insufficient amplitude data for comparison', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 14);
        title('Event Amplitude Comparison - INSUFFICIENT DATA');
        return;
    end
    
    % 2x2 layout
    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile(tl, 1);
    min_amp = min([wt_amps, mut_amps]);
    max_amp = max([wt_amps, mut_amps]);
    edges = linspace(min_amp, max_amp, 50);
    
    hold on;
    histogram(wt_amps, edges, 'FaceColor', colA, 'FaceAlpha', 0.7, ...
        'EdgeColor', 'none', 'DisplayName', sprintf([lblA ' (n=%d events)'], length(wt_amps)));
    histogram(mut_amps, edges, 'FaceColor', colB, 'FaceAlpha', 0.7, ...
        'EdgeColor', 'none', 'DisplayName', sprintf([lblB ' (n=%d events)'], length(mut_amps)));
    
    xlabel('Event Amplitude (dF/F)');
    ylabel('Event Count');
    title('Event Amplitude Distribution');
    legend('Location', 'best');
    grid on;
    
    % Add mean lines
    wt_mean = mean(wt_amps);
    mut_mean = mean(mut_amps);
    y_lim = ylim;
    plot([wt_mean, wt_mean], y_lim, '--', 'Color', colA, 'LineWidth', 2, 'HandleVisibility', 'off');
    plot([mut_mean, mut_mean], y_lim, '--', 'Color', colB, 'LineWidth', 2, 'HandleVisibility', 'off');
    
    nexttile(tl, 2);
    % Sample data for box plot if too large (seeded, so the figure is reproducible)
    rs = RandStream('mt19937ar', 'Seed', viz.random_seed);
    wt_sample = datasample(rs, wt_amps, min(viz.box_max_points, length(wt_amps)), 'Replace', false);
    mut_sample = datasample(rs, mut_amps, min(viz.box_max_points, length(mut_amps)), 'Replace', false);
    
    all_amps = [wt_sample, mut_sample];
    groups = [ones(1, length(wt_sample)), 2*ones(1, length(mut_sample))];
    boxplot(all_amps, groups, 'Labels', {lblA, lblB}, 'Colors', [colA; colB], 'Symbol', '+');
    ylabel('Event Amplitude (dF/F)');
    title(sprintf('Amplitude Box Plot Comparison (box: <=%d events/group; means: all events)', ...
        viz.box_max_points));
    grid on;

    % Mean amplitude labels (means of ALL events, not the box-plot sample)
    hold on;
    label_box_means([wt_mean, mut_mean], {lblA, lblB}, [colA; colB], ...
        max(all_amps) - min(all_amps), '%.4f dF/F');
    
    % Add statistical annotation
    if isfield(comparison, 'amplitude_pvalue') && ~isnan(comparison.amplitude_pvalue)
        y_max = max(all_amps);
        y_range = max(all_amps) - min(all_amps);
        y_text = y_max + y_range * 0.15;
        
        if comparison.amplitude_pvalue < 0.001
            p_text = 'p < 0.001 ***';
        else
            p_text = sprintf('p = %.4f', comparison.amplitude_pvalue);
        end
        text(1.5, y_text, p_text, 'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
        
        % Significance bracket
        line([1, 2], [y_text*0.9, y_text*0.9], 'Color', 'k', 'LineWidth', 1);
        line([1, 1], [y_text*0.9, y_text*0.85], 'Color', 'k', 'LineWidth', 1);
        line([2, 2], [y_text*0.9, y_text*0.85], 'Color', 'k', 'LineWidth', 1);
        
        ylim([0, y_text + y_range*0.05]);
    end
    
    nexttile(tl, 3);
    [wt_f, wt_x] = ecdf(wt_amps);
    [mut_f, mut_x] = ecdf(mut_amps);
    
    plot(wt_x, wt_f, '-', 'Color', colA, 'LineWidth', 2, 'DisplayName', lblA);
    hold on;
    plot(mut_x, mut_f, '-', 'Color', colB, 'LineWidth', 2, 'DisplayName', lblB);
    
    xlabel('Event Amplitude (dF/F)');
    ylabel('Cumulative Probability');
    title('Cumulative Distribution Functions');
    legend('Location', 'best');
    grid on;
    
    nexttile(tl, 4); axis off;
    
    summary_text = {
        'AMPLITUDE COMPARISON SUMMARY';
        '';
        sprintf([lblA ' (n = %d events):'], length(wt_amps));
        sprintf('  Mean: %.4f ± %.4f dF/F', wt_mean, std(wt_amps));
        sprintf('  Median: %.4f dF/F', median(wt_amps));
        sprintf('  Range: %.4f - %.4f dF/F', min(wt_amps), max(wt_amps));
        '';
        sprintf([lblB ' (n = %d events):'], length(mut_amps));
        sprintf('  Mean: %.4f ± %.4f dF/F', mut_mean, std(mut_amps));
        sprintf('  Median: %.4f dF/F', median(mut_amps));
        sprintf('  Range: %.4f - %.4f dF/F', min(mut_amps), max(mut_amps));
        '';
        'STATISTICAL TEST:';
        'Mann-Whitney U test';
    };
    
    if isfield(comparison, 'amplitude_pvalue') && ~isnan(comparison.amplitude_pvalue)
        if comparison.amplitude_pvalue < 0.001
            summary_text{end+1} = 'p < 0.001 ***';
        else
            summary_text{end+1} = sprintf('p = %.4f', comparison.amplitude_pvalue);
        end
    end
    
    text(0.05, 0.95, summary_text, 'FontSize', 10, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', 'Units', 'normalized', 'FontName', 'FixedWidth');
    
    title(tl, ['Event Amplitude Comparison: ' lblA ' vs ' lblB], 'FontSize', 16, 'FontWeight', 'bold');
end

function label_box_means(means, labels, colors, y_range, fmt)
    % Mean diamond + boxed "<group>: <mean>" label over each box (x = 1, 2, ...)
    % in the group's colour; labels staggered so neighbouring boxes don't overlap.
    for k = 1:numel(means)
        plot(k, means(k), 'd', 'MarkerSize', 8, 'MarkerEdgeColor', colors(k, :), ...
            'MarkerFaceColor', colors(k, :), 'LineWidth', 2);
        text(k, means(k) + y_range * (0.02 + 0.03 * k), sprintf(['%s: ' fmt], labels{k}, means(k)), ...
            'HorizontalAlignment', 'center', 'BackgroundColor', 'white', ...
            'EdgeColor', colors(k, :), 'FontSize', 10, 'Interpreter', 'none');
    end
end

function fig = create_file_summary_final(wt_data, mut_data, comparison)
    [lblA, lblB, colA, colB] = spont_group_labels(wt_data, mut_data);   % group labels/colours
    % FILE SUMMARY: Per-file comparisons and overall summary
    
    fig = figure('Name', ['File-Level Summary: ' lblA ' vs ' lblB], 'Position', [400, 100, 1200, 800]);
    
    % 2x2 layout
    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile(tl, 1);
    if ~isempty(wt_data.file_summaries) && ~isempty(mut_data.file_summaries)
        wt_events_per_file = arrayfun(@(x) x.num_events, wt_data.file_summaries);
        mut_events_per_file = arrayfun(@(x) x.num_events, mut_data.file_summaries);
        
        hold on;
        bar(1:length(wt_events_per_file), wt_events_per_file, 0.4, ...
            'FaceColor', colA, 'DisplayName', lblA);
        bar((1:length(mut_events_per_file)) + 0.4, mut_events_per_file, 0.4, ...
            'FaceColor', colB, 'DisplayName', lblB);
        
        xlabel('File Index');
        ylabel('Total Events');
        title('Total Events per File');
        legend('Location', 'best');
        grid on;
    end
    
    nexttile(tl, 2);
    if ~isempty(wt_data.file_summaries) && ~isempty(mut_data.file_summaries)
        wt_active_per_file = arrayfun(@(x) x.active_rois, wt_data.file_summaries);
        mut_active_per_file = arrayfun(@(x) x.active_rois, mut_data.file_summaries);
        
        hold on;
        bar(1:length(wt_active_per_file), wt_active_per_file, 0.4, ...
            'FaceColor', colA, 'DisplayName', lblA);
        bar((1:length(mut_active_per_file)) + 0.4, mut_active_per_file, 0.4, ...
            'FaceColor', colB, 'DisplayName', lblB);
        
        xlabel('File Index');
        ylabel('Active ROIs');
        title('Active ROIs per File');
        legend('Location', 'best');
        grid on;
    end
    
    nexttile(tl, 3);
    if ~isempty(wt_data.file_summaries) && ~isempty(mut_data.file_summaries)
        wt_activity_fractions = arrayfun(@(x) 100 * x.active_rois / x.num_rois, wt_data.file_summaries);
        mut_activity_fractions = arrayfun(@(x) 100 * x.active_rois / x.num_rois, mut_data.file_summaries);
        
        hold on;
        bar(1:length(wt_activity_fractions), wt_activity_fractions, 0.4, ...
            'FaceColor', colA, 'DisplayName', lblA);
        bar((1:length(mut_activity_fractions)) + 0.4, mut_activity_fractions, 0.4, ...
            'FaceColor', colB, 'DisplayName', lblB);
        
        xlabel('File Index');
        ylabel('Active ROI Fraction (%)');
        title('ROI Recruitment per File');
        legend('Location', 'best');
        grid on;
    end
    
    nexttile(tl, 4); axis off;
    
    summary_text = {
        'FILE-LEVEL SUMMARY';
        '';
        sprintf('Dataset Size:');
        sprintf(['  ' lblA ': %d files, %d ROIs, %d events'], ...
            wt_data.num_files, wt_data.total_rois, wt_data.total_events);
        sprintf(['  ' lblB ': %d files, %d ROIs, %d events'], ...
            mut_data.num_files, mut_data.total_rois, mut_data.total_events);
        '';
        sprintf('Per-File Averages:');
        sprintf(['  ' lblA ' events/file: %.0f ± %.0f'], ...
            mean(arrayfun(@(x) x.num_events, wt_data.file_summaries)), ...
            std(arrayfun(@(x) x.num_events, wt_data.file_summaries)));
        sprintf(['  ' lblB ' events/file: %.0f ± %.0f'], ...
            mean(arrayfun(@(x) x.num_events, mut_data.file_summaries)), ...
            std(arrayfun(@(x) x.num_events, mut_data.file_summaries)));
        '';
        sprintf(['  ' lblA ' active ROIs/file: %.0f ± %.0f'], ...
            mean(arrayfun(@(x) x.active_rois, wt_data.file_summaries)), ...
            std(arrayfun(@(x) x.active_rois, wt_data.file_summaries)));
        sprintf(['  ' lblB ' active ROIs/file: %.0f ± %.0f'], ...
            mean(arrayfun(@(x) x.active_rois, mut_data.file_summaries)), ...
            std(arrayfun(@(x) x.active_rois, mut_data.file_summaries)));
        '';
        sprintf('Statistical Tests:');
        ['  Frequency: ' p_label(comparison, 'frequency_pvalue')];
        ['  Amplitude: ' p_label(comparison, 'amplitude_pvalue')];
    };
    
    text(0.05, 0.95, summary_text, 'FontSize', 10, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', 'Units', 'normalized', 'FontName', 'FixedWidth');
    
    title(tl, ['File-Level Summary: ' lblA ' vs ' lblB], 'FontSize', 16, 'FontWeight', 'bold');
end

function file_results = process_files(entries, config)
    % Run spont_pipeline once per dataset file (from memory, with that file's
    % Recordings.csv timing). A failure is reported, not silently dropped, and
    % leaves an empty slot.
    file_results = cell(1, numel(entries));
    for i = 1:numel(entries)
        try
            result = spont_pipeline(entries(i), struct('createPlots', false, ...
                'verbose', false, 'config', config));
            if isfield(result, 'fileResults') && ~isempty(result.fileResults)
                file_results{i} = result.fileResults{1};
            else
                fprintf('  [%d/%d] NO RESULT %s\n', i, numel(entries), entries(i).name);
            end
        catch ME
            fprintf('  [%d/%d] FAILED %s: %s\n', i, numel(entries), entries(i).name, ME.message);
        end
    end
end

function condition_data = aggregate_condition_data_integrated(file_results, file_paths, condition_name)
    condition_data = struct();
    condition_data.condition = condition_name;
    condition_data.num_files = length(file_results);
    
    all_roi_frequencies = []; all_individual_events = []; all_roi_event_counts = [];
    total_rois = 0; total_events = 0; file_summaries = [];
    
    for i = 1:length(file_results)
        if ~isempty(file_results{i}) && isfield(file_results{i}, 'event_stats')
            es = file_results{i}.event_stats;
            all_roi_frequencies = [all_roi_frequencies, es.frequency_hz];
            all_individual_events = [all_individual_events, es.all_event_peaks];
            all_roi_event_counts = [all_roi_event_counts, es.events_per_roi];
            total_rois = total_rois + file_results{i}.numROIs;
            total_events = total_events + es.total_events;
            
            [~, fname] = fileparts(file_paths{i});
            fs = struct('filename', fname, 'num_rois', file_results{i}.numROIs, ...
                       'num_events', es.total_events, 'active_rois', sum(es.events_per_roi > 0), ...
                       'mean_frequency', mean(es.frequency_hz(es.frequency_hz > 0)));
            if isnan(fs.mean_frequency), fs.mean_frequency = 0; end
            file_summaries = [file_summaries; fs];
        end
    end
    
    condition_data.all_roi_frequencies = all_roi_frequencies;
    condition_data.all_individual_events = all_individual_events;
    condition_data.all_roi_event_counts = all_roi_event_counts;
    condition_data.file_summaries = file_summaries;
    condition_data.total_rois = total_rois;
    condition_data.total_events = total_events;
    
    % Calculate stats
    freqs = all_roi_frequencies;
    active_freqs = freqs(freqs > 0);
    amps = all_individual_events;
    
    condition_data.stats = struct();
    condition_data.stats.frequency.active_rois.mean = mean(active_freqs);
    condition_data.stats.frequency.active_rois.std = std(active_freqs);
    condition_data.stats.frequency.active_rois.count = length(active_freqs);
    condition_data.stats.frequency.active_rois.fraction = length(active_freqs) / length(freqs);
    condition_data.stats.amplitude.mean = mean(amps);
    condition_data.stats.amplitude.std = std(amps);
    condition_data.stats.amplitude.count = length(amps);
end

function comparison = compare_conditions_statistical_integrated(wt_data, mut_data)
    [lblA, lblB, colA, colB] = spont_group_labels(wt_data, mut_data);   % group labels/colours
    comparison = struct('test_type', 'Mann-Whitney U (ranksum)', 'alpha', 0.05);
    
    wt_freqs = wt_data.all_roi_frequencies(wt_data.all_roi_frequencies > 0);
    mut_freqs = mut_data.all_roi_frequencies(mut_data.all_roi_frequencies > 0);
    
    if ~isempty(wt_freqs) && ~isempty(mut_freqs)
        [comparison.frequency_pvalue, comparison.frequency_h] = ranksum(wt_freqs, mut_freqs);
        fprintf(['    Frequency: ' lblA '=%.4f Hz, ' lblB '=%.4f Hz, p=%.6f\n'], ...
            mean(wt_freqs), mean(mut_freqs), comparison.frequency_pvalue);
    end
    
    wt_amps = wt_data.all_individual_events;
    mut_amps = mut_data.all_individual_events;
    
    if ~isempty(wt_amps) && ~isempty(mut_amps)
        [comparison.amplitude_pvalue, comparison.amplitude_h] = ranksum(wt_amps, mut_amps);
        fprintf(['    Amplitude: ' lblA '=%.4f dF/F, ' lblB '=%.4f dF/F, p=%.6f\n'], ...
            mean(wt_amps), mean(mut_amps), comparison.amplitude_pvalue);
    end
    
    comparison.summary = struct();
    comparison.summary.wt_n_files = wt_data.num_files;
    comparison.summary.mut_n_files = mut_data.num_files;
    comparison.summary.wt_n_rois = wt_data.total_rois;
    comparison.summary.mut_n_rois = mut_data.total_rois;
    comparison.summary.wt_n_events = wt_data.total_events;
    comparison.summary.mut_n_events = mut_data.total_events;
    
    p_values = [comparison.frequency_pvalue, comparison.amplitude_pvalue];
    comparison.summary.significant_tests = sum(p_values < 0.05, 'omitnan');
    comparison.summary.total_tests = sum(~isnan(p_values));
end

function print_condition_summary_integrated(results)
    [lblA, lblB, colA, colB] = spont_group_labels(results);   % group labels/colours
    fprintf('\n=== INTEGRATED CONDITION COMPARISON SUMMARY ===\n');
    wt = results.wt_data; mut = results.mut_data; comp = results.comparison;
    
    fprintf('Dataset Overview:\n');
    fprintf(['  ' lblA ': %d files, %d ROIs, %d events\n'], wt.num_files, wt.total_rois, wt.total_events);
    fprintf(['  ' lblB ': %d files, %d ROIs, %d events\n'], mut.num_files, mut.total_rois, mut.total_events);
    
    fprintf('\nEvent Frequency (Active ROIs only):\n');
    fprintf(['  ' lblA ': %.4f ± %.4f Hz (n=%d ROIs)\n'], ...
        wt.stats.frequency.active_rois.mean, wt.stats.frequency.active_rois.std, wt.stats.frequency.active_rois.count);
    fprintf(['  ' lblB ': %.4f ± %.4f Hz (n=%d ROIs)\n'], ...
        mut.stats.frequency.active_rois.mean, mut.stats.frequency.active_rois.std, mut.stats.frequency.active_rois.count);
    
    if isfield(comp, 'frequency_pvalue')
        fprintf('  Mann-Whitney U (ROI-pooled): %s\n', p_label(comp, 'frequency_pvalue'));
    end
    
    fprintf('\nActive ROI Fraction:\n');
    fprintf(['  ' lblA ': %.1f%% (%d/%d ROIs)\n'], ...
        100 * wt.stats.frequency.active_rois.fraction, wt.stats.frequency.active_rois.count, wt.total_rois);
    fprintf(['  ' lblB ': %.1f%% (%d/%d ROIs)\n'], ...
        100 * mut.stats.frequency.active_rois.fraction, mut.stats.frequency.active_rois.count, mut.total_rois);
    
    fprintf('\n===============================================\n');
end

function save_condition_plots_integrated(plot_handles, folder_path)
    % Optional plot saving
end

function txt = p_label(comparison, field)
    % p-value text from the actual test result (never a fixed string).
    if ~isfield(comparison, field) || isnan(comparison.(field))
        txt = 'p = n/a';
    elseif comparison.(field) < 0.001
        txt = 'p < 0.001 ***';
    elseif comparison.(field) < 0.01
        txt = sprintf('p = %.3f **', comparison.(field));
    elseif comparison.(field) < 0.05
        txt = sprintf('p = %.3f *', comparison.(field));
    else
        txt = sprintf('p = %.3f (ns)', comparison.(field));
    end
end

function value = getfield_safe(s, field, default)
    % Safely get field with default value
    if isfield(s, field)
        value = s.(field);
    else
        value = default;
    end
end

function bins = activity_categories(edges)
    % {label, @(x) in-bin} rows from ascending integer edges: bin k holds
    % edges(k) .. edges(k+1)-1 events, the last bin edges(end) or more.
    edges = edges(:)';
    if isempty(edges) || any(diff(edges) <= 0) || any(edges ~= round(edges)) || edges(1) < 1
        error('spont_batch_analysis:activityBins', ...
            'config.visualization.spont.activity_bin_edges must be ascending integers >= 1.');
    end
    bins = cell(numel(edges), 2);
    for k = 1:numel(edges)
        lo = edges(k);
        if k == numel(edges)
            bins(k, :) = {sprintf('%d+ events', lo), @(x) x >= lo};
        else
            hi = edges(k + 1) - 1;
            if hi == lo
                lbl = sprintf('%d event%s', lo, repmat('s', 1, lo ~= 1));
            else
                lbl = sprintf('%d-%d events', lo, hi);
            end
            bins(k, :) = {lbl, @(x) x >= lo & x <= hi};
        end
    end
end
