%% ========================================================================
%% RUN_PIPELINE - single entry point for 1AP / 2AP / HFS / Spont analysis
%% ========================================================================
% Everything about the recordings comes from the data folder's Recordings.csv:
% the paradigm (Experiment column), frame interval, frame counts and stim
% frames. The folder comes from the loader:
%   - If `dataset` is in the workspace (from load_dataset or a previous run),
%     it is used as-is: no dialog, no re-reading of the CSVs.
%   - Otherwise a folder dialog opens at config.data_root and the folder is
%     loaded into `dataset`.
% After the run, `dataset` stays in the workspace for the next run when
% config.keep_dataset is true, or is cleared when it is false (the next run
% then asks for a folder again). If a run errors, `dataset` is kept either
% way so it can be re-run without reloading.
% To switch folders: clear dataset (or run load_dataset again).
%
% Paradigm: config.dataType = 'auto' uses the Experiment column (a folder that
% mixes paradigms is an error); set '1AP' | '2AP' | 'HFS' | 'Spont' to insist on one.
%
% Files recorded at different frame intervals (e.g. 5 ms and 10 ms PPF in one
% folder) are analysed as separate sets; with more than one set, each set's
% figures and tables go in its own sub-folder (e.g. ...\10ms\).
%
% spont_batch_analysis() (Spont) already creates ALL its plots, so
% this script does not recreate any (that caused duplicate figures before).
%% ========================================================================

clearvars -except dataset; close all; clc;

setup_pipeline

%% --- Groups from tracenorm_config; folder + timing from the loader ----
cfg = tracenorm_config();

if exist('dataset', 'var') && isstruct(dataset) && isfield(dataset, 'files')
    fprintf('Using dataset in workspace (loaded %s):\n  %s\n', ...
        char(dataset.loaded_at), dataset.folder);
    fprintf('  (run "clear dataset" to pick a different folder)\n');
else
    dataset = load_dataset_folder(select_data_folder(cfg), cfg);
end
if isempty(dataset.files)
    error('run_pipeline:noFiles', 'No loadable files in %s (see dataset.summary).', dataset.folder);
end

%% --- Paradigm (from the metadata unless the config insists) -----------
if numel(dataset.dataType) > 1
    error('run_pipeline:mixedParadigm', ...
        'The folder mixes paradigms (%s); analyse them from separate folders.', ...
        strjoin(dataset.dataType, ', '));
end
dataType = dataset.dataType{1};
if ~strcmpi(cfg.dataType, 'auto') && ~strcmpi(cfg.dataType, dataType)
    error('run_pipeline:paradigm', ...
        'tracenorm_config.dataType is %s but Recordings.csv says this folder is %s.', ...
        cfg.dataType, dataType);
end
fprintf('Paradigm (from %s): %s\n', cfg.metadata.filename, dataType);

options = struct();
options.verbose              = true;
options.createPlots          = true;

%% --- Run analysis, one set per frame interval ------------------------
parts = split_dataset_by_interval(dataset, cfg);
base_dir = make_figures_dir(dataset.folder);
fprintf('Output directory:\n  %s\n', base_dir);
results_all = cell(1, numel(parts));

for p = 1:numel(parts)
    part = parts(p);
    fprintf('\n=== RUNNING %s ANALYSIS (%s, %d files) ===\n', dataType, part.label, numel(part.files));
    switch dataType
        case '1AP'
            res = oneap_batch_analysis(part, options);
            tables = {res.oneap_metrics_table, 'oneap_per_trial_metrics.csv'; ...
                      res.oneap_roi_table,     'oneap_per_roi_release.csv'; ...
                      res.oneap_coverslip_table, 'oneap_per_coverslip.csv'};
        case '2AP'
            % Paired-pulse (PPF): PPR vs ISI (per coverslip) + method-comparison montage
            res = twoap_batch_analysis(part, options);
            tables = {res.twoap_metrics_table,   'twoap_per_roi_metrics.csv'; ...
                      res.twoap_coverslip_table, 'twoap_per_coverslip.csv'; ...
                      res.twoap_sync_async_table, 'twoap_sync_async_per_coverslip.csv'; ...
                      res.isi_checks,            'twoap_isi_checks.csv'};
        case 'HFS'
            % High-frequency train: per-pulse phasic / tonic + train charge, plateau,
            % decay, from per-coverslip average traces
            res = hfs_batch_analysis(part, options);
            tables = {res.hfs_roi_table,       'hfs_per_roi_metrics.csv'; ...
                      res.hfs_pulse_table,     'hfs_per_pulse.csv'; ...
                      res.hfs_coverslip_table, 'hfs_per_coverslip.csv'; ...
                      res.file_checks,         'hfs_file_checks.csv'};
            % No per-pulse table when single pulses are not resolvable (e.g. 25 ms frames)
            tables = tables(~cellfun(@(t) width(t) == 0, tables(:, 1)), :);
        case 'Spont'
            % Spontaneous: genotype x drug groups; core + 7 additional plots per
            % genotype pair within each drug
            res = spont_batch_analysis(part, options);
            tables = {};
    end
    results_all{p} = res;
    if cfg.visualization.y_min_zero
        apply_ymin_zero();
    end
    info = [];
    if cfg.visualization.show_experiment_info
        info = experiment_info(dataset.folder, part.files, cfg);
    end
    frame_ms = [];
    if cfg.visualization.show_frame_time
        frame_ms = median(arrayfun(@(f) f.meta.frame_interval_ms, part.files));
    end
    annotate_figure_info(info, frame_ms);

    % Each frame-interval set gets its own sub-folder when there are several.
    if numel(parts) > 1
        output_dir = fullfile(base_dir, part.label);
        if ~isfolder(output_dir), mkdir(output_dir); end
    else
        output_dir = base_dir;
    end
    if cfg.output.save_figures
        n_saved = save_open_figures(output_dir, cfg);
        fprintf('Saved %d figures (%s) to:\n  %s\n', n_saved, cfg.output.figure_format, output_dir);
    end
    if ~cfg.output.save_data, tables = {}; end
    for k = 1:size(tables, 1)
        out_csv = fullfile(output_dir, tables{k, 2});
        try
            writetable(tables{k, 1}, out_csv);
            fprintf('Saved %s\n', tables{k, 2});
        catch ME
            fprintf('Could not write %s: %s\n', tables{k, 2}, ME.message);
        end
    end
    if p < numel(parts)
        close all;   % keep only the last set's figures open
    end
end
writetable(dataset.summary, fullfile(base_dir, 'dataset_validation.csv'));
results = [results_all{:}];   % one element per frame-interval set

%% --- Report ----------------------------------------------------------
fprintf('\n=== PIPELINE COMPLETE ===\n');
list_open_figures();
fprintf('Figures and tables saved to:\n  %s\n', base_dir);
fprintf('To close all figures: close all\n\n');

%% --- Keep or clear the loaded data (config.keep_dataset) -------------
if cfg.keep_dataset
    fprintf('Dataset kept in workspace (config.keep_dataset = true); the next run reuses it.\n');
else
    clear dataset
    fprintf('Dataset cleared (config.keep_dataset = false); the next run asks for a folder.\n');
end
