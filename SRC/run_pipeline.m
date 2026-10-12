%% ========================================================================
%% RUN_PIPELINE - single entry point for 1AP / 2AP / Spont analysis
%% ========================================================================
% Analysis mode comes from tracenorm_config.m (config.dataType). The data
% folder comes from the loader:
%   - If `dataset` is in the workspace (from load_dataset or a previous run),
%     it is used as-is: no dialog, no re-reading of the CSVs.
%   - Otherwise a folder dialog opens at config.data_root, and for 1AP the
%     folder (data + Recordings.csv) is loaded into `dataset`, which stays in
%     the workspace for the next run.
% To switch folders: clear dataset (or run load_dataset again).
%
% 1AP: timing (frames, frame interval, stim frames) comes from Recordings.csv.
% 2AP / Spont: only the folder comes from `dataset`; they still read their
% files themselves with the legacy exposure in tracenorm_config.m.
%
% spont_batch_analysis() (Spont) already creates ALL its plots, so
% this script does not recreate any (that caused duplicate figures before).
%% ========================================================================

clearvars -except dataset; close all; clc;

setup_pipeline

%% --- Mode + groups from tracenorm_config; folder from the loader -----
cfg = tracenorm_config();

has_dataset = exist('dataset', 'var') && isstruct(dataset) && isfield(dataset, 'files');
if has_dataset
    fprintf('Using dataset in workspace (loaded %s):\n  %s\n', ...
        char(dataset.loaded_at), dataset.folder);
    fprintf('  (run "clear dataset" to pick a different folder)\n');
    folder_path = dataset.folder;
else
    folder_path = select_data_folder(cfg);
end

options = struct();
options.verbose              = true;
options.createPlots          = true;

%% --- Run analysis ----------------------------------------------------
fprintf('\n=== RUNNING ANALYSIS ===\n');
if strcmpi(cfg.dataType, '1AP')
    if ~has_dataset
        dataset = load_dataset_folder(folder_path, cfg);
    end
    results = oneap_batch_analysis(dataset, options);
elseif strcmpi(cfg.dataType, '2AP')
    % Paired-pulse (PPF): PPR vs ISI (per coverslip) + method-comparison montage
    % (legacy timing: exposure_ms / stim1_frame from tracenorm_config)
    results = twoap_batch_analysis(folder_path, options);
else
    % Spontaneous: genotype x drug groups; core + 7 additional plots per
    % genotype pair within each drug (legacy timing)
    options.frame_rate           = cfg.frame_rate;             % Hz
    options.recording_duration_s = cfg.recording_duration_s;   % seconds
    results = spont_batch_analysis(folder_path, options);
end

fprintf('\n=== Analysis complete ===\n');
fprintf('Open figures: %d\n', numel(findall(0, 'Type', 'figure')));

%% --- Save all open figures -------------------------------------------
output_dir = make_figures_dir(folder_path);
fprintf('Created output directory:\n  %s\n', output_dir);
n_saved = save_open_figures(output_dir);
fprintf('Saved %d figures.\n', n_saved);

%% --- Export per-ROI metrics table -------------------------------------
if strcmpi(cfg.dataType, '1AP') && isfield(results, 'oneap_metrics_table')
    tables = {results.oneap_metrics_table, 'oneap_per_trial_metrics.csv', 'per-ROI-per-trial'; ...
              results.oneap_roi_table,     'oneap_per_roi_release.csv',   'per-ROI release'; ...
              dataset.summary,             'dataset_validation.csv',      'dataset validation'};
    for k = 1:size(tables, 1)
        out_csv = fullfile(output_dir, tables{k, 2});
        try
            writetable(tables{k, 1}, out_csv);
            fprintf('Saved %s table to:\n  %s\n', tables{k, 3}, out_csv);
        catch ME
            fprintf('Could not write %s: %s\n', tables{k, 2}, ME.message);
        end
    end
elseif strcmpi(cfg.dataType, '2AP') && isfield(results, 'twoap_metrics_table')
    metrics_csv = fullfile(output_dir, 'twoap_per_roi_metrics.csv');
    try
        writetable(results.twoap_metrics_table, metrics_csv);
        fprintf('Saved per-ROI 2AP metrics to:\n  %s\n', metrics_csv);
    catch ME
        fprintf('Could not write metrics table: %s\n', ME.message);
    end
end

%% --- Report ----------------------------------------------------------
fprintf('\n=== PIPELINE COMPLETE ===\n');
list_open_figures();
fprintf('Figures saved to:\n  %s\n', output_dir);
fprintf('To close all figures: close all\n\n');