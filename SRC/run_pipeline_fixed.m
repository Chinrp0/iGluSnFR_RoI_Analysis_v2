%% ========================================================================
%% RUN_PIPELINE_FIXED - WT vs R213W spontaneous activity analysis
%% ========================================================================
% Single entry point for the analysis.
%
% integrated_batch_analysis_final() already creates ALL plots (the 4 core
% comparison plots in its Step 5 + the 7 additional plots in its Step 6).
% This script therefore does NOT recreate any plots. Recreating them in a
% separate Step 4 was the cause of the duplicate-figure problem
% (~26 figures instead of ~12).
%
% Acquisition parameters (frame_rate, recording_duration_s) are read from
% tracenorm_config() so there is ONE source of truth shared by the event
% detectors (inside main_pipeline) and the temporal plot modules (IEI,
% raster, cumulative). If those two ever disagree, frequencies and timing
% axes silently land on different clocks.
%
% To change acquisition rate, edit tracenorm_config.m only.
%% ========================================================================

clear; close all; clc;

setup_pipeline

%% --- Single source of truth: mode, exposure, groups AND data folder --
% Analysis mode ('1AP', '2AP' or 'Spont'), exposure and the data folder are all
% set in tracenorm_config.m (config.dataType, config.exposure_ms,
% config.folder_path). That one file selects the whole analysis path below.
cfg = tracenorm_config();
folder_path = cfg.folder_path;

options = struct();
options.frame_rate           = cfg.frame_rate;             % Hz
options.recording_duration_s = cfg.recording_duration_s;   % seconds
options.verbose              = true;
options.createPlots          = true;

fprintf('\nMode: %s | Acquisition (from tracenorm_config): %d Hz, %.1f s (%d frames expected)\n', ...
    cfg.dataType, cfg.frame_rate, cfg.recording_duration_s, cfg.expected_frames);

%% --- Run analysis ----------------------------------------------------
fprintf('\n=== RUNNING ANALYSIS ===\n');
if strcmpi(cfg.dataType, '1AP')
    results = oneap_batch_analysis(folder_path, options);
elseif strcmpi(cfg.dataType, '2AP')
    % Paired-pulse (PPF): PPR vs ISI (per coverslip) + method-comparison montage
    results = twoap_batch_analysis(folder_path, options);
else
    % Spontaneous: creates core + all 7 additional plots
    results = integrated_batch_analysis_final(folder_path, options);
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
    metrics_csv = fullfile(output_dir, 'oneap_per_roi_metrics.csv');
    try
        writetable(results.oneap_metrics_table, metrics_csv);
        fprintf('Saved per-ROI 1AP metrics to:\n  %s\n', metrics_csv);
    catch ME
        fprintf('Could not write metrics table: %s\n', ME.message);
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