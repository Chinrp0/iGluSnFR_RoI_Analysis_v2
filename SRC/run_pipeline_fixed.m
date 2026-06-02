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

%% --- Data folder ------------------------------------------------------
% EDIT THIS to point at your CSV folder:
folder_path = 'E:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglu3fast_NGR\1AP\GPU_Processed_Images_1AP\5_raw_mean';


%% --- Single source of truth for acquisition parameters ---------------
% Analysis mode ('Spont' or '1AP') is set in tracenorm_config.m via
% config.dataType. That one switch selects the whole analysis path below.
cfg = tracenorm_config();

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
else
    % Spontaneous: creates core + all 7 additional plots
    results = integrated_batch_analysis_final(folder_path, options);
end

fprintf('\n=== Analysis complete ===\n');
fprintf('Open figures: %d\n', numel(findall(0, 'Type', 'figure')));

%% --- Save all open figures -------------------------------------------
output_dir = make_output_dir(folder_path);
save_all_figures(output_dir);

%% --- 1AP: export per-ROI metrics table -------------------------------
if strcmpi(cfg.dataType, '1AP') && isfield(results, 'oneap_metrics_table')
    metrics_csv = fullfile(output_dir, 'oneap_per_roi_metrics.csv');
    try
        writetable(results.oneap_metrics_table, metrics_csv);
        fprintf('Saved per-ROI 1AP metrics to:\n  %s\n', metrics_csv);
    catch ME
        fprintf('Could not write metrics table: %s\n', ME.message);
    end
end

%% --- Report ----------------------------------------------------------
fprintf('\n=== PIPELINE COMPLETE ===\n');
list_open_figures();
fprintf('Figures saved to:\n  %s\n', output_dir);
fprintf('To close all figures: close all\n\n');


%% ========================================================================
%% Local helper functions
%% ========================================================================
function output_dir = make_output_dir(folder_path)
    % Create a dated, auto-incremented output directory one level above the
    % data folder.
    [parent, data_name] = fileparts(folder_path);
    [grandparent, ~]     = fileparts(parent);

    date_str = datestr(now, 'yyyymmdd');
    base     = sprintf('%s_figures_%s', data_name, date_str);

    run_num = 1;
    while exist(fullfile(grandparent, sprintf('%s_%d', base, run_num)), 'dir')
        run_num = run_num + 1;
    end

    output_dir = fullfile(grandparent, sprintf('%s_%d', base, run_num));
    mkdir(output_dir);
    fprintf('Created output directory (run %d):\n  %s\n', run_num, output_dir);
end

function save_all_figures(output_dir)
    % Save every open figure as both PNG (viewing) and FIG (editing).
    figs = findall(0, 'Type', 'figure');
    fprintf('Saving %d figures...\n', numel(figs));

    saved = 0;
    for i = 1:numel(figs)
        fig = figs(i);

        if ~isempty(fig.Name)
            name = fig.Name;
        else
            name = sprintf('Figure_%d', fig.Number);
        end

        % Clean filename
        name = strrep(name, ' ', '_');
        name = regexprep(name, '[^\w\-]', '');
        stem = sprintf('%d_%s', fig.Number, name);

        try
            saveas(fig, fullfile(output_dir, [stem '.png']));
            saveas(fig, fullfile(output_dir, [stem '.fig']));
            saved = saved + 1;
        catch ME
            fprintf('  Could not save %s: %s\n', stem, ME.message);
        end
    end

    fprintf('Saved %d/%d figures.\n', saved, numel(figs));
end