%% ========================================================================
%% IMPROVED RUN_PIPELINE - Fixed frame_rate issues
%% ========================================================================
% This version properly passes frame_rate to all plotting modules
% Run this instead of the original run_pipeline.m

%% Step 1: Clean workspace and setup
clear all
close all  % Close any existing figures
clc

setup_pipeline

%% Step 2: Set your data folder path
folder_path = 'E:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglusnfr4f_NGR\Spont\Spont_2025-12-29--1\E_raw_mean';

%% Step 3: Run the integrated batch analysis (creates core plots)
fprintf('\n=== RUNNING MAIN ANALYSIS ===\n');
results = integrated_batch_analysis_final(folder_path);

fprintf('\n=== Main analysis complete ===\n');
fprintf('Initial figure count: %d\n', length(findall(0, 'Type', 'figure')));

%% Step 4: Create additional plots WITH proper frame_rate parameter
fprintf('\n=== CREATING ADDITIONAL PLOTS ===\n');

% Set up options with frame_rate
options = struct();
if isfield(results, 'processing_info') && isfield(results.processing_info, 'frame_rate')
    options.frame_rate = results.processing_info.frame_rate;
else
    options.frame_rate = 100;  % Default
end

% Create plots that failed due to missing frame_rate
fprintf('Creating plots with frame_rate = %d Hz...\n', options.frame_rate);

try
    fprintf('  1. Inter-Event Interval comparison...\n');
    fig_iei = plot_iei_comparison(results, options);
    fprintf('     ✓ Created (Figure %d)\n', fig_iei.Number);
catch ME
    fprintf('     ✗ Failed: %s\n', ME.message);
end

try
    fprintf('  2. Amplitude-Frequency correlation...\n');
    fig_ampfreq = plot_amplitude_frequency_correlation(results, options);
    fprintf('     ✓ Created (Figure %d)\n', fig_ampfreq.Number);
catch ME
    fprintf('     ✗ Failed: %s\n', ME.message);
end

try
    fprintf('  3. Event Timing raster...\n');
    fig_raster = plot_event_raster(results, options);
    fprintf('     ✓ Created (Figure %d)\n', fig_raster.Number);
catch ME
    fprintf('     ✗ Failed: %s\n', ME.message);
end

try
    fprintf('  4. Cumulative Events over time...\n');
    fig_cumulative = plot_cumulative_events(results, options);
    fprintf('     ✓ Created (Figure %d)\n', fig_cumulative.Number);
catch ME
    fprintf('     ✗ Failed: %s\n', ME.message);
end

try
    fprintf('  5. ROI Trace comparison...\n');
    fig_traces = plot_condition_roi_traces(results, options);
    fprintf('     ✓ Created (Figure %d)\n', fig_traces.Number);
catch ME
    fprintf('     ✗ Failed: %s\n', ME.message);
end

%% Step 5: Save all figures
fprintf('\n=== SAVING FIGURES ===\n');

% Create output directory based on data folder
[parent_folder, data_folder_name] = fileparts(folder_path);
output_dir = fullfile(parent_folder, [data_folder_name, '_figures']);

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Created output directory: %s\n', output_dir);
else
    fprintf('Using existing output directory: %s\n', output_dir);
end

% Get all open figures
all_figures = findall(0, 'Type', 'figure');
fprintf('Found %d figures to save\n', length(all_figures));

% Save each figure
saved_count = 0;
for i = 1:length(all_figures)
    fig = all_figures(i);

    % Get figure name (use figure title or default name)
    if ~isempty(fig.Name)
        fig_name = fig.Name;
    elseif ~isempty(fig.Children) && isprop(fig.Children(1), 'Title')
        fig_name = get(fig.Children(1).Title, 'String');
        if iscell(fig_name)
            fig_name = fig_name{1};
        end
    else
        fig_name = sprintf('Figure_%d', fig.Number);
    end

    % Clean filename (remove special characters)
    fig_name = strrep(fig_name, ' ', '_');
    fig_name = strrep(fig_name, ':', '-');
    fig_name = strrep(fig_name, '/', '_');
    fig_name = strrep(fig_name, '\', '_');
    fig_name = regexprep(fig_name, '[^\w\-]', '');

    % Save as both PNG (for viewing) and FIG (for editing)
    png_path = fullfile(output_dir, [fig_name, '.png']);
    fig_path = fullfile(output_dir, [fig_name, '.fig']);

    try
        % Save as PNG (high resolution)
        saveas(fig, png_path);
        fprintf('  ✓ Saved: %s.png\n', fig_name);

        % Save as FIG (MATLAB figure file)
        saveas(fig, fig_path);
        fprintf('  ✓ Saved: %s.fig\n', fig_name);

        saved_count = saved_count + 1;
    catch ME
        fprintf('  ✗ Failed to save %s: %s\n', fig_name, ME.message);
    end
end

fprintf('\nSuccessfully saved %d/%d figures to:\n  %s\n', saved_count, length(all_figures), output_dir);

%% Step 6: Summary
fprintf('\n=== PIPELINE COMPLETE ===\n');
fprintf('Final figure count: %d\n', length(findall(0, 'Type', 'figure')));
fprintf('Figures saved to: %s\n', output_dir);
fprintf('\nTo see a list of all figures: list_open_figures()\n');
fprintf('To close all figures: close all\n\n');

%% Optional: List all open figures
fprintf('Listing all open figures:\n');
list_open_figures();
