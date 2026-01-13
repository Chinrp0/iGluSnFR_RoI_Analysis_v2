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

%% Step 5: Summary
fprintf('\n=== PIPELINE COMPLETE ===\n');
fprintf('Final figure count: %d\n', length(findall(0, 'Type', 'figure')));
fprintf('\nTo see a list of all figures: list_open_figures()\n');
fprintf('To close all figures: close all\n\n');

%% Optional: List all open figures
fprintf('Listing all open figures:\n');
list_open_figures();
