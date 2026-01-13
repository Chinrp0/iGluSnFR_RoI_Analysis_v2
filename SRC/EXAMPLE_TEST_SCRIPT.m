%% ========================================================================
%% EXAMPLE SCRIPT: Testing New Visualization Modules
%% ========================================================================
%% This script demonstrates how to use all 6 new visualization modules
%% Copy and paste sections as needed for your analysis
%% ========================================================================

%% Section 1: Setup and Configuration
% Navigate to your SRC directory
cd('C:\Users\xdach\OneDrive - Johns Hopkins\Maher_Lab\Protocols\Matlab_scripts\gluSnFR\iGluSnFR_RoI_Analysis_v2\iGluSnFR_RoI_Analysis_v2\SRC')

% Run setup to add all modules to path
setup_pipeline

% Define your data folder
folder_path = 'D:\Data\GluSnFR\Ms\2025-10-25_CP_Ms-hipp_Snfr4f-NGR_Doc2b_1mMCa_resave\spont\GPU_SNR_Processed\5_raw_mean';

%% Section 2: Run Analysis with ALL Visualizations (AUTOMATIC MODE)
fprintf('\n=== Running COMPLETE analysis with all plots ===\n');

% Configure options for 10ms exposure data
options = struct();
options.frame_rate = 100;  % 10ms exposure = 100 Hz
options.recording_duration_s = 30;  % 3000 frames / 100 Hz = 30s
options.verbose = true;
options.createPlots = true;

% Enable all new visualizations
options.create_iei_plot = true;
options.create_amp_freq_correlation = true;
options.create_event_raster = true;
options.create_biological_variability = true;
options.create_recruitment_curves = true;
options.create_cumulative_events = true;
options.create_roi_traces = true;

% Run the complete integrated analysis
tic;
results = integrated_batch_analysis_final(folder_path, options);
analysis_time = toc;

fprintf('\nTotal analysis time: %.1f seconds\n', analysis_time);
fprintf('Created %d core plots + %d additional plots\n', ...
    length(fieldnames(results.plot_handles)), ...
    length(fieldnames(results.additional_plots)));

%% Section 3: Access and Inspect Results
fprintf('\n=== Inspecting Results ===\n');

% Core statistics
fprintf('WT Statistics:\n');
fprintf('  Files: %d\n', results.wt_data.num_files);
fprintf('  Total ROIs: %d\n', results.wt_data.total_rois);
fprintf('  Active ROIs: %d (%.1f%%)\n', ...
    results.wt_data.stats.frequency.active_rois.count, ...
    100 * results.wt_data.stats.frequency.active_rois.fraction);
fprintf('  Mean frequency: %.4f ± %.4f Hz\n', ...
    results.wt_data.stats.frequency.active_rois.mean, ...
    results.wt_data.stats.frequency.active_rois.std);
fprintf('  Total events: %d\n', results.wt_data.total_events);

fprintf('\nR213W Statistics:\n');
fprintf('  Files: %d\n', results.mut_data.num_files);
fprintf('  Total ROIs: %d\n', results.mut_data.total_rois);
fprintf('  Active ROIs: %d (%.1f%%)\n', ...
    results.mut_data.stats.frequency.active_rois.count, ...
    100 * results.mut_data.stats.frequency.active_rois.fraction);
fprintf('  Mean frequency: %.4f ± %.4f Hz\n', ...
    results.mut_data.stats.frequency.active_rois.mean, ...
    results.mut_data.stats.frequency.active_rois.std);
fprintf('  Total events: %d\n', results.mut_data.total_events);

fprintf('\nStatistical Comparison:\n');
fprintf('  Frequency p-value: %.6f\n', results.comparison.frequency_pvalue);
fprintf('  Amplitude p-value: %.6f\n', results.comparison.amplitude_pvalue);

%% Section 4: Access Individual Plot Handles
fprintf('\n=== Accessing Plot Handles ===\n');

% Core plots
core_plots = results.plot_handles;
fprintf('Core plots available:\n');
disp(fieldnames(core_plots));

% Additional plots
additional_plots = results.additional_plots;
fprintf('\nAdditional plots available:\n');
disp(fieldnames(additional_plots));

%% Section 5: Standalone Mode - Create Individual Plots
fprintf('\n=== Creating Individual Plots (STANDALONE MODE) ===\n');

% First run analysis WITHOUT plots (if you haven't already)
% options_no_plots = struct('createPlots', false, 'frame_rate', 100);
% results = integrated_batch_analysis_final(folder_path, options_no_plots);

% Then create individual plots as needed
plot_options = struct('frame_rate', 100, 'recording_duration_s', 30);

fprintf('Creating IEI plot...\n');
fig_iei = plot_iei_comparison(results, plot_options);

fprintf('Creating amplitude-frequency correlation...\n');
fig_amp_freq = plot_amplitude_frequency_correlation(results, plot_options);

fprintf('Creating event raster...\n');
fig_raster = plot_event_raster(results, plot_options);

fprintf('Creating biological variability analysis...\n');
fig_variability = plot_biological_variability(results);

fprintf('Creating ROI recruitment curves...\n');
fig_recruitment = plot_roi_recruitment_curves(results);

fprintf('Creating cumulative events over time...\n');
fig_cumulative = plot_cumulative_events_over_time(results, plot_options);

fprintf('All standalone plots created successfully!\n');

%% Section 6: Selective Plotting - Choose Specific Plots
fprintf('\n=== Creating ONLY Selected Plots ===\n');

% Configure to create only specific plots
options_selective = struct();
options_selective.frame_rate = 100;
options_selective.recording_duration_s = 30;
options_selective.verbose = true;
options_selective.createPlots = true;

% Choose which plots to create
options_selective.create_iei_plot = true;  % YES
options_selective.create_amp_freq_correlation = true;  % YES
options_selective.create_event_raster = false;  % NO (skip - memory intensive)
options_selective.create_biological_variability = true;  % YES
options_selective.create_recruitment_curves = false;  % NO (skip)
options_selective.create_cumulative_events = true;  % YES
options_selective.create_roi_traces = true;  % YES

% Run with selective plotting
results_selective = integrated_batch_analysis_final(folder_path, options_selective);

fprintf('Created %d additional plots (selective mode)\n', ...
    length(fieldnames(results_selective.additional_plots)));

%% Section 7: Customizing Plot Options
fprintf('\n=== Creating Plots with Custom Options ===\n');

% IEI plot with custom options
iei_options = struct();
iei_options.frame_rate = 100;
iei_options.max_iei_s = 15;  % Show IEIs up to 15 seconds
iei_options.bin_size_ms = 100;  % 100ms bins instead of default 50ms

fig_iei_custom = plot_iei_comparison(results, iei_options);
fprintf('Created IEI plot with custom options\n');

% Raster plot with custom options
raster_options = struct();
raster_options.frame_rate = 100;
raster_options.recording_duration_s = 30;
raster_options.max_rois_per_plot = 50;  % Show only 50 ROIs (less memory)

fig_raster_custom = plot_event_raster(results, raster_options);
fprintf('Created raster plot with 50 ROIs\n');

%% Section 8: Saving Results and Plots
fprintf('\n=== Saving Results ===\n');

% Create output directory
output_dir = fullfile(folder_path, 'analysis_results');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Save results structure
save(fullfile(output_dir, 'analysis_results.mat'), 'results', '-v7.3');
fprintf('Saved results to: %s\n', fullfile(output_dir, 'analysis_results.mat'));

% Save individual figures
figure_names = fieldnames(results.additional_plots);
for i = 1:length(figure_names)
    fig_handle = results.additional_plots.(figure_names{i});
    if ishandle(fig_handle)
        saveas(fig_handle, fullfile(output_dir, sprintf('%s.png', figure_names{i})));
        saveas(fig_handle, fullfile(output_dir, sprintf('%s.fig', figure_names{i})));
    end
end

fprintf('Saved %d figures to: %s\n', length(figure_names), output_dir);

%% Section 9: Quick Data Export for External Analysis
fprintf('\n=== Exporting Key Metrics ===\n');

% Export frequency data
wt_freqs = results.wt_data.all_roi_frequencies(results.wt_data.all_roi_frequencies > 0);
mut_freqs = results.mut_data.all_roi_frequencies(results.mut_data.all_roi_frequencies > 0);

writematrix([wt_freqs(:); mut_freqs(:)], ...
    fullfile(output_dir, 'frequencies.csv'));

% Export amplitude data
wt_amps = results.wt_data.all_individual_events;
mut_amps = results.mut_data.all_individual_events;

writematrix([wt_amps(:); mut_amps(:)], ...
    fullfile(output_dir, 'amplitudes.csv'));

fprintf('Exported frequency and amplitude data\n');

%% Section 10: Generating Publication Summary
fprintf('\n=== Publication Summary ===\n');

% Create a summary table
summary_table = table();
summary_table.Condition = {'WT'; 'R213W'};
summary_table.Files = [results.wt_data.num_files; results.mut_data.num_files];
summary_table.TotalROIs = [results.wt_data.total_rois; results.mut_data.total_rois];
summary_table.ActiveROIs = [results.wt_data.stats.frequency.active_rois.count; ...
                            results.mut_data.stats.frequency.active_rois.count];
summary_table.ActivePercent = [100 * results.wt_data.stats.frequency.active_rois.fraction; ...
                               100 * results.mut_data.stats.frequency.active_rois.fraction];
summary_table.MeanFreq_Hz = [results.wt_data.stats.frequency.active_rois.mean; ...
                             results.mut_data.stats.frequency.active_rois.mean];
summary_table.StdFreq_Hz = [results.wt_data.stats.frequency.active_rois.std; ...
                            results.mut_data.stats.frequency.active_rois.std];
summary_table.MeanAmp_dFoF = [results.wt_data.stats.amplitude.mean; ...
                              results.mut_data.stats.amplitude.mean];
summary_table.StdAmp_dFoF = [results.wt_data.stats.amplitude.std; ...
                             results.mut_data.stats.amplitude.std];
summary_table.TotalEvents = [results.wt_data.total_events; ...
                             results.mut_data.total_events];

% Display summary
disp(summary_table);

% Save summary table
writetable(summary_table, fullfile(output_dir, 'summary_table.csv'));
fprintf('Saved summary table to: %s\n', fullfile(output_dir, 'summary_table.csv'));

%% Section 11: Memory-Efficient Analysis (Large Datasets)
fprintf('\n=== Memory-Efficient Mode ===\n');

% If you have limited RAM, skip the raster plot
options_lowmem = struct();
options_lowmem.frame_rate = 100;
options_lowmem.recording_duration_s = 30;
options_lowmem.verbose = false;  % Less console output
options_lowmem.createPlots = true;
options_lowmem.create_event_raster = false;  % SKIP (uses most memory)
options_lowmem.create_cumulative_events = true;
options_lowmem.create_iei_plot = true;
options_lowmem.create_amp_freq_correlation = true;
options_lowmem.create_biological_variability = true;
options_lowmem.create_recruitment_curves = true;

fprintf('Running memory-efficient analysis (no raster)...\n');
results_lowmem = integrated_batch_analysis_final(folder_path, options_lowmem);

%% Section 12: Troubleshooting Tests
fprintf('\n=== Running Troubleshooting Tests ===\n');

% Test 1: Check that all new functions are accessible
functions_to_test = {
    'plot_iei_comparison', ...
    'plot_amplitude_frequency_correlation', ...
    'plot_event_raster', ...
    'plot_biological_variability', ...
    'plot_roi_recruitment_curves', ...
    'plot_cumulative_events_over_time'
};

fprintf('Testing function accessibility:\n');
for i = 1:length(functions_to_test)
    func_name = functions_to_test{i};
    if exist(func_name, 'file')
        fprintf('  ✓ %s\n', func_name);
    else
        fprintf('  ✗ %s NOT FOUND\n', func_name);
    end
end

% Test 2: Verify results structure
fprintf('\nVerifying results structure:\n');
required_fields = {'wt_data', 'mut_data', 'comparison', 'file_results', ...
                  'processing_info', 'plot_handles'};
for i = 1:length(required_fields)
    if isfield(results, required_fields{i})
        fprintf('  ✓ results.%s\n', required_fields{i});
    else
        fprintf('  ✗ results.%s MISSING\n', required_fields{i});
    end
end

% Test 3: Verify frame rate
fprintf('\nFrame rate verification:\n');
if isfield(results.processing_info, 'frame_rate')
    fprintf('  Frame rate: %d Hz\n', results.processing_info.frame_rate);
    if results.processing_info.frame_rate == 100
        fprintf('  ✓ Correct for 10ms exposure data\n');
    else
        fprintf('  ! Warning: Expected 100 Hz for 3000-frame data\n');
    end
else
    fprintf('  ! Frame rate not set in processing_info\n');
end

fprintf('\n=== All Tests Complete ===\n');

%% ========================================================================
%% END OF EXAMPLE SCRIPT
%% ========================================================================
%% You can now:
%% 1. Run entire script with Run (F5)
%% 2. Run individual sections with Run Section (Ctrl+Enter)
%% 3. Copy specific sections for your analysis workflow
%% 4. Modify options to suit your needs
%% ========================================================================
