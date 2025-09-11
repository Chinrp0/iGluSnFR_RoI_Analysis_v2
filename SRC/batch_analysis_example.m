%% BATCH_CONDITION_ANALYSIS_EXAMPLE
% Example script showing how to use the batch condition analysis functions
% to compare WT vs R213W glutamate release events
%
% This script demonstrates the complete workflow for:
% 1. Processing multiple CSV files using existing pipeline
% 2. Automatically classifying files as WT or R213W
% 3. Aggregating event data by condition  
% 4. Statistical comparison between conditions
% 5. Creating comprehensive comparison visualizations
%
% REQUIREMENTS:
% - Existing pipeline functions (main_pipeline, corrected_schmitt_detector, etc.)
% - CSV files with 'WT' or 'R213W' in filenames for auto-classification
% - MATLAB Statistics and Machine Learning Toolbox for statistical tests

%% === SETUP ===
clear; close all; clc;

% Verify pipeline is available
if ~exist('main_pipeline', 'file')
    error('main_pipeline function not found. Please run setup_pipeline() first.');
end

fprintf('=== Batch Condition Analysis Example ===\n');
fprintf('Comparing WT vs R213W glutamate release events\n\n');

%% === CONFIGURATION ===

% Set your data folder path here
data_folder = '/path/to/your/csv/files';  % UPDATE THIS PATH

% If testing with a specific folder, uncomment and modify:
% data_folder = '/Users/yourname/Data/GlutamateRelease';

% Analysis options
options = struct();
options.verbose = true;                % Print detailed progress
options.createPlots = true;            % Generate comparison plots  
options.savePlots = true;              % Save plots to files
options.useParallel = false;           % Set to true if Parallel Computing Toolbox available
options.continueOnError = true;        % Continue if individual files fail

fprintf('Data folder: %s\n', data_folder);
fprintf('Options: verbose=%s, plots=%s, save=%s\n', ...
    mat2str(options.verbose), mat2str(options.createPlots), mat2str(options.savePlots));

%% === RUN BATCH ANALYSIS ===

try
    fprintf('\n=== Starting Batch Analysis ===\n');
    tic;
    
    % Run the complete batch condition analysis
    results = batch_condition_analysis(data_folder, options);
    
    analysis_time = toc;
    fprintf('\n=== Analysis Complete ===\n');
    fprintf('Total analysis time: %.2f seconds\n', analysis_time);
    
catch ME
    fprintf('\nERROR during batch analysis:\n');
    fprintf('  %s\n', ME.message);
    if ~isempty(ME.stack)
        fprintf('  In function: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
    end
    return;
end

%% === EXPLORE RESULTS ===

fprintf('\n=== Results Summary ===\n');

% Display basic information
fprintf('Files processed:\n');
fprintf('  WT: %d files\n', length(results.processing_info.wt_files));
fprintf('  R213W: %d files\n', length(results.processing_info.mut_files));

% Display key findings
wt = results.wt_data;
mut = results.mut_data;
comp = results.comparison;

fprintf('\nKey Findings:\n');
fprintf('1. Active ROI Fraction:\n');
fprintf('   WT: %.1f%% (%d/%d ROIs)\n', ...
    100 * wt.stats.frequency.active_rois.fraction, ...
    wt.stats.frequency.active_rois.count, wt.total_rois);
fprintf('   R213W: %.1f%% (%d/%d ROIs)\n', ...
    100 * mut.stats.frequency.active_rois.fraction, ...
    mut.stats.frequency.active_rois.count, mut.total_rois);

if isfield(comp, 'activity_pvalue') && ~isnan(comp.activity_pvalue)
    fprintf('   Statistical difference: p = %.4f\n', comp.activity_pvalue);
end

fprintf('\n2. Event Frequency (Active ROIs):\n');
fprintf('   WT: %.4f ± %.4f Hz\n', ...
    wt.stats.frequency.active_rois.mean, wt.stats.frequency.active_rois.std);
fprintf('   R213W: %.4f ± %.4f Hz\n', ...
    mut.stats.frequency.active_rois.mean, mut.stats.frequency.active_rois.std);

if isfield(comp, 'frequency_pvalue') && ~isnan(comp.frequency_pvalue)
    fprintf('   Statistical difference: p = %.4f\n', comp.frequency_pvalue);
end

fprintf('\n3. Event Amplitude (All Events):\n');
fprintf('   WT: %.4f ± %.4f dF/F (%d events)\n', ...
    wt.stats.amplitude.mean, wt.stats.amplitude.std, wt.stats.amplitude.count);
fprintf('   R213W: %.4f ± %.4f dF/F (%d events)\n', ...
    mut.stats.amplitude.mean, mut.stats.amplitude.std, mut.stats.amplitude.count);

if isfield(comp, 'amplitude_pvalue') && ~isnan(comp.amplitude_pvalue)
    fprintf('   Statistical difference: p = %.4f\n', comp.amplitude_pvalue);
end

%% === ADDITIONAL ANALYSIS EXAMPLES ===

fprintf('\n=== Additional Analysis Options ===\n');

% Example 1: Export summary table to CSV
if isfield(results, 'comparison')
    try
        timestamp = datestr(now, 'yyyymmdd_HHMMSS');
        summary_file = sprintf('condition_comparison_summary_%s.csv', timestamp);
        summary_table = export_condition_summary(results, summary_file);
        
        if ~isempty(summary_table)
            fprintf('1. Summary table exported to: %s\n', summary_file);
            
            % Display table preview
            fprintf('   Table preview:\n');
            disp(summary_table(1:5, :));  % Show first 5 rows
        end
    catch
        fprintf('1. Summary table export failed\n');
    end
end

% Example 2: Access individual file results
fprintf('2. Individual file results available:\n');
fprintf('   results.file_results{i} contains per-file event statistics\n');

% Show example of accessing per-file data
if ~isempty(results.file_results)
    valid_files = find(~cellfun(@isempty, results.file_results));
    if ~isempty(valid_files)
        example_file = results.file_results{valid_files(1)};
        fprintf('   Example (first file): %d ROIs, %d events\n', ...
            example_file.numROIs, example_file.event_stats.total_events);
    end
end

% Example 3: Plot handles for further customization  
if isfield(results, 'plot_handles') && ~isempty(results.plot_handles)
    plot_names = fieldnames(results.plot_handles);
    fprintf('3. Generated plots (%d):\n', length(plot_names));
    for i = 1:length(plot_names)
        fprintf('   %s\n', plot_names{i});
    end
    fprintf('   Access with: results.plot_handles.frequency_comparison\n');
end

%% === SAVE WORKSPACE (OPTIONAL) ===

% Uncomment to save complete results for later analysis
% timestamp = datestr(now, 'yyyymmdd_HHMMSS');
% results_file = sprintf('batch_analysis_results_%s.mat', timestamp);
% save(results_file, 'results', 'options');
% fprintf('\n4. Complete results saved to: %s\n', results_file);

%% === CUSTOM ANALYSIS EXAMPLES ===

fprintf('\n=== Custom Analysis Examples ===\n');

% Example: Find most active ROI in each condition
if ~isempty(wt.all_roi_frequencies) && ~isempty(mut.all_roi_frequencies)
    wt_max_freq = max(wt.all_roi_frequencies);
    mut_max_freq = max(mut.all_roi_frequencies);
    fprintf('Most active ROI frequency:\n');
    fprintf('  WT: %.4f Hz\n', wt_max_freq);
    fprintf('  R213W: %.4f Hz\n', mut_max_freq);
end

% Example: Count high-amplitude events (top 10%)
if ~isempty(wt.all_individual_events) && ~isempty(mut.all_individual_events)
    wt_top10_thresh = prctile(wt.all_individual_events, 90);
    mut_top10_thresh = prctile(mut.all_individual_events, 90);
    
    wt_high_amp_count = sum(wt.all_individual_events > wt_top10_thresh);
    mut_high_amp_count = sum(mut.all_individual_events > mut_top10_thresh);
    
    fprintf('High-amplitude events (top 10%%):\n');
    fprintf('  WT: %d events > %.4f dF/F\n', wt_high_amp_count, wt_top10_thresh);
    fprintf('  R213W: %d events > %.4f dF/F\n', mut_high_amp_count, mut_top10_thresh);
end

fprintf('\n=== Analysis Complete ===\n');
fprintf('Results stored in ''results'' variable\n');
fprintf('Plots displayed and optionally saved to files\n');

%% === QUICK START INSTRUCTIONS ===
% 
% TO USE THIS SCRIPT:
% 1. Update the 'data_folder' variable to point to your CSV files
% 2. Ensure CSV filenames contain 'WT' or 'R213W' for auto-classification  
% 3. Run the script
% 4. Examine the generated plots and printed statistics
% 5. Access detailed results in the 'results' variable
%
% EXPECTED CSV FILENAME PATTERNS:
% - 'CP_Snfr4-NGR_Doc2b-WT_Cs1-c1_spont-01.csv' (classified as WT)
% - 'CP_Snfr4-NGR_Doc2b-R213W_Cs1-c1_spont-01.csv' (classified as R213W) 
%
% OUTPUTS:
% - results.wt_data: Aggregated WT condition data
% - results.mut_data: Aggregated R213W condition data  
% - results.comparison: Statistical comparison results
% - results.plot_handles: Figure handles for generated plots
% - Condition comparison plots (4 comprehensive figures)
% - Optional: CSV summary table and saved plot files