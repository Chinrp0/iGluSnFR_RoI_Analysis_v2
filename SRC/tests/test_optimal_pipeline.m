%% TEST OPTIMAL PIPELINE INTEGRATION
% Test the optimal event detector integrated into your main pipeline
% Compares performance and results with current method

clear; clc; close all;

fprintf('=== Testing Optimal Pipeline Integration ===\n');

%% === Setup ===
addpath(genpath(pwd));
test_file = 'D:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglusnfr4f_NGR\Spont\GPU_SNR_Processed\5_raw_mean\CP_Snfr4-NGR_Doc2b-R213W_Cs1-c1_spont-01_mean.csv';

% Test options
test_options = struct();
test_options.verbose = true;
test_options.createPlots = true;
test_options.savePlots = false;

%% === Test 1: Current Pipeline ===
fprintf('\n1. Testing CURRENT pipeline...\n');
tic;
results_current = main_pipeline(test_file, test_options);
time_current = toc;

fprintf('   Current pipeline: %.3f seconds\n', time_current);
if ~isempty(results_current.fileResults)
    current_events = results_current.fileResults{1}.event_stats.total_events;
    current_active = results_current.fileResults{1}.event_stats.rois_with_events;
    fprintf('   Current results: %d events, %d active ROIs\n', current_events, current_active);
end

%% === Test 2: Optimal Pipeline (Manual Integration) ===
fprintf('\n2. Testing OPTIMAL pipeline with manual integration...\n');

% Load and process data manually with optimal detector
loader = csv_loader_v2();
[data, metadata] = loader.loadSingleFile(test_file);
config = tracenorm_config();
config.verbose = false;  % Reduce verbosity for cleaner output

% Setup parallel pool
parpool_obj = gcp('nocreate');
if isempty(parpool_obj)
    fprintf('   Starting parallel pool...\n');
    parpool('local', min(4, feature('numcores')));
end

tic;
% Standard processing
[baseline, outlier_mask, baseline_stats] = baseline_detector(data, config);
[dfof_data, dfof_stats] = dfof_calculator(data, baseline, config);

% OPTIMAL event detection
[event_mask_optimal, event_stats_optimal] = optimal_schmitt_event_detector(dfof_data, config, baseline_stats);

% Use fixed plotting with optimal results
plot_handles_optimal = baseline_plotter(data, baseline, dfof_data, outlier_mask, ...
    baseline_stats, dfof_stats, metadata, config, event_stats_optimal);

time_optimal = toc;

fprintf('   Optimal pipeline: %.3f seconds\n', time_optimal);
fprintf('   Optimal results: %d events, %d active ROIs\n', ...
    event_stats_optimal.total_events, event_stats_optimal.rois_with_events);

%% === Test 3: Performance Comparison ===
fprintf('\n=== PERFORMANCE COMPARISON ===\n');
fprintf('Metric                    Current      Optimal      Improvement\n');
fprintf('------------------------  -----------  -----------  -----------\n');
fprintf('Processing time           %.3f s      %.3f s      %.1fx %s\n', ...
    time_current, time_optimal, time_current/time_optimal, ...
    ternary(time_optimal < time_current, 'faster', 'slower'));

if ~isempty(results_current.fileResults)
    current_file = results_current.fileResults{1};
    
    fprintf('Total events              %6d       %6d       %.1fx more\n', ...
        current_events, event_stats_optimal.total_events, ...
        event_stats_optimal.total_events / current_events);
    
    fprintf('Active ROIs               %6d       %6d       %.1fx more\n', ...
        current_active, event_stats_optimal.rois_with_events, ...
        event_stats_optimal.rois_with_events / current_active);
    
    % Compare noise estimation
    current_thresholds = current_file.event_stats.upper_thresholds;
    optimal_thresholds = event_stats_optimal.noise_metrics.upper_thresholds;
    
    threshold_improvement = mean(current_thresholds ./ optimal_thresholds, 'omitnan');
    fprintf('Sensitivity (1/threshold) %.6f       %.6f       %.1fx better\n', ...
        mean(current_thresholds, 'omitnan'), mean(optimal_thresholds, 'omitnan'), threshold_improvement);
end

%% === Test 4: Noise Estimation Comparison ===
fprintf('\n=== NOISE ESTIMATION COMPARISON ===\n');

if ~isempty(results_current.fileResults)
    current_file = results_current.fileResults{1};
    
    % Compare noise estimation methods
    current_noise = current_file.event_stats.upper_thresholds / 3.0;  % Back-calculate noise
    optimal_noise = event_stats_optimal.noise_metrics.noise_std;
    
    fprintf('Current method (entire trace):\n');
    fprintf('  Mean noise: %.5f ± %.5f dF/F\n', mean(current_noise, 'omitnan'), std(current_noise, 'omitnan'));
    
    fprintf('Optimal method (adaptive):\n');
    fprintf('  Mean noise: %.5f ± %.5f dF/F\n', mean(optimal_noise, 'omitnan'), std(optimal_noise, 'omitnan'));
    
    % Method breakdown
    methods_used = event_stats_optimal.noise_metrics.methods_used;
    fprintf('  Methods used: %s\n', strjoin(methods_used, ', '));
    
    for method = methods_used
        method_rois = strcmp(event_stats_optimal.noise_metrics.noise_method, method);
        method_count = sum(method_rois);
        method_percentage = 100 * method_count / length(optimal_noise);
        fprintf('    %s: %d ROIs (%.1f%%)\n', method{1}, method_count, method_percentage);
    end
    
    % Overall sensitivity improvement
    sensitivity_improvement = mean(current_noise ./ optimal_noise, 'omitnan');
    fprintf('  Overall sensitivity improvement: %.1fx\n', sensitivity_improvement);
end

%% === Test 5: ROI-Level Comparison ===
fprintf('\n=== ROI-LEVEL COMPARISON (Sample) ===\n');
fprintf('ROI   Current Events   Optimal Events   Current Noise   Optimal Noise   Method\n');
fprintf('---   --------------   --------------   -------------   -------------   --------\n');

sample_rois = [1, 50, 100, 200, 500, 1000];  % Sample ROIs
sample_rois = sample_rois(sample_rois <= size(dfof_data, 2));

for roi = sample_rois
    if ~isempty(results_current.fileResults)
        current_events_roi = current_file.event_stats.events_per_roi(roi);
        current_noise_roi = current_file.event_stats.upper_thresholds(roi) / 3.0;
    else
        current_events_roi = 0;
        current_noise_roi = 0;
    end
    
    optimal_events_roi = event_stats_optimal.events_per_roi(roi);
    optimal_noise_roi = event_stats_optimal.noise_metrics.noise_std(roi);
    method_used = event_stats_optimal.noise_metrics.noise_method{roi};
    
    fprintf('%3d   %8d         %8d         %.6f      %.6f      %s\n', ...
        roi, current_events_roi, optimal_events_roi, current_noise_roi, optimal_noise_roi, method_used);
end

%% === Test 6: Integration Instructions ===
fprintf('\n=== INTEGRATION INSTRUCTIONS ===\n');
fprintf('To integrate optimal detection into your main pipeline:\n\n');

fprintf('1. REPLACE the schmitt_event_detector call in main_pipeline.m:\n');
fprintf('   OLD: [event_mask, event_stats] = schmitt_event_detector(dfof_data, config);\n');
fprintf('   NEW: [event_mask, event_stats] = optimal_schmitt_event_detector(dfof_data, config, baseline_stats);\n\n');

fprintf('2. UPDATE your tracenorm_config.m to include optimal settings:\n');
fprintf('   config.optimal_detection = true;\n');
fprintf('   config.use_parallel_noise_estimation = true;\n');
fprintf('   config.enhanced_event_characterization = true;\n\n');

fprintf('3. ENSURE parallel pool is available:\n');
fprintf('   if isempty(gcp(''nocreate''))\n');
fprintf('       parpool(''local'', min(4, feature(''numcores'')));\n');
fprintf('   end\n\n');

fprintf('4. ACCESS enhanced noise metrics:\n');
fprintf('   noise_metrics = results.fileResults{1}.event_stats.noise_metrics;\n');
fprintf('   noise_levels = noise_metrics.noise_std;\n');
fprintf('   methods_used = noise_metrics.noise_method;\n');
fprintf('   signal_quality = noise_metrics.signal_quality;\n\n');

%% === Test 7: Memory and Resource Usage ===
fprintf('=== RESOURCE USAGE ANALYSIS ===\n');

% Check parallel pool status
pool = gcp('nocreate');
if ~isempty(pool)
    fprintf('Parallel pool: %d workers active\n', pool.NumWorkers);
else
    fprintf('Parallel pool: Not active\n');
end

% Memory usage estimate
data_size_mb = numel(dfof_data) * 8 / (1024^2);  % 8 bytes per double
fprintf('Data size: %.1f MB\n', data_size_mb);

% Toolbox usage
toolboxes = {'Parallel Computing Toolbox', 'Signal Processing Toolbox', 'Statistics and Machine Learning Toolbox'};
licenses = {'Distrib_Computing_Toolbox', 'Signal_Toolbox', 'Statistics_Toolbox'};

fprintf('Required toolboxes:\n');
for i = 1:length(toolboxes)
    if license('test', licenses{i})
        fprintf('  ✓ %s: Available\n', toolboxes{i});
    else
        fprintf('  ✗ %s: Not available\n', toolboxes{i});
    end
end

%% === Summary and Recommendations ===
fprintf('\n=== SUMMARY AND RECOMMENDATIONS ===\n');

if event_stats_optimal.total_events > current_events
    improvement_factor = event_stats_optimal.total_events / current_events;
    fprintf('✅ OPTIMAL DETECTION SUCCESSFUL!\n');
    fprintf('   %.1fx more events detected\n', improvement_factor);
    fprintf('   %.1fx more active ROIs\n', event_stats_optimal.rois_with_events / current_active);
    fprintf('   Enhanced noise characterization with multiple methods\n');
    fprintf('   Processing time: %.3f s (%.1fx %s than current)\n', ...
        time_optimal, abs(time_optimal/time_current), ...
        ternary(time_optimal < time_current, 'faster', 'slower'));
    
    fprintf('\nRECOMMENDATION: ⭐ Integrate optimal detection into main pipeline\n');
    
    if time_optimal > time_current * 2
        fprintf('   Consider: Reduce parallel workers if speed is critical\n');
    end
else
    fprintf('⚠️  INVESTIGATION NEEDED\n');
    fprintf('   Optimal detection found fewer events than current method\n');
    fprintf('   This may indicate over-conservative thresholds\n');
    fprintf('   Recommend: Review noise estimation parameters\n');
end

fprintf('\n✅ Integration testing complete!\n');

%% === Helper Function ===
function result = ternary(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end