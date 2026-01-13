%% CONFIG PARAMETER VERIFICATION
% Verify that your config changes are actually being used by the functions
% Shows where parameters are stored and how to access them

clear; clc;

fprintf('=== Config Parameter Verification ===\n');

%% === Step 1: Check Current Config ===
fprintf('\n1. Current configuration settings:\n');

config = tracenorm_config();

% Display key parameters
fprintf('   Frame rate: %d Hz\n', config.frame_rate);
fprintf('   Rolling window: %.2f s (%d frames)\n', ...
    config.rolling_window_sec, config.rolling_window_frames);
fprintf('   Outlier threshold: %.1f σ\n', config.outlier_threshold_sigma);

fprintf('\n   Schmitt trigger settings:\n');
fprintf('   Upper threshold: %.1f σ\n', config.event_detection.upper_threshold_sigma);
fprintf('   Lower threshold: %.1f σ\n', config.event_detection.lower_threshold_sigma);
fprintf('   Decay extension: %d frames\n', config.event_detection.decay_extension_frames);
fprintf('   Rise extension: %d frames\n', config.event_detection.rise_extension_frames);

%% === Step 2: Test Parameter Changes ===
fprintf('\n2. Testing parameter modifications:\n');

% Backup original config
original_config = config;

% Modify some parameters
fprintf('   Modifying parameters for testing...\n');
config.event_detection.upper_threshold_sigma = 4.0;  % More conservative
config.event_detection.lower_threshold_sigma = 2.0;   % Higher lower threshold
config.event_detection.decay_extension_frames = 5;    % Longer extension

fprintf('   NEW Upper threshold: %.1f σ\n', config.event_detection.upper_threshold_sigma);
fprintf('   NEW Lower threshold: %.1f σ\n', config.event_detection.lower_threshold_sigma);
fprintf('   NEW Decay extension: %d frames\n', config.event_detection.decay_extension_frames);

%% === Step 3: Test with Modified Config ===
fprintf('\n3. Testing event detection with modified config:\n');

% Load sample data
addpath(genpath(pwd));
test_file = 'D:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglusnfr4f_NGR\Spont\GPU_SNR_Processed\5_raw_mean\CP_Snfr4-NGR_Doc2b-R213W_Cs1-c1_spont-01_mean.csv';

loader = csv_loader_v2();
[data, metadata] = loader.loadSingleFile(test_file);

[baseline, outlier_mask, baseline_stats] = baseline_detector(data, config);
[dfof_data, dfof_stats] = dfof_calculator(data, baseline, config);

% Test with ORIGINAL config
fprintf('   Testing with ORIGINAL config (3.0σ/1.5σ)...\n');
[event_mask_orig, event_stats_orig] = optimal_schmitt_event_detector(dfof_data, original_config, baseline_stats);

% Test with MODIFIED config
fprintf('   Testing with MODIFIED config (4.0σ/2.0σ)...\n');
[event_mask_mod, event_stats_mod] = optimal_schmitt_event_detector(dfof_data, config, baseline_stats);

%% === Step 4: Compare Results ===
fprintf('\n4. Comparing results to verify parameter usage:\n');
fprintf('Configuration          Total Events    Active ROIs    Mean Threshold\n');
fprintf('--------------------   ------------    -----------    --------------\n');
fprintf('Original (3.0σ/1.5σ)   %8d        %8d       %.5f\n', ...
    event_stats_orig.total_events, event_stats_orig.rois_with_events, ...
    mean(event_stats_orig.upper_thresholds, 'omitnan'));

fprintf('Modified (4.0σ/2.0σ)   %8d        %8d       %.5f\n', ...
    event_stats_mod.total_events, event_stats_mod.rois_with_events, ...
    mean(event_stats_mod.upper_thresholds, 'omitnan'));

% Calculate differences
event_ratio = event_stats_mod.total_events / event_stats_orig.total_events;
threshold_ratio = mean(event_stats_mod.upper_thresholds, 'omitnan') / mean(event_stats_orig.upper_thresholds, 'omitnan');

fprintf('\nParameter effect verification:\n');
fprintf('   Threshold ratio (modified/original): %.2f (should be ~1.33 for 4σ/3σ)\n', threshold_ratio);
fprintf('   Event ratio (modified/original): %.2f (should be <1.0 for higher thresholds)\n', event_ratio);

if abs(threshold_ratio - 4.0/3.0) < 0.1
    fprintf('   ✅ CONFIRMED: Threshold parameters are being used correctly\n');
else
    fprintf('   ❌ WARNING: Threshold parameters may not be applied correctly\n');
end

if event_ratio < 1.0
    fprintf('   ✅ CONFIRMED: Higher thresholds reduce event detection as expected\n');
else
    fprintf('   ❌ WARNING: Higher thresholds should reduce events\n');
end

%% === Step 5: Check Where Parameters Are Stored ===
fprintf('\n5. Parameter storage in results:\n');

% Run full pipeline to see where parameters end up
results = main_pipeline(test_file, struct('verbose', false, 'createPlots', false));

if ~isempty(results.fileResults)
    fileResult = results.fileResults{1};
    
    fprintf('   Parameters stored in results structure:\n');
    
    % Check event_stats for threshold info
    if isfield(fileResult, 'event_stats')
        fprintf('   ✓ fileResult.event_stats.upper_thresholds: Available [1 x %d]\n', ...
            length(fileResult.event_stats.upper_thresholds));
        fprintf('   ✓ fileResult.event_stats.lower_thresholds: Available [1 x %d]\n', ...
            length(fileResult.event_stats.lower_thresholds));
        
        % Show actual values
        sample_roi = 100;
        if length(fileResult.event_stats.upper_thresholds) >= sample_roi
            fprintf('     Example ROI %d: Upper=%.5f, Lower=%.5f\n', sample_roi, ...
                fileResult.event_stats.upper_thresholds(sample_roi), ...
                fileResult.event_stats.lower_thresholds(sample_roi));
        end
    end
    
    % Check for noise metrics (if using optimal detector)
    if isfield(fileResult.event_stats, 'noise_metrics')
        fprintf('   ✓ fileResult.event_stats.noise_metrics: Enhanced metrics available\n');
        fprintf('     - noise_std: [1 x %d] noise levels per ROI\n', ...
            length(fileResult.event_stats.noise_metrics.noise_std));
        fprintf('     - noise_method: {1 x %d} methods used per ROI\n', ...
            length(fileResult.event_stats.noise_metrics.noise_method));
        fprintf('     - signal_quality: [1 x %d] SNR values per ROI\n', ...
            length(fileResult.event_stats.noise_metrics.signal_quality));
        
        % Show method distribution
        methods_used = fileResult.event_stats.noise_metrics.methods_used;
        fprintf('     - Methods used: %s\n', strjoin(methods_used, ', '));
    else
        fprintf('   - fileResult.event_stats.noise_metrics: Not available (using basic detector)\n');
    end
    
    % Check if config is stored anywhere
    fprintf('\n   Configuration preservation:\n');
    if isfield(fileResult, 'config')
        fprintf('   ✓ fileResult.config: Full config preserved\n');
    else
        fprintf('   - fileResult.config: Not stored (parameters only in computed results)\n');
    end
end

%% === Step 6: Practical Parameter Access ===
fprintf('\n6. How to access and verify your parameters:\n');

fprintf('\n   A. From your config file:\n');
fprintf('      config = tracenorm_config();\n');
fprintf('      upper_sigma = config.event_detection.upper_threshold_sigma;  %% %.1f\n', ...
    config.event_detection.upper_threshold_sigma);
fprintf('      lower_sigma = config.event_detection.lower_threshold_sigma;  %% %.1f\n', ...
    config.event_detection.lower_threshold_sigma);

fprintf('\n   B. From pipeline results:\n');
fprintf('      results = main_pipeline(your_file);\n');
fprintf('      thresholds = results.fileResults{1}.event_stats.upper_thresholds;\n');
fprintf('      mean_threshold = mean(thresholds);  %% Actual thresholds used\n');

fprintf('\n   C. Verify parameter application:\n');
fprintf('      noise_levels = results.fileResults{1}.event_stats.noise_metrics.noise_std;\n');
fprintf('      expected_threshold = config.event_detection.upper_threshold_sigma * noise_levels;\n');
fprintf('      actual_threshold = results.fileResults{1}.event_stats.upper_thresholds;\n');
fprintf('      parameter_match = isequal(expected_threshold, actual_threshold);\n');

%% === Step 7: Config Modification Examples ===
fprintf('\n7. Examples of config modifications:\n');

fprintf('\n   To make detection MORE sensitive (lower thresholds):\n');
fprintf('      config.event_detection.upper_threshold_sigma = 2.5;  %% Default: 3.0\n');
fprintf('      config.event_detection.lower_threshold_sigma = 1.0;  %% Default: 1.5\n');

fprintf('\n   To make detection LESS sensitive (higher thresholds):\n');
fprintf('      config.event_detection.upper_threshold_sigma = 4.0;  %% Default: 3.0\n');
fprintf('      config.event_detection.lower_threshold_sigma = 2.5;  %% Default: 1.5\n');

fprintf('\n   To change event extension (kinetics):\n');
fprintf('      config.event_detection.decay_extension_frames = 5;   %% Default: 3\n');
fprintf('      config.event_detection.rise_extension_frames = 2;    %% Default: 1\n');

fprintf('\n   To save config in results (add to main_pipeline.m):\n');
fprintf('      fileResult.config_used = config;  %% Store config for reference\n');

%% === Step 8: Parameter Verification Summary ===
fprintf('\n=== PARAMETER VERIFICATION SUMMARY ===\n');

if abs(threshold_ratio - 4.0/3.0) < 0.1 && event_ratio < 1.0
    fprintf('✅ CONFIG PARAMETERS ARE WORKING CORRECTLY\n');
    fprintf('   - Threshold changes are applied\n');
    fprintf('   - Event detection responds as expected\n');
    fprintf('   - Parameters are accessible in results\n');
else
    fprintf('⚠️  POTENTIAL CONFIG ISSUES DETECTED\n');
    fprintf('   - Check if modified config is passed to functions\n');
    fprintf('   - Verify function calls use your config, not default\n');
    fprintf('   - Ensure no caching of old parameters\n');
end

fprintf('\nWhere your parameters end up:\n');
fprintf('   1. Applied during detection: event_stats.upper_thresholds, lower_thresholds\n');
fprintf('   2. Noise estimation: event_stats.noise_metrics (if using optimal detector)\n');
fprintf('   3. Config not stored by default (only computed results)\n');

fprintf('\nTo permanently store config in results, add this line to main_pipeline.m:\n');
fprintf('   fileResult.config_used = config;\n');

fprintf('\n✅ Parameter verification complete!\n');