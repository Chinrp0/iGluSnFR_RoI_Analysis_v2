%% MINIMAL BASELINE TEST - Test individual components step by step
% This script tests each part of the baseline calculation separately
% to isolate any remaining issues

clear; clc; close all;

fprintf('=== Minimal Baseline Calculation Test ===\n');

% Test file path
test_file = 'D:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglusnfr4f_NGR\Spont\GPU_SNR_Processed\5_raw_mean\CP_Snfr4-NGR_Doc2b-R213W_Cs1-c1_spont-01_mean.csv';

% Add paths
addpath(genpath(pwd));

%% Step 1: Test CSV Loading
fprintf('\n1. Testing CSV loading...\n');
try
    loader = csv_loader_v2();
    [data, metadata] = loader.loadSingleFile(test_file);
    fprintf('  ✓ CSV loaded: %d frames × %d ROIs\n', size(data, 1), size(data, 2));
catch ME
    fprintf('  ✗ CSV loading failed: %s\n', ME.message);
    return;
end

%% Step 2: Test Configuration
fprintf('\n2. Testing configuration...\n');
try
    config = tracenorm_config();
    fprintf('  ✓ Config loaded: %.1fs window, %.1fσ threshold\n', ...
        config.rolling_window_sec, config.outlier_threshold_sigma);
catch ME
    fprintf('  ✗ Config loading failed: %s\n', ME.message);
    return;
end

%% Step 3: Test Baseline Detection
fprintf('\n3. Testing baseline detection...\n');
try
    [baseline, outlier_mask, baseline_stats] = baseline_detector(data, config);
    fprintf('  ✓ Baseline calculated: %.2f%% outliers detected\n', ...
        100 * baseline_stats.mean_outlier_fraction);
    fprintf('    Good ROIs: %.1f%%, Transport ROIs: %d\n', ...
        100 * baseline_stats.fraction_good_rois, baseline_stats.num_transport_rois);
catch ME
    fprintf('  ✗ Baseline detection failed: %s\n', ME.message);
    fprintf('    Error location: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
    return;
end

%% Step 4: Test dF/F Calculation
fprintf('\n4. Testing dF/F calculation...\n');
try
    [dfof_data, dfof_stats] = dfof_calculator(data, baseline, config);
    fprintf('  ✓ dF/F calculated: Mean SNR = %.2f\n', dfof_stats.mean_snr);
    fprintf('    Active ROIs: %d (%.1f%%), Total events: %d\n', ...
        dfof_stats.event_summary.rois_with_events, ...
        100 * dfof_stats.event_summary.fraction_active_rois, ...
        dfof_stats.event_summary.total_events);
catch ME
    fprintf('  ✗ dF/F calculation failed: %s\n', ME.message);
    fprintf('    Error location: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
    return;
end

%% Step 5: Check Data Structures
fprintf('\n5. Checking data structures...\n');
fprintf('  baseline_stats fields: %s\n', strjoin(fieldnames(baseline_stats), ', '));
fprintf('  dfof_stats fields: %s\n', strjoin(fieldnames(dfof_stats), ', '));

% Check for the problematic snr field
if isfield(dfof_stats, 'snr')
    fprintf('  ✓ SNR field exists in dfof_stats\n');
else
    fprintf('  ✗ SNR field missing from dfof_stats\n');
end

%% Step 6: Test Simple Plotting (No GUI plots)
fprintf('\n6. Testing basic plotting...\n');
try
    % Create a simple figure to test basic plotting
    fig = figure('Name', 'Basic Test Plot', 'Position', [100, 100, 800, 400]);
    
    % Select first 3 ROIs for testing
    sample_rois = 1:min(3, size(data, 2));
    time_vector = (1:size(data, 1)) / config.frame_rate;
    
    for i = 1:length(sample_rois)
        roi_idx = sample_rois(i);
        
        subplot(3, 1, i);
        plot(time_vector, data(:, roi_idx), 'k-', 'LineWidth', 0.8);
        hold on;
        plot(time_vector, baseline(:, roi_idx), 'r-', 'LineWidth', 1.2);
        
        % Mark outliers
        outliers = outlier_mask(:, roi_idx);
        if any(outliers)
            plot(time_vector(outliers), data(outliers, roi_idx), 'o', ...
                'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor', [1 0.5 0], 'MarkerSize', 4);
        end
        
        xlabel('Time (s)');
        ylabel('Fluorescence');
        title(sprintf('ROI %d', roi_idx));
        grid on;
    end
    
    sgtitle('Baseline Calculation Test - Sample ROIs');
    fprintf('  ✓ Basic plotting successful\n');
    
catch ME
    fprintf('  ✗ Basic plotting failed: %s\n', ME.message);
    fprintf('    Error location: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
end

%% Step 7: Summary
fprintf('\n=== Test Summary ===\n');
fprintf('Data size: %d frames × %d ROIs (%.1f MB)\n', ...
    size(data, 1), size(data, 2), numel(data) * 4 / (1024^2));
fprintf('Baseline quality: %.1f%% good ROIs\n', 100 * baseline_stats.fraction_good_rois);
fprintf('dF/F quality: %.1f%% good ROIs\n', 100 * dfof_stats.fraction_good_rois);
fprintf('Event detection: %d active ROIs with %d total events\n', ...
    dfof_stats.event_summary.rois_with_events, dfof_stats.event_summary.total_events);

if baseline_stats.num_transport_rois > 0
    fprintf('Transport artifacts: %d ROIs detected\n', baseline_stats.num_transport_rois);
end

fprintf('\n✅ All core functions working! Issue may be in full plotting pipeline.\n');
fprintf('Next: Try main pipeline with createPlots = false\n');