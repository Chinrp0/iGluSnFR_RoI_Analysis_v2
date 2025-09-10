%% TEST_BASELINE_SINGLE_FILE - Test iterative rolling median baseline calculation
% Quick script to test the new baseline calculation on a single CSV file
% 
% This script demonstrates:
% 1. Single file processing with baseline calculation
% 2. Vectorized processing of 1000+ ROIs 
% 3. Quality control visualization
% 4. Performance benchmarking

clear; clc; close all;

%% === Configuration ===
% UPDATE THIS PATH to point to your test CSV file
test_file = 'path/to/your/test_file.csv';  % CHANGE THIS

% Test options
test_options = struct();
test_options.verbose = true;           % Show detailed progress
test_options.createPlots = true;       % Generate baseline validation plots
test_options.savePlots = false;        % Set to true to save plots to disk

fprintf('=== Testing Baseline Calculation on Single File ===\n');

%% === Setup Pipeline ===
% Ensure all modules are accessible
if ~exist('csv_loader_v2', 'file')
    fprintf('Setting up pipeline...\n');
    setup_pipeline();
end

%% === Verify Test File ===
if ~exist(test_file, 'file')
    fprintf('ERROR: Test file not found: %s\n', test_file);
    fprintf('Please update the test_file variable with the path to your CSV file.\n');
    fprintf('\nExample CSV files in your data should have:\n');
    fprintf('  - 1200 frames (rows)\n');
    fprintf('  - 1000+ ROIs (columns)\n');
    fprintf('  - 40 Hz sampling rate\n');
    return;
end

%% === Test Single File Processing ===
fprintf('\nTesting baseline calculation...\n');

% Time the entire process
tic;
results = main_pipeline(test_file, test_options);
total_time = toc;

%% === Display Results ===
if ~isempty(results.fileResults)
    fileResult = results.fileResults{1};  % Get the single file result
    
    fprintf('\n=== Baseline Calculation Results ===\n');
    fprintf('File: %s\n', fileResult.filename);
    fprintf('Data: %d frames × %d ROIs (%.1f MB)\n', ...
        fileResult.numFrames, fileResult.numROIs, fileResult.dataSize_MB);
    
    % Baseline quality
    fprintf('\nBaseline Quality:\n');
    fprintf('  Good ROIs: %.1f%% (%d/%d)\n', ...
        100 * fileResult.quality.fraction_good_baseline, ...
        round(fileResult.quality.fraction_good_baseline * fileResult.numROIs), ...
        fileResult.numROIs);
    fprintf('  Transport ROIs: %d\n', fileResult.quality.transport_rois);
    fprintf('  Mean outliers: %.2f%% per ROI\n', ...
        100 * fileResult.baseline_stats.mean_outlier_fraction);
    
    % dF/F quality  
    fprintf('\ndF/F Quality:\n');
    fprintf('  Good ROIs: %.1f%% (%d/%d)\n', ...
        100 * fileResult.quality.fraction_good_dfof, ...
        round(fileResult.quality.fraction_good_dfof * fileResult.numROIs), ...
        fileResult.numROIs);
    fprintf('  Active ROIs: %d (%.1f%%)\n', ...
        fileResult.quality.active_rois, ...
        100 * fileResult.quality.active_rois / fileResult.numROIs);
    fprintf('  Mean SNR: %.2f\n', fileResult.quality.mean_snr);
    fprintf('  Total events detected: %d\n', fileResult.dfof_stats.event_summary.total_events);
    
    % Performance
    fprintf('\nPerformance:\n');
    fprintf('  Baseline calculation: %.3f s\n', fileResult.timing.baseline_time);
    fprintf('  dF/F calculation: %.3f s\n', fileResult.timing.dfof_time);
    fprintf('  Visualization: %.3f s\n', fileResult.timing.plot_time);
    fprintf('  Total processing: %.3f s\n', total_time);
    fprintf('  Processing rate: %.1f ROIs/s\n', fileResult.numROIs / total_time);
    
    % Vectorization efficiency check
    if fileResult.numROIs > 1000
        fprintf('  ✓ Successfully processed %d ROIs in vectorized mode\n', fileResult.numROIs);
    end
    
    %% === Access Calculated Data ===
    fprintf('\n=== Data Access Examples ===\n');
    fprintf('Access baseline data:\n');
    fprintf('  baseline = results.fileResults{1}.baseline;  %% [1200 x %d]\n', fileResult.numROIs);
    fprintf('Access dF/F data:\n');
    fprintf('  dfof = results.fileResults{1}.dfof_data;     %% [1200 x %d]\n', fileResult.numROIs);
    fprintf('Access outlier mask:\n');
    fprintf('  outliers = results.fileResults{1}.outlier_mask;  %% [1200 x %d logical]\n', fileResult.numROIs);
    
    %% === Configuration Info ===
    config = tracenorm_config();
    fprintf('\n=== Baseline Algorithm Configuration ===\n');
    fprintf('Method: %s\n', config.baseline_method);
    fprintf('Rolling window: %.1f s (%d frames at %d Hz)\n', ...
        config.rolling_window_sec, config.rolling_window_frames, config.frame_rate);
    fprintf('Outlier threshold: %.1f standard deviations\n', config.outlier_threshold_sigma);
    fprintf('Max iterations: %d\n', config.max_iterations);
    
    %% === Quality Recommendations ===
    fprintf('\n=== Quality Assessment ===\n');
    
    if fileResult.quality.fraction_good_baseline > 0.8
        fprintf('✓ Baseline quality: GOOD (>80%% valid ROIs)\n');
    elseif fileResult.quality.fraction_good_baseline > 0.6
        fprintf('⚠ Baseline quality: MODERATE (60-80%% valid ROIs)\n');
    else
        fprintf('✗ Baseline quality: POOR (<60%% valid ROIs)\n');
        fprintf('  Consider adjusting outlier threshold or checking data quality\n');
    end
    
    if fileResult.quality.mean_snr > 3
        fprintf('✓ Signal quality: GOOD (SNR > 3)\n');
    elseif fileResult.quality.mean_snr > 2
        fprintf('⚠ Signal quality: MODERATE (SNR 2-3)\n');
    else
        fprintf('✗ Signal quality: POOR (SNR < 2)\n');
        fprintf('  Low signal-to-noise ratio detected\n');
    end
    
    if fileResult.quality.active_rois > 0
        fprintf('✓ Event detection: %d ROIs show activity\n', fileResult.quality.active_rois);
    else
        fprintf('⚠ Event detection: No clear events detected\n');
        fprintf('  Check if threshold settings are appropriate\n');
    end
    
else
    fprintf('ERROR: Failed to process file\n');
    if ~isempty(results.errors.failedFiles)
        fprintf('Error: %s\n', results.errors.errorMessages{1});
    end
end

%% === Next Steps ===
fprintf('\n=== Next Steps ===\n');
fprintf('1. Review the baseline validation plots\n');
fprintf('2. Check sample traces for baseline quality\n');
fprintf('3. Adjust parameters in tracenorm_config.m if needed\n');
fprintf('4. Test on additional files using: results = main_pipeline(folder_path)\n');
fprintf('5. Use baseline data for downstream analysis\n');

fprintf('\n=== Test Complete ===\n');