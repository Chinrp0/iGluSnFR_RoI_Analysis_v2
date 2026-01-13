%% TEST_BASELINE_NO_PLOTS - Test pipeline without plotting
% This will isolate whether the issue is in core processing or visualization

clear; clc; close all;

%% === Configuration ===
test_file = "D:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglusnfr4f_NGR\Spont\GPU_SNR_Processed\5_raw_mean\CP_Snfr4-NGR_Doc2b-R213W_Cs1-c1_spont-01_mean.csv";

% Test options - DISABLE PLOTTING
test_options = struct();
test_options.verbose = true;
test_options.createPlots = false;    % DISABLE PLOTS TO TEST CORE FUNCTIONALITY
test_options.savePlots = false;

fprintf('=== Testing Core Pipeline Without Plots ===\n');

%% === Setup ===
if ~exist('csv_loader_v2', 'file')
    setup_pipeline();
end

%% === Test Core Processing ===
fprintf('\nTesting core pipeline...\n');

tic;
try
    results = main_pipeline(test_file, test_options);
    total_time = toc;
    
    fprintf('\n✅ SUCCESS: Core pipeline completed without plotting\n');
    
    % Display results
    if ~isempty(results.fileResults)
        fileResult = results.fileResults{1};
        
        fprintf('\n=== Core Processing Results ===\n');
        fprintf('File: %s\n', fileResult.filename);
        fprintf('Data: %d frames × %d ROIs (%.1f MB)\n', ...
            fileResult.numFrames, fileResult.numROIs, fileResult.dataSize_MB);
        
        fprintf('\nBaseline Quality: %.1f%% good ROIs\n', ...
            100 * fileResult.quality.fraction_good_baseline);
        fprintf('dF/F Quality: %.1f%% good ROIs\n', ...
            100 * fileResult.quality.fraction_good_dfof);
        fprintf('Mean SNR: %.2f\n', fileResult.quality.mean_snr);
        fprintf('Active ROIs: %d\n', fileResult.quality.active_rois);
        
        fprintf('\nProcessing Time: %.3f s\n', total_time);
        fprintf('Processing Rate: %.1f ROIs/s\n', fileResult.numROIs / total_time);
        
        fprintf('\n✅ All core functionality confirmed working!\n');
        fprintf('✅ Data structures are properly created\n');
        fprintf('✅ Summary statistics working\n');
        
    else
        fprintf('\n❌ No results returned\n');
    end
    
catch ME
    fprintf('\n❌ Core pipeline failed: %s\n', ME.message);
    fprintf('Stack trace:\n');
    for i = 1:length(ME.stack)
        fprintf('  %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
    end
end

fprintf('\n=== Next Steps ===\n');
fprintf('If this test PASSES:\n');
fprintf('  → Issue is in the plotting pipeline (baseline_plotter.m)\n');
fprintf('  → Core processing works perfectly\n');
fprintf('  → We need to debug the visualization functions\n');
fprintf('\nIf this test FAILS:\n');
fprintf('  → Issue is in core data processing\n');
fprintf('  → Need to debug the main pipeline logic\n');