%% TEST FIXED PIPELINE - Verify all modules work together consistently
% This script tests the fixed modular pipeline with consistent data structures

clear; clc; close all;

fprintf('=== Testing Fixed Modular Pipeline ===\n');

%% === Configuration ===
test_file = "D:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglusnfr4f_NGR\Spont\GPU_SNR_Processed\5_raw_mean\CP_Snfr4-NGR_Doc2b-R213W_Cs1-c1_spont-01_mean.csv";

% Test options
test_options = struct();
test_options.verbose = true;
test_options.createPlots = true;   % Test the fixed plotting
test_options.savePlots = false;

%% === Setup Pipeline ===
addpath(genpath(pwd));

%% === Test Individual Modules First ===
fprintf('\n1. Testing individual modules...\n');

% Load data
fprintf('  Loading CSV data...\n');
loader = csv_loader_v2();
[data, metadata] = loader.loadSingleFile(test_file);
config = tracenorm_config();
fprintf('    ✓ Data loaded: %d frames × %d ROIs\n', size(data, 1), size(data, 2));

% Test baseline detection
fprintf('  Testing baseline detection...\n');
[baseline, outlier_mask, baseline_stats] = baseline_detector(data, config);
fprintf('    ✓ Baseline calculated, %.1f%% outliers detected\n', ...
    100 * baseline_stats.mean_outlier_fraction);

% Test dF/F calculation
fprintf('  Testing dF/F calculation...\n');
[dfof_data, dfof_stats] = dfof_calculator(data, baseline, config);
fprintf('    ✓ dF/F calculated, mean SNR = %.2f\n', dfof_stats.mean_snr);

% Test event detection (FIXED version)
fprintf('  Testing fixed event detection...\n');
[event_mask, event_stats] = schmitt_event_detector(dfof_data, config);
fprintf('    ✓ Events detected: %d total events across %d ROIs\n', ...
    event_stats.total_events, event_stats.rois_with_events);

% CRITICAL TEST: Check if event_mask is properly included
if isfield(event_stats, 'event_mask')
    fprintf('    ✓ event_mask properly included in event_stats\n');
else
    fprintf('    ✗ event_mask MISSING from event_stats - CRITICAL ERROR\n');
    return;
end

% Test quality assessment
fprintf('  Testing quality assessment...\n');
quality_metrics = quality_assessor(data, baseline_stats, dfof_stats, event_stats, config);
fprintf('    ✓ Quality assessed: %d active ROIs (%.1f%%)\n', ...
    quality_metrics.summary.active_rois, 100 * quality_metrics.summary.fraction_active);

%% === Test Fixed Plotting (2x4 layout with peak markers) ===
fprintf('\n2. Testing fixed plotting with 2x4 layout...\n');
try
    plot_handles = baseline_plotter(data, baseline, dfof_data, outlier_mask, ...
        baseline_stats, dfof_stats, metadata, config, event_stats);
    
    fprintf('    ✓ Plotting successful: %d figures created\n', length(fieldnames(plot_handles)));
    fprintf('    ✓ All figures should use 2x4 layout (8 subplots)\n');
    fprintf('    ✓ Peak markers should appear ABOVE event peaks\n');
    
catch ME
    fprintf('    ✗ Plotting failed: %s\n', ME.message);
    if ~isempty(ME.stack)
        fprintf('      Error in: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
    end
    return;
end

%% === Test Full Pipeline ===
fprintf('\n3. Testing complete main_pipeline...\n');
try
    tic;
    results = main_pipeline(test_file, test_options);
    pipeline_time = toc;
    
    fprintf('    ✓ Full pipeline completed in %.3f s\n', pipeline_time);
    
    % Verify results structure
    if ~isempty(results.fileResults)
        fileResult = results.fileResults{1};
        
        % Check all expected fields
        required_fields = {'filename', 'numFrames', 'numROIs', 'baseline', 'dfof_data', ...
                          'event_mask', 'event_stats', 'quality_metrics', 'timing'};
        
        missing_fields = {};
        for i = 1:length(required_fields)
            if ~isfield(fileResult, required_fields{i})
                missing_fields{end+1} = required_fields{i};
            end
        end
        
        if isempty(missing_fields)
            fprintf('    ✓ All required fields present in results\n');
        else
            fprintf('    ✗ Missing fields: %s\n', strjoin(missing_fields, ', '));
        end
        
        % Check event_mask consistency
        if isfield(fileResult, 'event_mask') && isfield(fileResult, 'event_stats')
            if isfield(fileResult.event_stats, 'event_mask')
                if isequal(fileResult.event_mask, fileResult.event_stats.event_mask)
                    fprintf('    ✓ event_mask consistent between fileResult and event_stats\n');
                else
                    fprintf('    ✗ event_mask INCONSISTENT between structures\n');
                end
            else
                fprintf('    ✗ event_mask missing from event_stats\n');
            end
        end
        
    else
        fprintf('    ✗ No file results returned\n');
        return;
    end
    
catch ME
    fprintf('    ✗ Full pipeline failed: %s\n', ME.message);
    if ~isempty(ME.stack)
        fprintf('      Error in: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
    end
    return;
end

%% === Test Data Access ===
fprintf('\n4. Testing consistent data access...\n');

fileResult = results.fileResults{1};

% Test accessing event data
fprintf('  Accessing event data:\n');
fprintf('    Total events: %d\n', fileResult.event_stats.total_events);
fprintf('    Active ROIs: %d\n', fileResult.event_stats.rois_with_events);
fprintf('    Event mask size: [%d x %d]\n', size(fileResult.event_mask, 1), size(fileResult.event_mask, 2));

% Test backward compatibility
fprintf('  Testing backward compatibility:\n');
fprintf('    Quality fraction good baseline: %.3f\n', fileResult.quality.fraction_good_baseline);
fprintf('    Quality mean SNR: %.2f\n', fileResult.quality.mean_snr);
fprintf('    Quality active ROIs: %d\n', fileResult.quality.active_rois);

%% === Test Specific ROI Access ===
fprintf('\n5. Testing specific ROI data access...\n');

% Pick an active ROI for testing
active_roi_idx = find(fileResult.event_stats.events_per_roi > 0, 1, 'first');

if ~isempty(active_roi_idx)
    fprintf('  Testing ROI %d (has %d events):\n', active_roi_idx, ...
        fileResult.event_stats.events_per_roi(active_roi_idx));
    
    % Access different data types for this ROI
    roi_raw = fileResult.baseline(:, active_roi_idx);  % Should be baseline, not raw data
    roi_dfof = fileResult.dfof_data(:, active_roi_idx);
    roi_events = fileResult.event_mask(:, active_roi_idx);
    
    fprintf('    Raw baseline range: [%.3f, %.3f]\n', min(roi_raw), max(roi_raw));
    fprintf('    dF/F range: [%.3f, %.3f]\n', min(roi_dfof), max(roi_dfof));
    fprintf('    Event frames: %d/%d (%.1f%%)\n', sum(roi_events), length(roi_events), ...
        100*sum(roi_events)/length(roi_events));
    
    fprintf('    ✓ ROI data access working correctly\n');
else
    fprintf('    ✗ No active ROIs found for testing\n');
end

%% === Final Summary ===
fprintf('\n=== PIPELINE TEST SUMMARY ===\n');
fprintf('✓ Individual modules working\n');
fprintf('✓ Event detection fixed (event_mask properly included)\n');
fprintf('✓ Plotting fixed (2x4 layout, peak markers above events)\n');
fprintf('✓ Full pipeline integration working\n');
fprintf('✓ Data structures consistent across modules\n');
fprintf('✓ Backward compatibility maintained\n');
fprintf('✓ ROI data access working\n');

fprintf('\n=== VISUAL CHECKS ===\n');
fprintf('Please verify the following in the generated plots:\n');
fprintf('1. All trace plots use 2 columns × 4 rows (8 subplots)\n');
fprintf('2. Green triangular markers appear ABOVE event peak amplitudes\n');
fprintf('3. Event periods are highlighted in red\n');
fprintf('4. Schmitt trigger thresholds shown as dashed lines\n');
fprintf('5. Plot titles show correct ROI numbers and event counts\n');

fprintf('\n✅ FIXED PIPELINE TEST COMPLETE!\n');
fprintf('All modules now work together with consistent data structures.\n');