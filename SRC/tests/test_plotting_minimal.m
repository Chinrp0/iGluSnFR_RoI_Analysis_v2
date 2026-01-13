%% MINIMAL PLOTTING TEST - Test baseline_plotter step by step
% This isolates plotting issues from the core processing

clear; clc; close all;

fprintf('=== Testing Baseline Plotting Components ===\n');

%% Setup
addpath(genpath(pwd));
test_file = 'D:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglusnfr4f_NGR\Spont\GPU_SNR_Processed\5_raw_mean\CP_Snfr4-NGR_Doc2b-R213W_Cs1-c1_spont-01_mean.csv';

%% Load and process data (we know this works)
fprintf('Loading data and computing baseline...\n');
loader = csv_loader_v2();
[data, metadata] = loader.loadSingleFile(test_file);
config = tracenorm_config();

[baseline, outlier_mask, baseline_stats] = baseline_detector(data, config);
[dfof_data, dfof_stats] = dfof_calculator(data, baseline, config);

fprintf('✓ Data loaded and processed successfully\n');
fprintf('  Data: %d frames × %d ROIs\n', size(data, 1), size(data, 2));

%% Test 1: Can we create a simple figure?
fprintf('\nTest 1: Basic figure creation...\n');
try
    test_fig = figure('Name', 'Test Figure');
    plot(1:10, rand(1,10));
    title('Basic Test Plot');
    fprintf('✓ Basic figure creation works\n');
    close(test_fig);
catch ME
    fprintf('✗ Basic figure creation failed: %s\n', ME.message);
    return;
end

%% Test 2: Check required data structures
fprintf('\nTest 2: Checking data structures...\n');

required_baseline_fields = {'fraction_good_rois', 'num_transport_rois', 'outlier_fraction', 'valid_fraction'};
for i = 1:length(required_baseline_fields)
    field = required_baseline_fields{i};
    if isfield(baseline_stats, field)
        fprintf('  ✓ baseline_stats.%s exists\n', field);
    else
        fprintf('  ✗ baseline_stats.%s MISSING\n', field);
    end
end

required_dfof_fields = {'mean_snr', 'fraction_good_rois', 'snr', 'events_per_roi', 'event_summary'};
for i = 1:length(required_dfof_fields)
    field = required_dfof_fields{i};
    if isfield(dfof_stats, field)
        fprintf('  ✓ dfof_stats.%s exists\n', field);
    else
        fprintf('  ✗ dfof_stats.%s MISSING\n', field);
    end
end

%% Test 3: Test individual plotting functions
fprintf('\nTest 3: Testing plotting functions individually...\n');

% Test create_trace_overview (simplified version)
fprintf('  Testing trace overview creation...\n');
try
    time_vector = (1:size(data, 1)) / config.frame_rate;
    
    % Create simplified trace overview
    fig1 = figure('Name', 'Trace Overview Test', 'Position', [100, 100, 800, 600]);
    
    % Select first 3 ROIs
    sample_rois = 1:min(3, size(data, 2));
    
    for i = 1:length(sample_rois)
        roi_idx = sample_rois(i);
        
        subplot(3, 2, 2*i-1);
        plot(time_vector, data(:, roi_idx), 'k-', 'LineWidth', 0.8);
        hold on;
        plot(time_vector, baseline(:, roi_idx), 'r-', 'LineWidth', 1.2);
        xlabel('Time (s)');
        ylabel('Fluorescence');
        title(sprintf('ROI %d - Raw + Baseline', roi_idx));
        grid on;
        
        subplot(3, 2, 2*i);
        plot(time_vector, dfof_data(:, roi_idx), 'b-', 'LineWidth', 1.0);
        xlabel('Time (s)');
        ylabel('dF/F');
        title(sprintf('ROI %d - dF/F', roi_idx));
        grid on;
    end
    
    sgtitle('Simplified Trace Overview');
    fprintf('  ✓ Trace overview creation successful\n');
    
catch ME
    fprintf('  ✗ Trace overview failed: %s\n', ME.message);
    fprintf('    Error in: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
end

%% Test 4: Test validation plots
fprintf('  Testing validation plots...\n');
try
    fig2 = figure('Name', 'Validation Test', 'Position', [200, 100, 800, 600]);
    
    % Test histogram creation with baseline stats
    subplot(2, 2, 1);
    if isfield(baseline_stats, 'baseline_cv') && ~isempty(baseline_stats.baseline_cv)
        histogram(baseline_stats.baseline_cv, 20);
        title('Baseline CV');
        xlabel('Coefficient of Variation');
        ylabel('Count');
    else
        text(0.5, 0.5, 'baseline_cv data missing', 'HorizontalAlignment', 'center');
        title('Baseline CV - NO DATA');
    end
    
    subplot(2, 2, 2);
    if isfield(dfof_stats, 'snr') && ~isempty(dfof_stats.snr)
        histogram(dfof_stats.snr, 20);
        title('SNR Distribution');
        xlabel('SNR');
        ylabel('Count');
    else
        text(0.5, 0.5, 'SNR data missing', 'HorizontalAlignment', 'center');
        title('SNR - NO DATA');
    end
    
    subplot(2, 2, 3);
    if isfield(baseline_stats, 'outlier_fraction') && ~isempty(baseline_stats.outlier_fraction)
        histogram(baseline_stats.outlier_fraction, 20);
        title('Outlier Fraction');
        xlabel('Fraction');
        ylabel('Count');
    else
        text(0.5, 0.5, 'Outlier fraction missing', 'HorizontalAlignment', 'center');
        title('Outlier Fraction - NO DATA');
    end
    
    subplot(2, 2, 4);
    if isfield(dfof_stats, 'events_per_roi') && ~isempty(dfof_stats.events_per_roi)
        histogram(dfof_stats.events_per_roi, 0:max(dfof_stats.events_per_roi));
        title('Events per ROI');
        xlabel('Event Count');
        ylabel('ROI Count');
    else
        text(0.5, 0.5, 'Events per ROI missing', 'HorizontalAlignment', 'center');
        title('Events per ROI - NO DATA');
    end
    
    sgtitle('Validation Plots Test');
    fprintf('  ✓ Validation plots successful\n');
    
catch ME
    fprintf('  ✗ Validation plots failed: %s\n', ME.message);
    fprintf('    Error in: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
end

%% Test 5: Try calling actual baseline_plotter
fprintf('  Testing actual baseline_plotter function...\n');
try
    plot_handles = baseline_plotter(data, baseline, dfof_data, outlier_mask, ...
        baseline_stats, metadata, config);
    fprintf('  ✓ baseline_plotter completed successfully!\n');
    fprintf('    Created %d plots\n', length(fieldnames(plot_handles)));
    
catch ME
    fprintf('  ✗ baseline_plotter failed: %s\n', ME.message);
    fprintf('    Error in: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
    
    % Print more detailed error info
    if length(ME.stack) > 1
        fprintf('    Call stack:\n');
        for i = 1:min(3, length(ME.stack))
            fprintf('      %d: %s (line %d)\n', i, ME.stack(i).name, ME.stack(i).line);
        end
    end
end

fprintf('\n=== Plotting Test Summary ===\n');
fprintf('If all tests pass: Plotting should work in main pipeline\n');
fprintf('If baseline_plotter fails: We found the specific issue\n');
fprintf('Check the error messages above for specific problems\n');