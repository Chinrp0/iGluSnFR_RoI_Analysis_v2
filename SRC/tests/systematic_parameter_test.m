%% SYSTEMATIC PARAMETER TESTING
% Test upper threshold variations (3.0-5.0σ) and minimum frame lengths (3-6)
% Focus on ROIs 636 and 627

clear; clc; close all;

fprintf('=== Systematic Schmitt Trigger Parameter Testing ===\n');

%% === Load Data ===
test_file = 'D:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglusnfr4f_NGR\Spont\GPU_SNR_Processed\5_raw_mean\CP_Snfr4-NGR_Doc2b-R213W_Cs1-c1_spont-01_mean.csv';

try
    loader = csv_loader_v2();
    [data, metadata] = loader.loadSingleFile(test_file);
    config = tracenorm_config();
    
    [baseline, outlier_mask, baseline_stats] = baseline_detector(data, config);
    [dfof_data, dfof_stats] = dfof_calculator(data, baseline, config);
    
    fprintf('✓ Data loaded: %d frames × %d ROIs\n', size(data, 1), size(data, 2));
catch ME
    fprintf('✗ Error loading data: %s\n', ME.message);
    return;
end

%% === Test ROIs ===
test_rois = [636, 627];
time_vector = (1:size(dfof_data, 1)) / config.frame_rate;

%% === Parameter Ranges ===
% Fixed lower threshold at 1.5σ
lower_sigma = 1.5;

% Test upper threshold variations
upper_sigmas = [3.0, 3.5, 4.0, 4.5, 5.0];

% Test minimum frame length variations  
min_frame_lengths = [3, 4, 5, 6];

fprintf('\nTesting Parameters:\n');
fprintf('  Lower threshold: %.1fσ (fixed)\n', lower_sigma);
fprintf('  Upper thresholds: '); fprintf('%.1f ', upper_sigmas); fprintf('σ\n');
fprintf('  Min frame lengths: '); fprintf('%d ', min_frame_lengths); fprintf(' frames\n');

%% === Test Each ROI ===
for roi_idx = 1:length(test_rois)
    roi = test_rois(roi_idx);
    roi_trace = dfof_data(:, roi);
    
    fprintf('\n=== Testing ROI %d ===\n', roi);
    
    % Calculate noise and thresholds
    valid_data = roi_trace(~isnan(roi_trace));
    if length(valid_data) > 20
        baseline_70th = prctile(valid_data, 70);
        baseline_data = valid_data(valid_data <= baseline_70th);
        noise_std = mad(baseline_data, 1) * 1.4826;
    else
        noise_std = mad(valid_data, 1) * 1.4826;
    end
    
    lower_threshold = lower_sigma * noise_std;
    
    fprintf('ROI %d characteristics:\n', roi);
    fprintf('  Noise std: %.6f\n', noise_std);
    fprintf('  Lower threshold (%.1fσ): %.6f\n', lower_sigma, lower_threshold);
    fprintf('  Trace range: [%.6f, %.6f]\n', min(valid_data), max(valid_data));
    
    %% === Test 1: Upper Threshold Variation (fixed 3-frame minimum) ===
    fprintf('\n--- Test 1: Upper Threshold Variation (3-frame minimum) ---\n');
    
    upper_results = struct();
    fixed_min_frames = 3;
    
    for u = 1:length(upper_sigmas)
        upper_sigma = upper_sigmas(u);
        upper_threshold = upper_sigma * noise_std;
        
        events = test_schmitt_parameters(roi_trace, upper_threshold, lower_threshold, fixed_min_frames);
        event_starts = find(diff([false; events]) == 1);
        event_ends = find(diff([events; false]) == -1);
        
        upper_results(u).upper_sigma = upper_sigma;
        upper_results(u).upper_threshold = upper_threshold;
        upper_results(u).events = events;
        upper_results(u).starts = event_starts;
        upper_results(u).ends = event_ends;
        upper_results(u).count = length(event_starts);
        
        if ~isempty(event_starts)
            durations = event_ends - event_starts + 1;
            upper_results(u).avg_duration = mean(durations);
            upper_results(u).total_frames = sum(durations);
        else
            upper_results(u).avg_duration = 0;
            upper_results(u).total_frames = 0;
        end
        
        fprintf('  %.1fσ upper: %2d events, avg %.1f frames\n', ...
            upper_sigma, upper_results(u).count, upper_results(u).avg_duration);
    end
    
    %% === Test 2: Minimum Frame Length Variation (fixed 3.0σ upper) ===
    fprintf('\n--- Test 2: Minimum Frame Length Variation (3.0σ upper) ---\n');
    
    frame_results = struct();
    fixed_upper_sigma = 3.0;
    fixed_upper_threshold = fixed_upper_sigma * noise_std;
    
    for f = 1:length(min_frame_lengths)
        min_frames = min_frame_lengths(f);
        
        events = test_schmitt_parameters(roi_trace, fixed_upper_threshold, lower_threshold, min_frames);
        event_starts = find(diff([false; events]) == 1);
        event_ends = find(diff([events; false]) == -1);
        
        frame_results(f).min_frames = min_frames;
        frame_results(f).events = events;
        frame_results(f).starts = event_starts;
        frame_results(f).ends = event_ends;
        frame_results(f).count = length(event_starts);
        
        if ~isempty(event_starts)
            durations = event_ends - event_starts + 1;
            frame_results(f).avg_duration = mean(durations);
            frame_results(f).total_frames = sum(durations);
        else
            frame_results(f).avg_duration = 0;
            frame_results(f).total_frames = 0;
        end
        
        fprintf('  %d-frame min: %2d events, avg %.1f frames\n', ...
            min_frames, frame_results(f).count, frame_results(f).avg_duration);
    end
    
    %% === Visualization for this ROI ===
    fig = figure('Name', sprintf('Parameter Testing - ROI %d', roi), ...
        'Position', [100 + (roi_idx-1)*50, 100 + (roi_idx-1)*50, 1600, 1200]);
    
    % Plot 1: Upper threshold comparison
    subplot(3, 2, 1);
    plot(time_vector, roi_trace, 'k-', 'LineWidth', 1);
    hold on;
    
    % Show baseline (3.0σ) result
    baseline_idx = find(upper_sigmas == 3.0);
    if ~isempty(baseline_idx) && any(upper_results(baseline_idx).events)
        baseline_trace = roi_trace;
        baseline_trace(~upper_results(baseline_idx).events) = NaN;
        plot(time_vector, baseline_trace, 'r-', 'LineWidth', 2);
    end
    
    % Add thresholds
    line([time_vector(1), time_vector(end)], [upper_results(baseline_idx).upper_threshold, upper_results(baseline_idx).upper_threshold], ...
        'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.5);
    line([time_vector(1), time_vector(end)], [lower_threshold, lower_threshold], ...
        'Color', 'c', 'LineStyle', '--', 'LineWidth', 1.5);
    
    title(sprintf('ROI %d: Baseline (3.0σ/1.5σ, 3-frame)', roi));
    ylabel('dF/F'); legend('Trace', 'Events', 'Upper', 'Lower', 'Location', 'best');
    grid on;
    
    % Plot 2: Upper threshold event counts
    subplot(3, 2, 2);
    event_counts_upper = [upper_results.count];
    bar(event_counts_upper, 'FaceColor', [0.7 0.3 0.3], 'EdgeColor', 'k');
    set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%.1fσ', x), upper_sigmas, 'UniformOutput', false));
    ylabel('Event Count'); title('Upper Threshold Sensitivity');
    grid on;
    % Add count labels
    for i = 1:length(event_counts_upper)
        text(i, event_counts_upper(i) + 0.3, sprintf('%d', event_counts_upper(i)), ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
    
    % Plot 3: Most sensitive upper threshold (3.0σ)
    subplot(3, 2, 3);
    plot(time_vector, roi_trace, 'k-', 'LineWidth', 1);
    hold on;
    if any(upper_results(1).events)  % 3.0σ result
        sensitive_trace = roi_trace;
        sensitive_trace(~upper_results(1).events) = NaN;
        plot(time_vector, sensitive_trace, 'b-', 'LineWidth', 2);
    end
    
    line([time_vector(1), time_vector(end)], [upper_results(1).upper_threshold, upper_results(1).upper_threshold], ...
        'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.5);
    line([time_vector(1), time_vector(end)], [lower_threshold, lower_threshold], ...
        'Color', 'c', 'LineStyle', '--', 'LineWidth', 1.5);
    
    title(sprintf('Most Sensitive (%.1fσ/1.5σ): %d events', upper_sigmas(1), upper_results(1).count));
    ylabel('dF/F'); legend('Trace', 'Events', 'Upper', 'Lower', 'Location', 'best');
    grid on;
    
    % Plot 4: Frame length event counts
    subplot(3, 2, 4);
    event_counts_frames = [frame_results.count];
    bar(event_counts_frames, 'FaceColor', [0.3 0.7 0.3], 'EdgeColor', 'k');
    set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%d-frame', x), min_frame_lengths, 'UniformOutput', false));
    ylabel('Event Count'); title('Minimum Frame Length Effect');
    grid on;
    % Add count labels
    for i = 1:length(event_counts_frames)
        text(i, event_counts_frames(i) + 0.3, sprintf('%d', event_counts_frames(i)), ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
    
    % Plot 5: Strictest frame requirement
    subplot(3, 2, 5);
    plot(time_vector, roi_trace, 'k-', 'LineWidth', 1);
    hold on;
    strictest_idx = length(frame_results);  % 6-frame requirement
    if any(frame_results(strictest_idx).events)
        strict_trace = roi_trace;
        strict_trace(~frame_results(strictest_idx).events) = NaN;
        plot(time_vector, strict_trace, 'm-', 'LineWidth', 2);
    end
    
    line([time_vector(1), time_vector(end)], [fixed_upper_threshold, fixed_upper_threshold], ...
        'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.5);
    line([time_vector(1), time_vector(end)], [lower_threshold, lower_threshold], ...
        'Color', 'c', 'LineStyle', '--', 'LineWidth', 1.5);
    
    title(sprintf('Strictest (3.0σ/1.5σ, %d-frame): %d events', ...
        min_frame_lengths(end), frame_results(strictest_idx).count));
    xlabel('Time (s)'); ylabel('dF/F');
    legend('Trace', 'Events', 'Upper', 'Lower', 'Location', 'best');
    grid on;
    
    % Plot 6: Summary comparison
    subplot(3, 2, 6);
    
    % Create comparison data
    param_names = {};
    param_counts = [];
    param_colors = [];
    
    % Add upper threshold results
    for u = 1:length(upper_sigmas)
        param_names{end+1} = sprintf('%.1fσ/1.5σ', upper_sigmas(u));
        param_counts(end+1) = upper_results(u).count;
        param_colors(end+1, :) = [0.7, 0.3, 0.3];
    end
    
    % Add frame length results (excluding 3-frame since it's same as 3.0σ baseline)
    for f = 2:length(min_frame_lengths)
        param_names{end+1} = sprintf('3.0σ/%d-fr', min_frame_lengths(f));
        param_counts(end+1) = frame_results(f).count;
        param_colors(end+1, :) = [0.3, 0.7, 0.3];
    end
    
    b = bar(param_counts);
    b.FaceColor = 'flat';
    b.CData = param_colors;
    b.EdgeColor = 'k';
    
    set(gca, 'XTickLabel', param_names, 'XTickLabelRotation', 45);
    ylabel('Event Count'); title('Parameter Comparison Summary');
    grid on;
    
    % Add count labels
    for i = 1:length(param_counts)
        text(i, param_counts(i) + 0.3, sprintf('%d', param_counts(i)), ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 8);
    end
    
    sgtitle(sprintf('ROI %d Parameter Testing: Upper Threshold & Frame Length Effects', roi), ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    %% === Save results for this ROI ===
    roi_results(roi_idx).roi = roi;
    roi_results(roi_idx).noise_std = noise_std;
    roi_results(roi_idx).upper_tests = upper_results;
    roi_results(roi_idx).frame_tests = frame_results;
end

%% === Cross-ROI Summary ===
fprintf('\n=== CROSS-ROI PARAMETER SUMMARY ===\n');

% Summary table for upper thresholds
fprintf('\nUpper Threshold Effects (1.5σ lower, 3-frame minimum):\n');
fprintf('Upper σ    ROI 636    ROI 627    Combined\n');
fprintf('-------  ---------  ---------  ----------\n');
for u = 1:length(upper_sigmas)
    roi636_count = roi_results(1).upper_tests(u).count;
    roi627_count = roi_results(2).upper_tests(u).count;
    combined = roi636_count + roi627_count;
    fprintf(' %5.1f    %8d   %8d      %6d\n', upper_sigmas(u), roi636_count, roi627_count, combined);
end

% Summary table for frame lengths
fprintf('\nFrame Length Effects (3.0σ/1.5σ thresholds):\n');
fprintf('Frames     ROI 636    ROI 627    Combined\n');
fprintf('-------  ---------  ---------  ----------\n');
for f = 1:length(min_frame_lengths)
    roi636_count = roi_results(1).frame_tests(f).count;
    roi627_count = roi_results(2).frame_tests(f).count;
    combined = roi636_count + roi627_count;
    fprintf(' %5d    %8d   %8d      %6d\n', min_frame_lengths(f), roi636_count, roi627_count, combined);
end

%% === Recommendations ===
fprintf('\n=== PARAMETER RECOMMENDATIONS ===\n');

% Find optimal parameters
target_events_per_roi = [3, 12];  % Reasonable range

fprintf('Target: %d-%d events per ROI for biological realism\n', target_events_per_roi(1), target_events_per_roi(2));

% Check upper threshold options
fprintf('\nUpper threshold recommendations:\n');
for u = 1:length(upper_sigmas)
    roi636_ok = roi_results(1).upper_tests(u).count >= target_events_per_roi(1) && ...
                roi_results(1).upper_tests(u).count <= target_events_per_roi(2);
    roi627_ok = roi_results(2).upper_tests(u).count >= target_events_per_roi(1) && ...
                roi_results(2).upper_tests(u).count <= target_events_per_roi(2);
    
    if roi636_ok && roi627_ok
        fprintf('  ✓ %.1fσ upper: Both ROIs in target range\n', upper_sigmas(u));
    elseif roi636_ok || roi627_ok
        fprintf('  ~ %.1fσ upper: One ROI in target range\n', upper_sigmas(u));
    else
        fprintf('  ✗ %.1fσ upper: Both ROIs outside target range\n', upper_sigmas(u));
    end
end

% Check frame length options
fprintf('\nFrame length recommendations:\n');
for f = 1:length(min_frame_lengths)
    roi636_ok = roi_results(1).frame_tests(f).count >= target_events_per_roi(1) && ...
                roi_results(1).frame_tests(f).count <= target_events_per_roi(2);
    roi627_ok = roi_results(2).frame_tests(f).count >= target_events_per_roi(1) && ...
                roi_results(2).frame_tests(f).count <= target_events_per_roi(2);
    
    if roi636_ok && roi627_ok
        fprintf('  ✓ %d-frame minimum: Both ROIs in target range\n', min_frame_lengths(f));
    elseif roi636_ok || roi627_ok
        fprintf('  ~ %d-frame minimum: One ROI in target range\n', min_frame_lengths(f));
    else
        fprintf('  ✗ %d-frame minimum: Both ROIs outside target range\n', min_frame_lengths(f));
    end
end

fprintf('\n✅ Systematic parameter testing complete!\n');

%% === Helper Function ===
function events = test_schmitt_parameters(trace, upper_thresh, lower_thresh, min_frames)
    % Test Schmitt trigger with specific parameters
    
    events = false(size(trace));
    state = 'baseline';
    event_start = 0;
    validated_events = [];
    
    for frame = 1:length(trace)
        signal = trace(frame);
        if isnan(signal), continue; end
        
        if strcmp(state, 'baseline')
            if signal > upper_thresh
                state = 'validating';
                event_start = frame;
            end
            
        elseif strcmp(state, 'validating')
            frames_since_start = frame - event_start + 1;
            
            if signal < lower_thresh
                % Validation failed
                state = 'baseline';
                event_start = 0;
            elseif frames_since_start >= min_frames
                % Validation passed
                state = 'validated_event';
            end
            
        elseif strcmp(state, 'validated_event')
            if signal < lower_thresh
                event_end = frame - 1;
                validated_events(end+1, :) = [event_start, event_end];
                state = 'baseline';
                event_start = 0;
            end
        end
    end
    
    % Handle event extending to end
    if strcmp(state, 'validated_event') && event_start > 0
        validated_events(end+1, :) = [event_start, length(trace)];
    end
    
    % Create event mask with merging
    for i = 1:size(validated_events, 1)
        events(validated_events(i, 1):validated_events(i, 2)) = true;
    end
end
