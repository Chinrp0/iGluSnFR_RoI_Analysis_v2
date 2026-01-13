%% FOCUSED FRAME LENGTH COMPARISON
% Test 3.5σ/1.5σ thresholds with 3-frame vs 4-frame minimum duration
% ROIs 636 and 627

clear; clc; close all;

fprintf('=== Focused Frame Length Comparison ===\n');
fprintf('Parameters: 3.5σ upper / 1.5σ lower\n');
fprintf('Comparing: 3-frame vs 4-frame minimum duration\n');
fprintf('ROIs: 636 and 627\n\n');

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

%% === Fixed Parameters ===
test_rois = [636, 627];
upper_sigma = 3.5;
lower_sigma = 1.5;
frame_lengths = [3, 4];
time_vector = (1:size(dfof_data, 1)) / config.frame_rate;

%% === Process Each ROI ===
results = struct();

for roi_idx = 1:length(test_rois)
    roi = test_rois(roi_idx);
    roi_trace = dfof_data(:, roi);
    
    fprintf('=== Processing ROI %d ===\n', roi);
    
    % Calculate thresholds
    valid_data = roi_trace(~isnan(roi_trace));
    if length(valid_data) > 20
        baseline_70th = prctile(valid_data, 70);
        baseline_data = valid_data(valid_data <= baseline_70th);
        noise_std = mad(baseline_data, 1) * 1.4826;
    else
        noise_std = mad(valid_data, 1) * 1.4826;
    end
    
    upper_threshold = upper_sigma * noise_std;
    lower_threshold = lower_sigma * noise_std;
    
    fprintf('  Noise std: %.6f\n', noise_std);
    fprintf('  Upper threshold (%.1fσ): %.6f\n', upper_sigma, upper_threshold);
    fprintf('  Lower threshold (%.1fσ): %.6f\n', lower_sigma, lower_threshold);
    fprintf('  Trace range: [%.6f, %.6f]\n', min(valid_data), max(valid_data));
    
    % Test both frame lengths
    for frame_idx = 1:length(frame_lengths)
        min_frames = frame_lengths(frame_idx);
        
        fprintf('\n--- Testing %d-frame minimum ---\n', min_frames);
        
        % Apply Schmitt trigger with validation
        events = schmitt_with_frame_validation(roi_trace, upper_threshold, lower_threshold, ...
            min_frames, roi, min_frames == 3);  % Debug for 3-frame only
        
        event_starts = find(diff([false; events]) == 1);
        event_ends = find(diff([events; false]) == -1);
        
        fprintf('  Detected %d events\n', length(event_starts));
        for i = 1:min(5, length(event_starts))
            duration = event_ends(i) - event_starts(i) + 1;
            fprintf('    Event %d: [%d-%d] = %d frames (%.0f ms)\n', ...
                i, event_starts(i), event_ends(i), duration, duration * 1000 / config.frame_rate);
        end
        
        % Store results
        results(roi_idx, frame_idx).roi = roi;
        results(roi_idx, frame_idx).min_frames = min_frames;
        results(roi_idx, frame_idx).events = events;
        results(roi_idx, frame_idx).starts = event_starts;
        results(roi_idx, frame_idx).ends = event_ends;
        results(roi_idx, frame_idx).count = length(event_starts);
        results(roi_idx, frame_idx).upper_threshold = upper_threshold;
        results(roi_idx, frame_idx).lower_threshold = lower_threshold;
        
        if ~isempty(event_starts)
            durations = event_ends - event_starts + 1;
            results(roi_idx, frame_idx).avg_duration = mean(durations);
            results(roi_idx, frame_idx).total_frames = sum(durations);
        else
            results(roi_idx, frame_idx).avg_duration = 0;
            results(roi_idx, frame_idx).total_frames = 0;
        end
    end
end

%% === Create Comprehensive Visualization ===
fig = figure('Name', 'Frame Length Comparison: 3.5σ/1.5σ Thresholds', ...
    'Position', [100, 100, 1800, 1200]);

% Plot layout: 2 ROIs × 2 frame lengths = 4 main plots + 1 comparison
subplot_positions = [
    1, 2;  % ROI 636: 3-frame, 4-frame
    3, 4;  % ROI 627: 3-frame, 4-frame
];

for roi_idx = 1:length(test_rois)
    roi = test_rois(roi_idx);
    roi_trace = dfof_data(:, roi);
    
    for frame_idx = 1:length(frame_lengths)
        subplot_idx = subplot_positions(roi_idx, frame_idx);
        subplot(3, 2, subplot_idx);
        
        % Plot trace
        plot(time_vector, roi_trace, 'k-', 'LineWidth', 1.0);
        hold on;
        
        % Highlight detected events
        events = results(roi_idx, frame_idx).events;
        if any(events)
            event_trace = roi_trace;
            event_trace(~events) = NaN;
            
            if frame_idx == 1
                color = 'b';  % Blue for 3-frame
                linewidth = 2.5;
            else
                color = 'r';  % Red for 4-frame
                linewidth = 2.0;
            end
            
            plot(time_vector, event_trace, color, 'LineWidth', linewidth);
        end
        
        % Add threshold lines
        upper_thresh = results(roi_idx, frame_idx).upper_threshold;
        lower_thresh = results(roi_idx, frame_idx).lower_threshold;
        
        line([time_vector(1), time_vector(end)], [upper_thresh, upper_thresh], ...
            'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.5);
        line([time_vector(1), time_vector(end)], [lower_thresh, lower_thresh], ...
            'Color', 'c', 'LineStyle', '--', 'LineWidth', 1.5);
        line([time_vector(1), time_vector(end)], [0, 0], ...
            'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 1);
        
        % Labels and title
        xlabel('Time (s)'); ylabel('dF/F');
        title(sprintf('ROI %d: %d-frame minimum (%d events)', ...
            roi, frame_lengths(frame_idx), results(roi_idx, frame_idx).count));
        
        legend('dF/F', 'Detected Events', 'Upper (3.5σ)', 'Lower (1.5σ)', 'Zero', ...
            'Location', 'best');
        grid on;
    end
end

% Comparison plots
% Plot 5: Side-by-side event count comparison
subplot(3, 2, 5);
roi_names = arrayfun(@(x) sprintf('ROI %d', x), test_rois, 'UniformOutput', false);
frame_3_counts = [results(1,1).count, results(2,1).count];
frame_4_counts = [results(1,2).count, results(2,2).count];

x = 1:length(test_rois);
width = 0.35;

b1 = bar(x - width/2, frame_3_counts, width, 'FaceColor', [0.3 0.6 0.9], 'EdgeColor', 'k');
hold on;
b2 = bar(x + width/2, frame_4_counts, width, 'FaceColor', [0.9 0.3 0.3], 'EdgeColor', 'k');

% Add count labels
for i = 1:length(frame_3_counts)
    text(i - width/2, frame_3_counts(i) + 0.2, sprintf('%d', frame_3_counts(i)), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    text(i + width/2, frame_4_counts(i) + 0.2, sprintf('%d', frame_4_counts(i)), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

set(gca, 'XTickLabel', roi_names);
ylabel('Event Count');
title('Event Count Comparison');
legend('3-frame minimum', '4-frame minimum', 'Location', 'best');
grid on;

% Plot 6: Duration comparison
subplot(3, 2, 6);
frame_3_durations = [results(1,1).avg_duration, results(2,1).avg_duration];
frame_4_durations = [results(1,2).avg_duration, results(2,2).avg_duration];

b1 = bar(x - width/2, frame_3_durations, width, 'FaceColor', [0.3 0.6 0.9], 'EdgeColor', 'k');
hold on;
b2 = bar(x + width/2, frame_4_durations, width, 'FaceColor', [0.9 0.3 0.3], 'EdgeColor', 'k');

% Add duration labels
for i = 1:length(frame_3_durations)
    if frame_3_durations(i) > 0
        text(i - width/2, frame_3_durations(i) + 0.1, sprintf('%.1f', frame_3_durations(i)), ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
    if frame_4_durations(i) > 0
        text(i + width/2, frame_4_durations(i) + 0.1, sprintf('%.1f', frame_4_durations(i)), ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
end

set(gca, 'XTickLabel', roi_names);
ylabel('Average Duration (frames)');
title('Average Event Duration');
legend('3-frame minimum', '4-frame minimum', 'Location', 'best');
grid on;

sgtitle('Frame Length Comparison: 3.5σ/1.5σ Thresholds', 'FontSize', 16, 'FontWeight', 'bold');

%% === Summary Table ===
fprintf('\n=== COMPARISON SUMMARY ===\n');
fprintf('Parameters: %.1fσ upper / %.1fσ lower thresholds\n', upper_sigma, lower_sigma);
fprintf('\n');
fprintf('ROI    Frame Min    Events    Avg Duration    Total Frames    Effect\n');
fprintf('-----  ----------  --------  --------------  ------------  ----------\n');

for roi_idx = 1:length(test_rois)
    roi = test_rois(roi_idx);
    
    % 3-frame results
    count_3 = results(roi_idx, 1).count;
    dur_3 = results(roi_idx, 1).avg_duration;
    total_3 = results(roi_idx, 1).total_frames;
    
    % 4-frame results  
    count_4 = results(roi_idx, 2).count;
    dur_4 = results(roi_idx, 2).avg_duration;
    total_4 = results(roi_idx, 2).total_frames;
    
    fprintf(' %3d      3-frame    %8d      %8.1f        %8d   Baseline\n', ...
        roi, count_3, dur_3, total_3);
    
    % Calculate effect
    if count_3 > 0
        count_change = count_4 - count_3;
        if count_change == 0
            effect_str = 'No change';
        elseif count_change > 0
            effect_str = sprintf('+%d events', count_change);
        else
            effect_str = sprintf('%d events', count_change);
        end
    else
        effect_str = 'N/A';
    end
    
    fprintf(' %3d      4-frame    %8d      %8.1f        %8d   %s\n', ...
        roi, count_4, dur_4, total_4, effect_str);
    fprintf('\n');
end

%% === Combined Analysis ===
total_3_frame = sum([results(:,1).count]);
total_4_frame = sum([results(:,2).count]);

fprintf('=== COMBINED ANALYSIS ===\n');
fprintf('Total events (both ROIs):\n');
fprintf('  3-frame minimum: %d events\n', total_3_frame);
fprintf('  4-frame minimum: %d events\n', total_4_frame);
fprintf('  Difference: %+d events\n', total_4_frame - total_3_frame);

if total_3_frame > 0
    percent_change = 100 * (total_4_frame - total_3_frame) / total_3_frame;
    fprintf('  Percent change: %+.1f%%\n', percent_change);
end

%% === Recommendations ===
fprintf('\n=== RECOMMENDATIONS ===\n');

if total_4_frame < total_3_frame
    fprintf('✓ 4-frame minimum reduces events by %d (%.1f%% reduction)\n', ...
        total_3_frame - total_4_frame, 100*(total_3_frame - total_4_frame)/total_3_frame);
    fprintf('  → Better noise rejection, fewer false positives\n');
elseif total_4_frame == total_3_frame
    fprintf('= 4-frame minimum shows no change in event count\n');
    fprintf('  → Both parameters detect same sustained events\n');
else
    fprintf('? 4-frame minimum increases events (unexpected)\n');
    fprintf('  → Check for detection artifacts\n');
end

% Check if both ROIs show consistent pattern
roi636_change = results(2,2).count - results(2,1).count;
roi627_change = results(1,2).count - results(1,1).count;

if sign(roi636_change) == sign(roi627_change) || (roi636_change == 0 && roi627_change == 0)
    fprintf('✓ Both ROIs show consistent response to 4-frame minimum\n');
else
    fprintf('⚠ ROIs show different responses - may need ROI-specific parameters\n');
end

fprintf('\n✅ Focused frame length comparison complete!\n');

%% === Helper Function ===
function events = schmitt_with_frame_validation(trace, upper_thresh, lower_thresh, min_frames, roi, show_debug)
    % Schmitt trigger with specific frame validation
    
    events = false(size(trace));
    state = 'baseline';
    event_start = 0;
    validated_events = [];
    
    if show_debug
        fprintf('    --- Debug for ROI %d (%d-frame minimum) ---\n', roi, min_frames);
    end
    
    for frame = 1:length(trace)
        signal = trace(frame);
        if isnan(signal), continue; end
        
        if strcmp(state, 'baseline')
            if signal > upper_thresh
                state = 'validating';
                event_start = frame;
                if show_debug
                    fprintf('    Potential START at frame %d: signal=%.6f > %.6f\n', ...
                        frame, signal, upper_thresh);
                end
            end
            
        elseif strcmp(state, 'validating')
            frames_since_start = frame - event_start + 1;
            
            if signal < lower_thresh
                if show_debug
                    fprintf('      VALIDATION FAILED at frame %d: signal=%.6f < %.6f after %d frames\n', ...
                        frame, signal, lower_thresh, frames_since_start);
                end
                state = 'baseline';
                event_start = 0;
            elseif frames_since_start >= min_frames
                state = 'validated_event';
                if show_debug
                    fprintf('      VALIDATION PASSED at frame %d: continuing event...\n', frame);
                end
            end
            
        elseif strcmp(state, 'validated_event')
            if signal < lower_thresh
                event_end = frame - 1;
                validated_events(end+1, :) = [event_start, event_end];
                duration = event_end - event_start + 1;
                if show_debug
                    fprintf('      EVENT END at frame %d: duration=%d frames\n', event_end, duration);
                end
                state = 'baseline';
                event_start = 0;
            end
        end
    end
    
    % Handle event extending to end
    if strcmp(state, 'validated_event') && event_start > 0
        validated_events(end+1, :) = [event_start, length(trace)];
        if show_debug
            fprintf('      EVENT extends to end: duration=%d frames\n', length(trace) - event_start + 1);
        end
    end
    
    % Create event mask
    for i = 1:size(validated_events, 1)
        events(validated_events(i, 1):validated_events(i, 2)) = true;
    end
    
    if show_debug
        fprintf('    Total validated events: %d\n', size(validated_events, 1));
    end
end
