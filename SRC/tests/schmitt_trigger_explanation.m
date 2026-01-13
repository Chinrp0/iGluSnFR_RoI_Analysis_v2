%% SCHMITT TRIGGER MECHANISM EXPLAINED
% Visual demonstration of how Schmitt trigger detection works
% Shows the difference between simple threshold and Schmitt trigger

clear; clc; close all;

fprintf('=== Schmitt Trigger Detection Mechanism ===\n');

%% === Load Real Data for Demonstration ===
addpath(genpath(pwd));
test_file = 'D:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglusnfr4f_NGR\Spont\GPU_SNR_Processed\5_raw_mean\CP_Snfr4-NGR_Doc2b-R213W_Cs1-c1_spont-01_mean.csv';

loader = csv_loader_v2();
[data, metadata] = loader.loadSingleFile(test_file);
config = tracenorm_config();

[baseline, outlier_mask, baseline_stats] = baseline_detector(data, config);
[dfof_data, dfof_stats] = dfof_calculator(data, baseline, config);

% Select an active ROI for demonstration
events_per_roi = sum(dfof_data > 0.02, 1);  % Quick event estimate
[~, active_idx] = sort(events_per_roi, 'descend');
demo_roi = active_idx(1);  % Most active ROI

roi_trace = dfof_data(:, demo_roi);
time_vector = (1:length(roi_trace)) / config.frame_rate;

fprintf('Using ROI %d for demonstration (has ~%d potential events)\n', demo_roi, events_per_roi(demo_roi));

%% === Calculate Thresholds ===
% Use same method as optimal detector
noise_std = mad(roi_trace, 1) * 1.4826;

% For demonstration, let's use slightly different parameters to show hysteresis clearly
upper_threshold = 3.0 * noise_std;  % Event START threshold
lower_threshold = 1.5 * noise_std;  % Event END threshold

fprintf('\nSchmitt trigger parameters for ROI %d:\n', demo_roi);
fprintf('  Noise level (σ): %.5f dF/F\n', noise_std);
fprintf('  Upper threshold (3.0σ): %.5f dF/F\n', upper_threshold);
fprintf('  Lower threshold (1.5σ): %.5f dF/F\n', lower_threshold);
fprintf('  Hysteresis gap: %.5f dF/F\n', upper_threshold - lower_threshold);

%% === Apply Different Detection Methods ===

% Method 1: Simple threshold (for comparison)
simple_events = roi_trace > upper_threshold;

% Method 2: Schmitt trigger (hysteresis)
schmitt_events = apply_schmitt_trigger_demo(roi_trace, upper_threshold, lower_threshold);

% Method 3: Extended Schmitt trigger (with kinetics extension)
extended_events = extend_events_demo(schmitt_events, 1, 3);  % 1 rise, 3 decay frames

fprintf('\nDetection comparison:\n');
fprintf('  Simple threshold: %d event frames\n', sum(simple_events));
fprintf('  Schmitt trigger: %d event frames\n', sum(schmitt_events));
fprintf('  Extended Schmitt: %d event frames\n', sum(extended_events));

%% === Create Visualization ===
create_schmitt_explanation_figure(roi_trace, time_vector, upper_threshold, lower_threshold, ...
    simple_events, schmitt_events, extended_events, demo_roi, metadata);

%% === Detailed Step-by-Step Analysis ===
fprintf('\n=== STEP-BY-STEP SCHMITT TRIGGER ANALYSIS ===\n');

% Find some example events for detailed analysis
schmitt_starts = find(diff([false; schmitt_events]) == 1);
schmitt_ends = find(diff([schmitt_events; false]) == -1);

fprintf('Found %d discrete Schmitt trigger events:\n', length(schmitt_starts));

% Analyze first few events in detail
num_to_analyze = min(3, length(schmitt_starts));

for e = 1:num_to_analyze
    start_frame = schmitt_starts(e);
    end_frame = schmitt_ends(e);
    start_time = time_vector(start_frame);
    end_time = time_vector(end_frame);
    duration = end_frame - start_frame + 1;
    
    % Find crossing points
    peak_frame = start_frame + find(roi_trace(start_frame:end_frame) == max(roi_trace(start_frame:end_frame)), 1) - 1;
    peak_amplitude = roi_trace(peak_frame);
    peak_time = time_vector(peak_frame);
    
    fprintf('\nEvent %d:\n', e);
    fprintf('  Start: Frame %d (%.2f s) - signal crosses above %.5f\n', start_frame, start_time, upper_threshold);
    fprintf('  Peak:  Frame %d (%.2f s) - amplitude %.5f dF/F\n', peak_frame, peak_time, peak_amplitude);
    fprintf('  End:   Frame %d (%.2f s) - signal drops below %.5f\n', end_frame, end_time, lower_threshold);
    fprintf('  Duration: %d frames (%.0f ms)\n', duration, duration * 1000 / config.frame_rate);
    fprintf('  Peak/Noise ratio: %.1fx above noise level\n', peak_amplitude / noise_std);
end

%% === What You See in Your Plots ===
fprintf('\n=== WHAT YOU SEE IN YOUR ACTUAL PLOTS ===\n');
fprintf('In your 2x4 subplot figures:\n\n');

fprintf('1. BLACK LINE: Original dF/F trace\n');
fprintf('   - Shows the actual fluorescence changes over time\n');
fprintf('   - Baseline should be around zero after normalization\n\n');

fprintf('2. RED HIGHLIGHTED REGIONS: Detected events\n');
fprintf('   - These are the periods where Schmitt trigger detected events\n');
fprintf('   - Starts when signal crosses ABOVE green dashed line (%.1fσ)\n', config.event_detection.upper_threshold_sigma);
fprintf('   - Ends when signal drops BELOW cyan dashed line (%.1fσ)\n', config.event_detection.lower_threshold_sigma);
fprintf('   - Includes kinetic extensions (%d frames before, %d frames after)\n', ...
        config.event_detection.rise_extension_frames, config.event_detection.decay_extension_frames);

fprintf('\n3. GREEN DASHED LINE: Upper threshold (%.1fσ)\n', config.event_detection.upper_threshold_sigma);
fprintf('   - Signal must cross ABOVE this line to START an event\n');
fprintf('   - Calculated as %.1f × noise_level for each ROI\n', config.event_detection.upper_threshold_sigma);
fprintf('   - Prevents false positives from small noise fluctuations\n\n');

fprintf('4. CYAN DASHED LINE: Lower threshold (%.1fσ)\n', config.event_detection.lower_threshold_sigma);
fprintf('   - Signal must drop BELOW this line to END an event\n');
fprintf('   - Calculated as %.1f × noise_level for each ROI\n', config.event_detection.lower_threshold_sigma);
fprintf('   - Ensures complete event capture including decay\n\n');

fprintf('5. ROI TITLE INFO:\n');
fprintf('   - "ROI X: Y events" = This ROI had Y discrete events detected\n');
fprintf('   - "σ=Z.ZZZZZ" = Estimated noise level (standard deviation)\n');
fprintf('   - "SNR=A.B" = Signal-to-noise ratio\n');
fprintf('   - "(iterative)" = Noise estimation method used\n\n');

%% === Why Schmitt Trigger vs Simple Threshold ===
fprintf('=== WHY SCHMITT TRIGGER vs SIMPLE THRESHOLD ===\n');
fprintf('Simple threshold problems:\n');
fprintf('  - Signal hovering around threshold creates multiple false events\n');
fprintf('  - Noise can cause brief threshold crossings\n');
fprintf('  - Doesn''t account for event kinetics (rise/decay)\n\n');

fprintf('Schmitt trigger advantages:\n');
fprintf('  - HYSTERESIS: Different start/stop thresholds prevent false events\n');
fprintf('  - KINETICS AWARE: Extends events for realistic sensor dynamics\n');
fprintf('  - NOISE ROBUST: Higher start threshold reduces false positives\n');
fprintf('  - BIOLOGY APPROPRIATE: Captures complete glutamate transients\n\n');

fprintf('In your high-noise ROIs:\n');
fprintf('  - Higher noise = higher thresholds\n');
fprintf('  - Only clear events above noise are detected\n');
fprintf('  - Red regions represent confident event detection\n\n');

fprintf('In your low-noise ROIs:\n');
fprintf('  - Lower noise = lower thresholds\n');
fprintf('  - More sensitive to smaller events\n');
fprintf('  - Can detect subtle glutamate release events\n');

fprintf('\n✅ This is why your optimal detector finds more events:\n');
fprintf('   Better noise estimation → more appropriate thresholds → better sensitivity!\n');

%% === Helper Functions ===

function schmitt_mask = apply_schmitt_trigger_demo(trace, upper_thresh, lower_thresh)
    % Demonstrate Schmitt trigger state machine
    
    numFrames = length(trace);
    schmitt_mask = false(numFrames, 1);
    
    % State machine: 'baseline' or 'in_event'
    state = 'baseline';
    event_start = 0;
    
    for frame = 1:numFrames
        signal = trace(frame);
        
        if isnan(signal)
            continue;
        end
        
        switch state
            case 'baseline'
                % Look for signal crossing ABOVE upper threshold
                if signal > upper_thresh
                    state = 'in_event';
                    event_start = frame;
                    fprintf('    Event starts at frame %d (signal=%.5f > %.5f)\n', ...
                        frame, signal, upper_thresh);
                end
                
            case 'in_event'
                % Look for signal falling BELOW lower threshold
                if signal < lower_thresh
                    % Event ends - mark all frames from start to here
                    schmitt_mask(event_start:frame) = true;
                    state = 'baseline';
                    fprintf('    Event ends at frame %d (signal=%.5f < %.5f)\n', ...
                        frame, signal, lower_thresh);
                    event_start = 0;
                end
        end
    end
    
    % Handle event extending to end of trace
    if strcmp(state, 'in_event') && event_start > 0
        schmitt_mask(event_start:end) = true;
        fprintf('    Event extends to end of trace\n');
    end
end

function extended_mask = extend_events_demo(event_mask, rise_frames, decay_frames)
    % Demonstrate event extension for sensor kinetics
    
    numFrames = length(event_mask);
    extended_mask = event_mask;
    
    if ~any(event_mask)
        return;
    end
    
    % Find event boundaries
    event_starts = find(diff([false; event_mask]) == 1);
    event_ends = find(diff([event_mask; false]) == -1);
    
    fprintf('\nExtending %d events for sensor kinetics:\n', length(event_starts));
    
    % Extend each event
    for i = 1:length(event_starts)
        start_frame = event_starts(i);
        end_frame = event_ends(i);
        
        % Extend backward for rise time
        extended_start = max(1, start_frame - rise_frames);
        
        % Extend forward for decay time
        extended_end = min(numFrames, end_frame + decay_frames);
        
        fprintf('  Event %d: [%d-%d] → [%d-%d] (+%d rise, +%d decay frames)\n', ...
            i, start_frame, end_frame, extended_start, extended_end, ...
            start_frame - extended_start, extended_end - end_frame);
        
        % Mark extended region
        extended_mask(extended_start:extended_end) = true;
    end
end

function create_schmitt_explanation_figure(trace, time_vector, upper_thresh, lower_thresh, ...
    simple_events, schmitt_events, extended_events, roi_idx, metadata)
    % Create comprehensive visualization of Schmitt trigger operation
    
    fig = figure('Name', sprintf('Schmitt Trigger Explanation - ROI %d', roi_idx), ...
        'Position', [100, 100, 1400, 1000]);
    
    % Main trace with all detection methods
    subplot(3, 1, 1);
    plot(time_vector, trace, 'k-', 'LineWidth', 1.2);
    hold on;
    
    % Add thresholds
    yline(upper_thresh, 'g--', sprintf('Upper (%.5f)', upper_thresh), 'LineWidth', 2);
    yline(lower_thresh, 'c--', sprintf('Lower (%.5f)', lower_thresh), 'LineWidth', 2);
    yline(0, 'k:', 'Alpha', 0.5);
    
    % Highlight different detection methods
    simple_trace = trace; simple_trace(~simple_events) = NaN;
    plot(time_vector, simple_trace, 'Color', [1 0.5 0], 'LineWidth', 3, 'DisplayName', 'Simple Threshold');
    
    schmitt_trace = trace; schmitt_trace(~schmitt_events) = NaN;
    plot(time_vector, schmitt_trace, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Schmitt Trigger');
    
    xlabel('Time (s)'); ylabel('dF/F');
    title(sprintf('ROI %d: Schmitt Trigger vs Simple Threshold - %s', roi_idx, metadata.filename));
    legend('dF/F', 'Upper Threshold', 'Lower Threshold', 'Zero', 'Simple', 'Schmitt', 'Location', 'best');
    grid on;
    
    % Zoomed view of first event
    subplot(3, 1, 2);
    if any(schmitt_events)
        % Find first event for zoom
        first_event_start = find(schmitt_events, 1, 'first');
        zoom_start = max(1, first_event_start - 50);
        zoom_end = min(length(trace), first_event_start + 150);
        zoom_indices = zoom_start:zoom_end;
        
        plot(time_vector(zoom_indices), trace(zoom_indices), 'k-', 'LineWidth', 1.5);
        hold on;
        
        yline(upper_thresh, 'g--', 'Upper', 'LineWidth', 2);
        yline(lower_thresh, 'c--', 'Lower', 'LineWidth', 2);
        yline(0, 'k:', 'Alpha', 0.5);
        
        % Show event progression
        schmitt_zoom = trace(zoom_indices);
        schmitt_zoom(~schmitt_events(zoom_indices)) = NaN;
        plot(time_vector(zoom_indices), schmitt_zoom, 'r-', 'LineWidth', 3);
        
        % Mark crossing points
        for i = zoom_start:zoom_end-1
            if ~schmitt_events(i) && schmitt_events(i+1)
                plot(time_vector(i+1), trace(i+1), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Event Start');
            end
            if schmitt_events(i) && ~schmitt_events(i+1)
                plot(time_vector(i), trace(i), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Event End');
            end
        end
        
        xlabel('Time (s)'); ylabel('dF/F');
        title('Zoomed View: Event Start/Stop Detection');
        legend('dF/F', 'Upper', 'Lower', 'Zero', 'Event', 'Start', 'End', 'Location', 'best');
        grid on;
    end
    
    % Comparison statistics
    subplot(3, 1, 3);
    axis off;
    
    comparison_text = {
        sprintf('DETECTION METHOD COMPARISON (ROI %d)', roi_idx);
        '';
        sprintf('Simple Threshold (signal > %.5f):', upper_thresh);
        sprintf('  Event frames: %d', sum(simple_events));
        sprintf('  Continuous regions: %d', length(find(diff([false; simple_events]) == 1)));
        '';
        sprintf('Schmitt Trigger (%.5f ↑ / %.5f ↓):', upper_thresh, lower_thresh);
        sprintf('  Event frames: %d', sum(schmitt_events));
        sprintf('  Discrete events: %d', length(find(diff([false; schmitt_events]) == 1)));
        '';
        sprintf('Extended Schmitt (+1 rise, +3 decay frames):');
        sprintf('  Event frames: %d', sum(extended_events));
        '';
        'KEY DIFFERENCES:';
        '• Simple: Any crossing creates event (noise sensitive)';
        '• Schmitt: Hysteresis prevents false events (noise robust)';
        '• Extended: Includes sensor kinetics (biologically accurate)';
        '';
        'RED REGIONS in your plots = Extended Schmitt trigger results';
    };
    
    text(0.05, 0.95, comparison_text, 'FontSize', 11, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', 'Units', 'normalized');
    
    sgtitle('How Schmitt Trigger Detection Works', 'FontSize', 16, 'FontWeight', 'bold');
end