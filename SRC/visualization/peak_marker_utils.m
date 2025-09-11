function varargout = peak_marker_utils(action, varargin)
    % PEAK_MARKER_UTILS - Utility functions for adding peak markers to plots
    % 
    % Usage:
    %   peaks = peak_marker_utils('detect', dfof_data, config)
    %   peak_marker_utils('add_to_plot', time_vector, data, peaks, plot_type)
    %   peak_marker_utils('plot_roi_with_peaks', roi_idx, data, baseline, dfof_data, peaks, time_vector)
    
    switch lower(action)
        case 'detect'
            varargout{1} = detect_peaks(varargin{:});
        case 'add_to_plot'
            add_peaks_to_plot(varargin{:});
        case 'plot_roi_with_peaks'
            plot_roi_with_peaks(varargin{:});
        case 'example'
            run_example();
        otherwise
            error('Unknown action: %s', action);
    end
end

function peak_locations = detect_peaks(dfof_data, config)
    % Detect peaks in dF/F data
    % Returns: cell array where each cell contains peak indices for that ROI
    
    if nargin < 2
        config = struct();
        config.frame_rate = 40;
    end
    
    [numFrames, numROIs] = size(dfof_data);
    peak_locations = cell(1, numROIs);
    
    % Detection parameters
    min_peak_height = 0.02;  % Minimum 2% dF/F
    min_peak_distance = round(0.5 * config.frame_rate);  % 0.5 second minimum distance
    noise_multiplier = 3;    % Peaks must be 3x above local noise
    
    fprintf('Detecting peaks in %d ROIs...\n', numROIs);
    
    for roi = 1:numROIs
        trace = dfof_data(:, roi);
        
        % Skip problematic traces
        if any(isnan(trace)) || std(trace, 'omitnan') > 0.2
            peak_locations{roi} = [];
            continue;
        end
        
        try
            % Find peaks using MATLAB's findpeaks
            [peak_vals, peak_locs] = findpeaks(trace, ...
                'MinPeakHeight', min_peak_height, ...
                'MinPeakDistance', min_peak_distance);
            
            % Additional filtering based on local noise
            if ~isempty(peak_locs)
                % Calculate local noise around each peak
                window_size = round(2 * config.frame_rate);  % 2-second window
                local_noise = movstd(trace, window_size, 'omitnan');
                
                % Keep peaks that are significantly above local noise
                valid_peaks = peak_vals > (noise_multiplier * local_noise(peak_locs));
                peak_locations{roi} = peak_locs(valid_peaks);
            else
                peak_locations{roi} = [];
            end
            
        catch ME
            % Fallback: simple threshold detection
            warning('Peak detection failed for ROI %d: %s. Using fallback.', roi, ME.message);
            threshold = max(min_peak_height, 3 * std(trace, 'omitnan'));
            above_threshold = trace > threshold;
            
            % Find start of threshold crossings
            threshold_starts = find(diff([0; above_threshold; 0]) == 1);
            peak_locations{roi} = threshold_starts;
        end
    end
    
    total_peaks = sum(cellfun(@length, peak_locations));
    fprintf('  Detected %d total peaks across all ROIs\n', total_peaks);
end

function add_peaks_to_plot(time_vector, data, peaks, plot_type)
    % Add peak markers to current plot
    % 
    % Inputs:
    %   time_vector - time points
    %   data - y-axis data (raw fluorescence or dF/F)
    %   peaks - peak indices for this ROI
    %   plot_type - 'raw' or 'dfof' (determines marker color)
    
    if isempty(peaks)
        return;
    end
    
    % Ensure we're adding to existing plot
    hold on;
    
    % Choose marker style based on plot type
    switch lower(plot_type)
        case 'raw'
            marker_color = [0 0.8 0];      % Green
            edge_color = [0 0.6 0];        % Darker green
            marker_symbol = '^';           % Upward triangle
            marker_size = 60;
            
        case 'dfof'
            marker_color = [1 0.5 0];      % Orange
            edge_color = [1 0.3 0];        % Darker orange  
            marker_symbol = '^';           % Upward triangle
            marker_size = 60;
            
        otherwise
            marker_color = [1 0 0];        % Red (default)
            edge_color = [0.8 0 0];        % Darker red
            marker_symbol = 'o';           % Circle
            marker_size = 50;
    end
    
    % Add markers
    scatter(time_vector(peaks), data(peaks), marker_size, marker_symbol, ...
        'MarkerFaceColor', marker_color, 'MarkerEdgeColor', edge_color, ...
        'LineWidth', 1.5);
end

function plot_roi_with_peaks(roi_idx, raw_data, baseline, dfof_data, peak_locations, time_vector)
    % Plot a single ROI with all traces and peak markers
    
    peaks = peak_locations{roi_idx};
    
    figure('Name', sprintf('ROI %d with Peak Detection', roi_idx), 'Position', [100, 100, 1200, 600]);
    
    % Raw trace with baseline and peaks
    subplot(2, 1, 1);
    plot(time_vector, raw_data(:, roi_idx), 'k-', 'LineWidth', 0.8, 'DisplayName', 'Raw');
    hold on;
    plot(time_vector, baseline(:, roi_idx), 'r-', 'LineWidth', 1.2, 'DisplayName', 'Baseline');
    
    % Add peak markers to raw trace
    add_peaks_to_plot(time_vector, raw_data(:, roi_idx), peaks, 'raw');
    
    xlabel('Time (s)');
    ylabel('Fluorescence');
    title(sprintf('ROI %d - Raw Trace with %d Detected Peaks', roi_idx, length(peaks)));
    legend('Raw', 'Baseline', sprintf('%d Peaks', length(peaks)), 'Location', 'best');
    grid on;
    
    % dF/F trace with peaks
    subplot(2, 1, 2);
    plot(time_vector, dfof_data(:, roi_idx), 'b-', 'LineWidth', 1.0, 'DisplayName', 'dF/F');
    hold on;
    
    % Reference lines
    yline(0, 'k--', 'Alpha', 0.5, 'DisplayName', 'Baseline');
    yline(0.02, 'g:', 'Alpha', 0.7, 'DisplayName', '2% Threshold');
    
    % Add peak markers to dF/F trace
    add_peaks_to_plot(time_vector, dfof_data(:, roi_idx), peaks, 'dfof');
    
    xlabel('Time (s)');
    ylabel('dF/F');
    title(sprintf('ROI %d - dF/F with Peak Detection (max: %.3f)', roi_idx, max(dfof_data(:, roi_idx))));
    legend('dF/F', 'Zero Line', '2% Threshold', 'Peaks', 'Location', 'best');
    grid on;
    
    % Print peak information
    if ~isempty(peaks)
        peak_times = time_vector(peaks);
        peak_amplitudes = dfof_data(peaks, roi_idx);
        
        fprintf('\nROI %d Peak Summary:\n', roi_idx);
        fprintf('  Total peaks: %d\n', length(peaks));
        fprintf('  Peak times (s): [%s]\n', sprintf('%.1f ', peak_times));
        fprintf('  Peak amplitudes (dF/F): [%s]\n', sprintf('%.3f ', peak_amplitudes));
        fprintf('  Mean amplitude: %.3f\n', mean(peak_amplitudes));
        fprintf('  Max amplitude: %.3f\n', max(peak_amplitudes));
    else
        fprintf('\nROI %d: No peaks detected\n', roi_idx);
    end
end

function run_example()
    % Example usage of peak marker utilities
    
    fprintf('=== Peak Marker Utilities Example ===\n');
    fprintf('This example shows how to use the peak marker utilities.\n\n');
    
    % Load your data first
    fprintf('1. Load your data:\n');
    fprintf('   addpath(genpath(pwd));\n');
    fprintf('   test_file = ''your_file.csv'';\n');
    fprintf('   loader = csv_loader_v2();\n');
    fprintf('   [data, metadata] = loader.loadSingleFile(test_file);\n');
    fprintf('   config = tracenorm_config();\n');
    fprintf('   [baseline, ~, ~] = baseline_detector(data, config);\n');
    fprintf('   [dfof_data, ~] = dfof_calculator(data, baseline, config);\n');
    fprintf('   time_vector = (1:size(data, 1)) / config.frame_rate;\n\n');
    
    fprintf('2. Detect peaks:\n');
    fprintf('   peaks = peak_marker_utils(''detect'', dfof_data, config);\n\n');
    
    fprintf('3. Plot a specific ROI with peaks:\n');
    fprintf('   roi_idx = 100;  %% Choose your ROI\n');
    fprintf('   peak_marker_utils(''plot_roi_with_peaks'', roi_idx, data, baseline, dfof_data, peaks, time_vector);\n\n');
    
    fprintf('4. Or add peaks to your own plots:\n');
    fprintf('   figure; plot(time_vector, dfof_data(:, roi_idx));\n');
    fprintf('   peak_marker_utils(''add_to_plot'', time_vector, dfof_data(:, roi_idx), peaks{roi_idx}, ''dfof'');\n\n');
    
    fprintf('Try running these commands with your data!\n');
end