%% TEST PEAK MARKERS - Quick test of peak detection and marking
% Run this to see peak markers on your spontaneous release data

clear; clc; close all;

fprintf('=== Testing Peak Markers on Spontaneous Release Data ===\n');

%% === Load Data ===
addpath(genpath(pwd));
test_file = 'D:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglusnfr4f_NGR\Spont\GPU_SNR_Processed\5_raw_mean\CP_Snfr4-NGR_Doc2b-R213W_Cs1-c1_spont-01_mean.csv';

fprintf('Loading data and processing...\n');
loader = csv_loader_v2();
[data, metadata] = loader.loadSingleFile(test_file);
config = tracenorm_config();

[baseline, outlier_mask, baseline_stats] = baseline_detector(data, config);
[dfof_data, dfof_stats] = dfof_calculator(data, baseline, config);

time_vector = (1:size(data, 1)) / config.frame_rate;

%% === Detect Peaks ===
fprintf('Detecting peaks in dF/F data...\n');

% Save the peak_marker_utils function first (copy from above)
save_peak_utils();

% Detect peaks
peaks = peak_marker_utils('detect', dfof_data, config);

% Count total peaks
total_peaks = sum(cellfun(@length, peaks));
active_rois = sum(cellfun(@length, peaks) > 0);

fprintf('Peak detection results:\n');
fprintf('  Active ROIs: %d/%d (%.1f%%)\n', active_rois, size(data, 2), 100*active_rois/size(data, 2));
fprintf('  Total peaks: %d\n', total_peaks);
fprintf('  Average peaks per active ROI: %.1f\n', total_peaks/active_rois);

%% === Show Most Active ROIs with Peak Markers ===
fprintf('\nCreating plots with peak markers...\n');

% Find most active ROIs
events_per_roi = cellfun(@length, peaks);
[~, most_active_idx] = sort(events_per_roi, 'descend');

% Show top 8 most active ROIs
num_to_show = 8;
selected_rois = most_active_idx(1:num_to_show);

% Create the plot
fig1 = figure('Name', 'Spontaneous Release with Peak Markers', 'Position', [50, 50, 1600, 1000]);

for i = 1:num_to_show
    roi_idx = selected_rois(i);
    roi_peaks = peaks{roi_idx};
    
    % Raw trace subplot
    subplot(4, 4, 2*i-1);
    plot(time_vector, data(:, roi_idx), 'k-', 'LineWidth', 0.8);
    hold on;
    plot(time_vector, baseline(:, roi_idx), 'r-', 'LineWidth', 1.2);
    
    % Add peak markers to raw trace
    if ~isempty(roi_peaks)
        peak_marker_utils('add_to_plot', time_vector, data(:, roi_idx), roi_peaks, 'raw');
    end
    
    xlabel('Time (s)');
    ylabel('Fluorescence');
    title(sprintf('ROI %d: %d peaks', roi_idx, length(roi_peaks)));
    grid on;
    
    if i == 1
        legend('Raw', 'Baseline', 'Peaks', 'Location', 'best', 'FontSize', 8);
    end
    
    % dF/F trace subplot
    subplot(4, 4, 2*i);
    plot(time_vector, dfof_data(:, roi_idx), 'b-', 'LineWidth', 1.0);
    hold on;
    yline(0, 'k--', 'Alpha', 0.5);
    yline(0.02, 'g:', 'Alpha', 0.7);
    
    % Add peak markers to dF/F trace
    if ~isempty(roi_peaks)
        peak_marker_utils('add_to_plot', time_vector, dfof_data(:, roi_idx), roi_peaks, 'dfof');
    end
    
    xlabel('Time (s)');
    ylabel('dF/F');
    title(sprintf('Max dF/F: %.3f', max(dfof_data(:, roi_idx))));
    grid on;
    
    if i == 1
        legend('dF/F', 'Zero', '2%', 'Peaks', 'Location', 'best', 'FontSize', 8);
    end
end

sgtitle(sprintf('Most Active ROIs with Peak Detection - %s', metadata.filename), 'FontSize', 14);

%% === Show Individual ROI Details ===
% Pick the most active ROI for detailed view
most_active_roi = selected_rois(1);
fprintf('\nShowing detailed view of most active ROI: %d\n', most_active_roi);

peak_marker_utils('plot_roi_with_peaks', most_active_roi, data, baseline, dfof_data, peaks, time_vector);

%% === Peak Statistics ===
fprintf('\n=== Peak Detection Statistics ===\n');
fprintf('Top 10 Most Active ROIs:\n');
fprintf('ROI\tPeaks\tMax dF/F\tPeak Rate (Hz)\n');
fprintf('---\t-----\t--------\t--------------\n');

for i = 1:min(10, length(selected_rois))
    roi_idx = selected_rois(i);
    num_peaks = length(peaks{roi_idx});
    max_dfof = max(dfof_data(:, roi_idx));
    recording_duration = length(time_vector) / config.frame_rate;
    peak_rate = num_peaks / recording_duration;
    
    fprintf('%d\t%d\t%.3f\t\t%.2f\n', roi_idx, num_peaks, max_dfof, peak_rate);
end

%% === Access Peak Data ===
fprintf('\n=== Accessing Peak Data ===\n');
fprintf('Peak locations are stored in cell array ''peaks'':\n');
fprintf('  peaks{roi_idx} contains frame indices of peaks for that ROI\n');
fprintf('  Example: peaks{%d} = [%s]\n', most_active_roi, sprintf('%d ', peaks{most_active_roi}(1:min(5, length(peaks{most_active_roi})))));

if length(peaks{most_active_roi}) > 5
    fprintf('           ... and %d more peaks\n', length(peaks{most_active_roi}) - 5);
end

fprintf('\nTo get peak times in seconds:\n');
fprintf('  peak_times = time_vector(peaks{roi_idx});\n');
fprintf('\nTo get peak amplitudes:\n');
fprintf('  peak_amplitudes = dfof_data(peaks{roi_idx}, roi_idx);\n');

%% === Summary ===
fprintf('\n=== SUMMARY ===\n');
fprintf('✓ Peak detection working!\n');
fprintf('✓ Peak markers added to plots (green triangles on raw, orange on dF/F)\n');
fprintf('✓ Found %d peaks across %d active ROIs\n', total_peaks, active_rois);
fprintf('✓ Most active ROI: %d with %d peaks\n', most_active_roi, length(peaks{most_active_roi}));

fprintf('\nNext steps:\n');
fprintf('1. Review the plotted peaks - do they look like real spontaneous release events?\n');
fprintf('2. Adjust peak detection parameters if needed\n');
fprintf('3. Use the updated baseline_plotter.m for your main pipeline\n');
fprintf('4. Use peak_marker_utils for custom plotting\n');

%% === Helper Function ===
function save_peak_utils()
    % Save the peak_marker_utils function to file if it doesn't exist
    if ~exist('peak_marker_utils.m', 'file')
        fprintf('Creating peak_marker_utils.m function...\n');
        
        % Write the function to file (simplified version)
        fid = fopen('peak_marker_utils.m', 'w');
        if fid > 0
            fprintf(fid, 'function varargout = peak_marker_utils(action, varargin)\n');
            fprintf(fid, '%% Quick peak detection and marking utility\n');
            fprintf(fid, 'switch lower(action)\n');
            fprintf(fid, '    case ''detect''\n');
            fprintf(fid, '        varargout{1} = detect_peaks_simple(varargin{:});\n');
            fprintf(fid, '    case ''add_to_plot''\n');
            fprintf(fid, '        add_peaks_simple(varargin{:});\n');
            fprintf(fid, '    case ''plot_roi_with_peaks''\n');
            fprintf(fid, '        plot_roi_simple(varargin{:});\n');
            fprintf(fid, 'end\n');
            fprintf(fid, 'end\n\n');
            
            fprintf(fid, 'function peaks = detect_peaks_simple(dfof_data, config)\n');
            fprintf(fid, '[~, numROIs] = size(dfof_data);\n');
            fprintf(fid, 'peaks = cell(1, numROIs);\n');
            fprintf(fid, 'min_height = 0.02;\n');
            fprintf(fid, 'min_dist = round(0.5 * config.frame_rate);\n');
            fprintf(fid, 'for roi = 1:numROIs\n');
            fprintf(fid, '    try\n');
            fprintf(fid, '        [~, locs] = findpeaks(dfof_data(:, roi), ''MinPeakHeight'', min_height, ''MinPeakDistance'', min_dist);\n');
            fprintf(fid, '        peaks{roi} = locs;\n');
            fprintf(fid, '    catch\n');
            fprintf(fid, '        peaks{roi} = [];\n');
            fprintf(fid, '    end\n');
            fprintf(fid, 'end\n');
            fprintf(fid, 'end\n\n');
            
            fprintf(fid, 'function add_peaks_simple(time_vec, data, peak_idx, type)\n');
            fprintf(fid, 'if isempty(peak_idx), return; end\n');
            fprintf(fid, 'hold on;\n');
            fprintf(fid, 'if strcmp(type, ''raw'')\n');
            fprintf(fid, '    scatter(time_vec(peak_idx), data(peak_idx), 60, ''^'', ''MarkerFaceColor'', [0 0.8 0], ''MarkerEdgeColor'', [0 0.6 0], ''LineWidth'', 1.5);\n');
            fprintf(fid, 'else\n');
            fprintf(fid, '    scatter(time_vec(peak_idx), data(peak_idx), 60, ''^'', ''MarkerFaceColor'', [1 0.5 0], ''MarkerEdgeColor'', [1 0.3 0], ''LineWidth'', 1.5);\n');
            fprintf(fid, 'end\n');
            fprintf(fid, 'end\n\n');
            
            fprintf(fid, 'function plot_roi_simple(roi_idx, raw_data, baseline, dfof_data, peaks, time_vector)\n');
            fprintf(fid, 'figure(''Name'', sprintf(''ROI %%d with Peaks'', roi_idx));\n');
            fprintf(fid, 'subplot(2,1,1); plot(time_vector, raw_data(:, roi_idx), ''k-''); hold on;\n');
            fprintf(fid, 'plot(time_vector, baseline(:, roi_idx), ''r-'');\n');
            fprintf(fid, 'add_peaks_simple(time_vector, raw_data(:, roi_idx), peaks{roi_idx}, ''raw'');\n');
            fprintf(fid, 'title(sprintf(''ROI %%d Raw (%%d peaks)'', roi_idx, length(peaks{roi_idx})));\n');
            fprintf(fid, 'subplot(2,1,2); plot(time_vector, dfof_data(:, roi_idx), ''b-''); hold on;\n');
            fprintf(fid, 'add_peaks_simple(time_vector, dfof_data(:, roi_idx), peaks{roi_idx}, ''dfof'');\n');
            fprintf(fid, 'title(''dF/F with Peaks'');\n');
            fprintf(fid, 'end\n');
            
            fclose(fid);
            fprintf('✓ Created peak_marker_utils.m\n');
        end
    end
end