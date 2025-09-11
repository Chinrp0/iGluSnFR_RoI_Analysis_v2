function fig = plot_condition_roi_traces(results, options)
    % PLOT_CONDITION_ROI_TRACES - Creates 2x4 subplot figure comparing WT vs R213W traces
    % 
    % Layout:
    % Row 1: Highest frequency ROI    (WT | R213W)
    % Row 2: Lowest frequency ROI     (WT | R213W)  
    % Row 3: Highest amplitude ROI    (WT | R213W)
    % Row 4: Lowest amplitude ROI     (WT | R213W)
    %
    % Args:
    %   results - Output from integrated_batch_analysis_final()
    %   options - Optional struct with plotting preferences
    
    if nargin < 2
        options = struct();
    end
    
    % Set defaults
    if ~isfield(options, 'frame_rate'), options.frame_rate = 40; end  % Hz
    if ~isfield(options, 'event_color'), options.event_color = [1, 0, 0]; end  % Red
    if ~isfield(options, 'trace_color'), options.trace_color = [0, 0, 0]; end  % Black
    if ~isfield(options, 'show_thresholds'), options.show_thresholds = true; end
    
    fprintf('Creating WT vs R213W ROI comparison traces...\n');
    
    %% === Extract and Find Target ROIs ===
    [target_rois, trace_data] = find_target_rois_and_data(results, options);
    
    %% === Create Figure ===
    fig = figure('Name', 'WT vs R213W ROI Comparison Traces', ...
        'Position', [100, 100, 1400, 1200]);
    
    % 4 rows x 2 columns layout
    rows = 4;
    cols = 2;
    
    row_labels = {'Highest Frequency', 'Lowest Frequency', 'Highest Amplitude', 'Lowest Amplitude'};
    col_labels = {'WT', 'R213W'};
    
    %% === Plot Each Comparison ===
    for row = 1:rows
        for col = 1:cols
            subplot_idx = (row - 1) * cols + col;
            subplot(rows, cols, subplot_idx);
            
            % Get data for this subplot
            condition = col_labels{col};
            category = row_labels{row};
            
            if col == 1  % WT
                roi_data = target_rois.wt.(get_category_field(category));
                file_data = trace_data.wt;
            else  % R213W
                roi_data = target_rois.mut.(get_category_field(category));
                file_data = trace_data.mut;
            end
            
            % Plot the trace
            plot_single_roi_trace(roi_data, file_data, options, condition, category);
        end
    end
    
    % Add overall title
    sgtitle('WT vs R213W: Representative ROI Traces with Events', ...
        'FontSize', 16, 'FontWeight', 'bold');
    
    fprintf('  Created 2x4 trace comparison figure\n');
end

function [target_rois, trace_data] = find_target_rois_and_data(results, options)
    % Find the target ROIs and extract their trace data
    
    fprintf('  Selecting target ROIs...\n');
    
    %% === Extract WT Data ===
    wt_data = results.wt_data;
    wt_files = results.file_results.wt_results;
    
    % Remove empty files
    wt_valid = wt_files(~cellfun(@isempty, wt_files));
    
    %% === Extract R213W Data ===
    mut_data = results.mut_data;
    mut_files = results.file_results.mut_results;
    
    % Remove empty files  
    mut_valid = mut_files(~cellfun(@isempty, mut_files));
    
    %% === Find Target ROIs for WT ===
    target_rois.wt = find_condition_target_rois(wt_data, wt_valid);
    
    %% === Find Target ROIs for R213W ===
    target_rois.mut = find_condition_target_rois(mut_data, mut_valid);
    
    %% === Extract Trace Data ===
    trace_data.wt = extract_trace_data_for_targets(target_rois.wt, wt_valid);
    trace_data.mut = extract_trace_data_for_targets(target_rois.mut, mut_valid);
    
    % Print selection summary
    fprintf('    WT targets: freq[%d,%d], amp[%d,%d]\n', ...
        target_rois.wt.highest_freq.global_roi, target_rois.wt.lowest_freq.global_roi, ...
        target_rois.wt.highest_amp.global_roi, target_rois.wt.lowest_amp.global_roi);
    fprintf('    R213W targets: freq[%d,%d], amp[%d,%d]\n', ...
        target_rois.mut.highest_freq.global_roi, target_rois.mut.lowest_freq.global_roi, ...
        target_rois.mut.highest_amp.global_roi, target_rois.mut.lowest_amp.global_roi);
end

function targets = find_condition_target_rois(condition_data, file_results)
    % Find the highest/lowest frequency and amplitude ROIs for one condition
    
    %% === Frequency-based Selection ===
    freqs = condition_data.all_roi_frequencies;
    active_freqs = freqs(freqs > 0);  % Only consider active ROIs
    
    if ~isempty(active_freqs)
        % Find highest frequency ROI
        [max_freq, max_freq_idx] = max(freqs);
        targets.highest_freq = find_roi_in_files(max_freq_idx, file_results);
        targets.highest_freq.value = max_freq;
        
        % Find lowest frequency ROI (among active ROIs)
        [min_active_freq, min_active_global_idx] = min(active_freqs);
        active_indices = find(freqs > 0);
        min_freq_idx = active_indices(min_active_global_idx);
        targets.lowest_freq = find_roi_in_files(min_freq_idx, file_results);
        targets.lowest_freq.value = min_active_freq;
    else
        % No active ROIs - use zeros
        targets.highest_freq = create_empty_target();
        targets.lowest_freq = create_empty_target();
    end
    
    %% === Amplitude-based Selection ===
    all_amplitudes = condition_data.all_individual_events;
    
    if ~isempty(all_amplitudes)
        % Find ROI with highest individual peak
        roi_max_peaks = zeros(1, length(freqs));
        roi_cumulative_idx = 0;
        
        for file_idx = 1:length(file_results)
            if ~isempty(file_results{file_idx}) && isfield(file_results{file_idx}, 'event_stats')
                es = file_results{file_idx}.event_stats;
                file_rois = file_results{file_idx}.numROIs;
                
                if isfield(es, 'mean_peak_amplitude_per_roi')
                    roi_max_peaks((roi_cumulative_idx+1):(roi_cumulative_idx+file_rois)) = ...
                        es.mean_peak_amplitude_per_roi;
                end
                roi_cumulative_idx = roi_cumulative_idx + file_rois;
            end
        end
        
        % Find highest amplitude ROI
        [max_amp, max_amp_idx] = max(roi_max_peaks);
        targets.highest_amp = find_roi_in_files(max_amp_idx, file_results);
        targets.highest_amp.value = max_amp;
        
        % Find lowest amplitude ROI (among those with events)
        active_amps = roi_max_peaks(roi_max_peaks > 0);
        if ~isempty(active_amps)
            [min_amp, min_amp_local_idx] = min(active_amps);
            active_amp_indices = find(roi_max_peaks > 0);
            min_amp_idx = active_amp_indices(min_amp_local_idx);
            targets.lowest_amp = find_roi_in_files(min_amp_idx, file_results);
            targets.lowest_amp.value = min_amp;
        else
            targets.lowest_amp = create_empty_target();
        end
    else
        targets.highest_amp = create_empty_target();
        targets.lowest_amp = create_empty_target();
    end
end

function roi_info = find_roi_in_files(global_roi_idx, file_results)
    % Find which file and local ROI index corresponds to global ROI index
    
    roi_cumulative = 0;
    
    for file_idx = 1:length(file_results)
        if ~isempty(file_results{file_idx})
            file_rois = file_results{file_idx}.numROIs;
            
            if global_roi_idx <= (roi_cumulative + file_rois)
                % Found the file!
                local_roi_idx = global_roi_idx - roi_cumulative;
                
                roi_info = struct();
                roi_info.file_idx = file_idx;
                roi_info.local_roi = local_roi_idx;
                roi_info.global_roi = global_roi_idx;
                roi_info.filename = file_results{file_idx}.filename;
                return;
            end
            
            roi_cumulative = roi_cumulative + file_rois;
        end
    end
    
    % If not found, create empty target
    roi_info = create_empty_target();
end

function target = create_empty_target()
    % Create empty target for missing data
    target = struct('file_idx', 0, 'local_roi', 0, 'global_roi', 0, ...
                   'filename', 'N/A', 'value', 0);
end

function trace_data = extract_trace_data_for_targets(targets, file_results)
    % Extract the actual dF/F traces and event masks for target ROIs
    
    categories = {'highest_freq', 'lowest_freq', 'highest_amp', 'lowest_amp'};
    trace_data = struct();
    
    for i = 1:length(categories)
        category = categories{i};
        target = targets.(category);
        
        if target.file_idx > 0 && target.file_idx <= length(file_results)
            file_result = file_results{target.file_idx};
            
            if ~isempty(file_result) && isfield(file_result, 'dfof_data')
                % Extract trace and event mask
                trace_data.(category).trace = file_result.dfof_data(:, target.local_roi);
                
                if isfield(file_result, 'event_mask')
                    trace_data.(category).event_mask = file_result.event_mask(:, target.local_roi);
                else
                    trace_data.(category).event_mask = false(size(trace_data.(category).trace));
                end
                
                % Get thresholds if available
                if isfield(file_result, 'event_stats') && isfield(file_result.event_stats, 'upper_thresholds')
                    trace_data.(category).upper_thresh = file_result.event_stats.upper_thresholds(target.local_roi);
                    trace_data.(category).lower_thresh = file_result.event_stats.lower_thresholds(target.local_roi);
                else
                    trace_data.(category).upper_thresh = NaN;
                    trace_data.(category).lower_thresh = NaN;
                end
                
                trace_data.(category).target_info = target;
            else
                trace_data.(category) = create_empty_trace_data();
            end
        else
            trace_data.(category) = create_empty_trace_data();
        end
    end
end

function empty_data = create_empty_trace_data()
    % Create empty trace data structure
    empty_data = struct();
    empty_data.trace = NaN(1200, 1);  % Default size
    empty_data.event_mask = false(1200, 1);
    empty_data.upper_thresh = NaN;
    empty_data.lower_thresh = NaN;
    empty_data.target_info = create_empty_target();
end

function plot_single_roi_trace(roi_data, file_data, options, condition, category)
    % Plot a single dF/F trace with event highlighting
    
    category_field = get_category_field(category);
    
    if ~isfield(file_data, category_field) || isempty(file_data.(category_field).trace)
        % No data available
        text(0.5, 0.5, 'No data', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'FontSize', 12);
        title(sprintf('%s\n%s', condition, category), 'FontSize', 11);
        return;
    end
    
    trace_info = file_data.(category_field);
    trace = trace_info.trace;
    events = trace_info.event_mask;
    
    % Create time vector
    numFrames = length(trace);
    time_vector = (1:numFrames) / options.frame_rate;
    
    % Plot main trace
    plot(time_vector, trace, 'Color', options.trace_color, 'LineWidth', 1.0);
    hold on;
    
    % Highlight events in red
    if any(events)
        event_trace = trace;
        event_trace(~events) = NaN;
        plot(time_vector, event_trace, 'Color', options.event_color, 'LineWidth', 2.5);
    end
    
    % Add thresholds if available and requested
    if options.show_thresholds && ~isnan(trace_info.upper_thresh)
        yline(trace_info.upper_thresh, 'g--', 'LineWidth', 1, 'Alpha', 0.7);
        yline(trace_info.lower_thresh, 'c--', 'LineWidth', 1, 'Alpha', 0.7);
    end
    
    % Add zero reference
    yline(0, 'k:', 'Alpha', 0.5, 'LineWidth', 0.8);
    
    % Formatting
    xlabel('Time (s)', 'FontSize', 10);
    ylabel('dF/F', 'FontSize', 10);
    
    % Create informative title
    target_info = trace_info.target_info;
    if target_info.file_idx > 0
        title_str = sprintf('%s\n%s\nROI %d (%.4f)', condition, category, ...
            target_info.global_roi, target_info.value);
    else
        title_str = sprintf('%s\n%s\nNo data', condition, category);
    end
    
    title(title_str, 'FontSize', 11, 'FontWeight', 'bold');
    grid on;
    
    % Set y-limits with some margin
    if ~all(isnan(trace))
        y_range = [min(trace, [], 'omitnan'), max(trace, [], 'omitnan')];
        if diff(y_range) > 0
            y_margin = diff(y_range) * 0.15;
            ylim([y_range(1) - y_margin, y_range(2) + y_margin]);
        end
    end
end

function field_name = get_category_field(category)
    % Convert category label to field name
    switch category
        case 'Highest Frequency'
            field_name = 'highest_freq';
        case 'Lowest Frequency'
            field_name = 'lowest_freq';
        case 'Highest Amplitude'
            field_name = 'highest_amp';
        case 'Lowest Amplitude'
            field_name = 'lowest_amp';
        otherwise
            error('Unknown category: %s', category);
    end
end