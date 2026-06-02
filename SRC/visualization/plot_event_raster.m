function fig = plot_event_raster(results, options)
    % PLOT_EVENT_RASTER - Event timing raster plots
    % Shows temporal distribution of events across ROIs
    %
    % Args:
    %   results - Output from integrated_batch_analysis_final()
    %   options - Optional struct with:
    %             .frame_rate (default: 100 Hz)
    %             .max_rois_per_plot (default: 100, subsample if more)
    %             .recording_duration_s (default: 30s)
    
    if nargin < 2
        options = struct();
    end
    
    % Set defaults
    if ~isfield(options, 'frame_rate'), options.frame_rate = 100; end
    if ~isfield(options, 'max_rois_per_plot'), options.max_rois_per_plot = 100; end
    if ~isfield(options, 'recording_duration_s'), options.recording_duration_s = 30; end
    
    fprintf('Creating event timing raster plots...\n');
    
    %% Extract event timing data
    [wt_raster, wt_stats] = extract_raster_data(results.file_results.wt_results, ...
        options.frame_rate, options.max_rois_per_plot);
    [mut_raster, mut_stats] = extract_raster_data(results.file_results.mut_results, ...
        options.frame_rate, options.max_rois_per_plot);
    
    fprintf('  WT: %d active ROIs (showing %d)\n', wt_stats.total_active_rois, wt_stats.displayed_rois);
    fprintf('  R213W: %d active ROIs (showing %d)\n', mut_stats.total_active_rois, mut_stats.displayed_rois);
    
    %% Create figure
    fig = figure('Name', 'Event Timing Raster: WT vs R213W', ...
        'Position', [100, 100, 1400, 1000]);
    
    %% Subplot 1: WT Raster
    subplot(3, 2, [1, 3]);
    plot_raster_subplot(wt_raster, 'WT', [0.2, 0.6, 1.0], options.recording_duration_s);
    
    %% Subplot 2: R213W Raster
    subplot(3, 2, [2, 4]);
    plot_raster_subplot(mut_raster, 'R213W', [1.0, 0.4, 0.2], options.recording_duration_s);
    
    %% Subplot 3: Population event rate over time (WT)
    subplot(3, 2, 5);
    plot_population_rate(wt_raster, wt_stats, 'WT', [0.2, 0.6, 1.0], options);
    
    %% Subplot 4: Population event rate over time (R213W)
    subplot(3, 2, 6);
    plot_population_rate(mut_raster, mut_stats, 'R213W', [1.0, 0.4, 0.2], options);
    
    sgtitle('Event Timing Raster Plots: WT vs R213W', 'FontSize', 16, 'FontWeight', 'bold');
    
    fprintf('  Created event timing raster plots\n');
end

function plot_raster_subplot(raster_data, condition_name, color, recording_duration_s)
    % Plot individual raster
    
    if isempty(raster_data)
        text(0.5, 0.5, 'No event data available', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'FontSize', 12);
        title(sprintf('%s Event Raster', condition_name));
        return;
    end
    
    hold on;
    
    % Plot each ROI's events
    for roi_idx = 1:length(raster_data)
        if ~isempty(raster_data{roi_idx})
            event_times = raster_data{roi_idx};
            y_pos = roi_idx * ones(size(event_times));
            plot(event_times, y_pos, '|', 'Color', color, 'MarkerSize', 8, 'LineWidth', 1);
        end
    end
    
    xlabel('Time (s)');
    ylabel('ROI Index');
    title(sprintf('%s Event Raster (n=%d ROIs)', condition_name, length(raster_data)));
    xlim([0, recording_duration_s]);
    ylim([0, length(raster_data) + 1]);
    grid on;
    box on;
end

function plot_population_rate(raster_data, stats, condition_name, color, options)
    % Plot population event rate over time
    
    if isempty(raster_data)
        text(0.5, 0.5, 'No event data available', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'FontSize', 12);
        title(sprintf('%s Population Rate', condition_name));
        return;
    end
    
    % Bin events into time windows
    bin_size_s = 1.0;  % 1 second bins
    time_edges = 0:bin_size_s:options.recording_duration_s;
    event_counts = zeros(1, length(time_edges) - 1);
    
    % Count events in each time bin
    for roi_idx = 1:length(raster_data)
        if ~isempty(raster_data{roi_idx})
            event_times = raster_data{roi_idx};
            counts = histcounts(event_times, time_edges);
            event_counts = event_counts + counts;
        end
    end
    
    % Convert to rate (events/s)
    event_rate = event_counts / bin_size_s;
    time_centers = time_edges(1:end-1) + bin_size_s / 2;
    
    % Plot as area
    area(time_centers, event_rate, 'FaceColor', color, 'FaceAlpha', 0.5, 'EdgeColor', color, 'LineWidth', 2);
    hold on;
    
    % Add mean line
    mean_rate = mean(event_rate);
    yline(mean_rate, '--', sprintf('Mean: %.1f events/s', mean_rate), ...
        'LineWidth', 2, 'Color', 'k', 'LabelHorizontalAlignment', 'left');
    
    xlabel('Time (s)');
    ylabel('Population Event Rate (events/s)');
    title(sprintf('%s Population Activity', condition_name));
    xlim([0, options.recording_duration_s]);
    grid on;
    
    % Add variability info
    cv_rate = std(event_rate) / mean(event_rate);
    text(0.98, 0.98, sprintf('CV = %.3f', cv_rate), 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', ...
        'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize', 9);
end

function [raster_data, stats] = extract_raster_data(file_results, frame_rate, max_rois)
    % Extract event timing data for raster plot
    
    % First pass: collect all active ROIs across all files
    all_active_rois = {};
    all_roi_event_counts = [];
    
    for file_idx = 1:length(file_results)
        if isempty(file_results{file_idx})
            continue;
        end
        
        file_result = file_results{file_idx};
        
        if ~isfield(file_result, 'event_mask') || isempty(file_result.event_mask)
            continue;
        end
        
        event_mask = file_result.event_mask;
        [num_frames, num_rois] = size(event_mask);
        
        % Process each ROI
        for roi = 1:num_rois
            roi_events = event_mask(:, roi);
            
            % Find event start times
            event_starts = find(diff([false; roi_events]) == 1);
            
            if ~isempty(event_starts)
                % Convert to seconds
                event_times_s = event_starts / frame_rate;
                all_active_rois{end+1} = event_times_s;
                all_roi_event_counts(end+1) = length(event_times_s);
            end
        end
    end
    
    total_active_rois = length(all_active_rois);
    
    % Subsample if too many ROIs
    if total_active_rois > max_rois
        % Sort by event count and sample across the range
        [~, sorted_idx] = sort(all_roi_event_counts, 'descend');
        
        % Sample evenly across sorted ROIs
        sample_idx = round(linspace(1, total_active_rois, max_rois));
        display_idx = sorted_idx(sample_idx);
        
        raster_data = all_active_rois(display_idx);
        displayed_rois = length(raster_data);
    else
        raster_data = all_active_rois;
        displayed_rois = total_active_rois;
    end
    
    stats = struct();
    stats.total_active_rois = total_active_rois;
    stats.displayed_rois = displayed_rois;
end
