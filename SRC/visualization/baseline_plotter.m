function plot_handles = baseline_plotter(raw_data, baseline, dfof_data, outlier_mask, stats, metadata, config)
    % BASELINE_PLOTTER - Visualize baseline calculation performance
    % Creates comprehensive plots to validate baseline detection quality
    %
    % Inputs:
    %   raw_data     - [frames x ROIs] original fluorescence data
    %   baseline     - [frames x ROIs] calculated baseline
    %   dfof_data    - [frames x ROIs] normalized dF/F traces  
    %   outlier_mask - [frames x ROIs] logical mask of detected outliers
    %   stats        - Statistics struct from baseline_detector
    %   metadata     - File metadata struct
    %   config       - Configuration struct
    %
    % Outputs:
    %   plot_handles - Struct containing figure handles for saving
    
    if nargin < 7
        config = tracenorm_config();
    end
    
    [numFrames, numROIs] = size(raw_data);
    time_vector = (1:numFrames) / config.frame_rate;  % Convert to seconds
    
    plot_handles = struct();
    
    fprintf('Creating baseline validation plots for %s...\n', metadata.filename);
    
    %% === Plot 1: Sample Trace Overview ===
    if config.plot_sample_traces
        plot_handles.trace_overview = create_trace_overview(raw_data, baseline, dfof_data, ...
            outlier_mask, time_vector, metadata, config);
    end
    
    %% === Plot 2: Baseline Quality Validation ===
    if config.plot_validation
        plot_handles.baseline_validation = create_baseline_validation(raw_data, baseline, ...
            stats, metadata, config);
    end
    
    %% === Plot 3: dF/F Quality Assessment ===
    plot_handles.dfof_assessment = create_dfof_assessment(dfof_data, stats, ...
        metadata, config);
    
    %% === Plot 4: Outlier Detection Summary ===
    plot_handles.outlier_summary = create_outlier_summary(outlier_mask, stats, ...
        time_vector, metadata, config);
    
    %% === Plot 5: Transport ROI Detection (if enabled) ===
    if config.detect_transport_rois && stats.num_transport_rois > 0
        plot_handles.transport_detection = create_transport_plots(raw_data, baseline, ...
            stats, time_vector, metadata, config);
    end
    
    fprintf('  Created %d validation plots\n', length(fieldnames(plot_handles)));
end

function fig = create_trace_overview(raw_data, baseline, dfof_data, outlier_mask, time_vector, metadata, config)
    % Show sample traces with baseline overlay and dF/F results
    
    [~, numROIs] = size(raw_data);
    
    % Select representative ROIs for plotting
    sample_indices = select_sample_rois(raw_data, config.num_sample_traces);
    
    fig = figure('Name', sprintf('Baseline Overview - %s', metadata.filename), ...
        'Position', [100, 100, 1200, 800]);
    
    num_plots = length(sample_indices);
    rows = ceil(sqrt(num_plots));
    cols = ceil(num_plots / rows);
    
    for i = 1:num_plots
        roi_idx = sample_indices(i);
        
        % Raw trace + baseline subplot
        subplot(rows, cols*2, 2*i-1);
        
        % Plot raw data
        plot(time_vector, raw_data(:, roi_idx), 'k-', 'LineWidth', 0.8, 'DisplayName', 'Raw');
        hold on;
        
        % Plot baseline
        plot(time_vector, baseline(:, roi_idx), 'r-', 'LineWidth', 1.2, 'DisplayName', 'Baseline');
        
        % Highlight outliers
        outlier_points = outlier_mask(:, roi_idx);
        if any(outlier_points)
            scatter(time_vector(outlier_points), raw_data(outlier_points, roi_idx), ...
                30, 'o', 'MarkerFaceColor', 'orange', 'MarkerEdgeColor', 'orange', ...
                'DisplayName', 'Outliers');
        end
        
        xlabel('Time (s)');
        ylabel('Fluorescence');
        title(sprintf('ROI %d - Raw + Baseline', roi_idx));
        legend('Location', 'best');
        grid on;
        
        % dF/F subplot
        subplot(rows, cols*2, 2*i);
        plot(time_vector, dfof_data(:, roi_idx), 'b-', 'LineWidth', 1.0);
        xlabel('Time (s)');
        ylabel('dF/F');
        title(sprintf('ROI %d - Normalized', roi_idx));
        grid on;
        
        % Add zero line for reference
        yline(0, 'k--', 'Alpha', 0.5);
    end
    
    sgtitle(sprintf('Baseline Calculation Overview - %s', metadata.filename), ...
        'FontSize', 14, 'FontWeight', 'bold');
end

function fig = create_baseline_validation(raw_data, baseline, stats, metadata, config)
    % Validate baseline calculation quality
    
    fig = figure('Name', sprintf('Baseline Validation - %s', metadata.filename), ...
        'Position', [200, 100, 1000, 700]);
    
    % Subplot 1: Baseline stability (CV distribution)
    subplot(2, 3, 1);
    histogram(stats.baseline_cv, 50, 'EdgeColor', 'none', 'FaceColor', [0.3, 0.7, 0.9]);
    xlabel('Baseline CV');
    ylabel('Count');
    title('Baseline Stability');
    xline(0.2, 'r--', 'High Noise Threshold');
    grid on;
    
    % Subplot 2: Signal-to-baseline ratio
    subplot(2, 3, 2);
    histogram(stats.signal_to_baseline, 50, 'EdgeColor', 'none', 'FaceColor', [0.7, 0.9, 0.3]);
    xlabel('Signal/Baseline Ratio');
    ylabel('Count');
    title('Signal Strength');
    grid on;
    
    % Subplot 3: Valid data fraction
    subplot(2, 3, 3);
    histogram(stats.valid_fraction, 50, 'EdgeColor', 'none', 'FaceColor', [0.9, 0.7, 0.3]);
    xlabel('Fraction Valid Data');
    ylabel('Count');
    title('Data Validity');
    xline(0.8, 'r--', 'Minimum Threshold');
    grid on;
    
    % Subplot 4: Outlier fraction per ROI
    subplot(2, 3, 4);
    histogram(stats.outlier_fraction, 50, 'EdgeColor', 'none', 'FaceColor', [0.9, 0.3, 0.7]);
    xlabel('Outlier Fraction');
    ylabel('Count');
    title('Outlier Detection');
    grid on;
    
    % Subplot 5: Transport slope distribution
    subplot(2, 3, 5);
    histogram(stats.transport_slopes, 50, 'EdgeColor', 'none', 'FaceColor', [0.7, 0.3, 0.9]);
    xlabel('Baseline Slope');
    ylabel('Count');
    title('Transport Detection');
    xline(config.transport_slope_threshold, 'r--', 'Transport Threshold');
    grid on;
    
    % Subplot 6: Quality summary
    subplot(2, 3, 6);
    quality_categories = {'Good ROIs', 'Low Valid Data', 'High Noise', 'Transport'};
    quality_counts = [
        sum(~(stats.quality_flags.low_valid_data | stats.quality_flags.high_noise | stats.quality_flags.potential_transport)),
        sum(stats.quality_flags.low_valid_data),
        sum(stats.quality_flags.high_noise), 
        sum(stats.quality_flags.potential_transport)
    ];
    
    pie(quality_counts, quality_categories);
    title('ROI Quality Summary');
    
    sgtitle(sprintf('Baseline Validation - %s (%.1f%% Good ROIs)', ...
        metadata.filename, 100 * stats.fraction_good_rois), ...
        'FontSize', 14, 'FontWeight', 'bold');
end

function fig = create_dfof_assessment(dfof_data, stats, metadata, config)
    % Assess dF/F calculation quality
    
    fig = figure('Name', sprintf('dF/F Assessment - %s', metadata.filename), ...
        'Position', [300, 100, 1000, 600]);
    
    % Subplot 1: dF/F distribution
    subplot(2, 3, 1);
    histogram(dfof_data(:), 100, 'EdgeColor', 'none', 'FaceColor', [0.3, 0.7, 0.9]);
    xlabel('dF/F');
    ylabel('Count');
    title('dF/F Distribution');
    grid on;
    
    % Subplot 2: SNR distribution
    subplot(2, 3, 2);
    histogram(stats.snr, 50, 'EdgeColor', 'none', 'FaceColor', [0.7, 0.9, 0.3]);
    xlabel('Signal-to-Noise Ratio');
    ylabel('Count');
    title('SNR Distribution');
    xline(2, 'r--', 'Minimum SNR');
    grid on;
    
    % Subplot 3: Dynamic range
    subplot(2, 3, 3);
    histogram(stats.dynamic_range, 50, 'EdgeColor', 'none', 'FaceColor', [0.9, 0.7, 0.3]);
    xlabel('Dynamic Range');
    ylabel('Count');
    title('Signal Dynamic Range');
    grid on;
    
    % Subplot 4: Events per ROI
    subplot(2, 3, 4);
    histogram(stats.events_per_roi, 0:max(stats.events_per_roi), ...
        'EdgeColor', 'none', 'FaceColor', [0.9, 0.3, 0.7]);
    xlabel('Events per ROI');
    ylabel('Count');
    title('Event Detection');
    grid on;
    
    % Subplot 5: Max dF/F per ROI
    subplot(2, 3, 5);
    histogram(stats.dfof_max, 50, 'EdgeColor', 'none', 'FaceColor', [0.7, 0.3, 0.9]);
    xlabel('Max dF/F');
    ylabel('Count');
    title('Peak Response Strength');
    grid on;
    
    % Subplot 6: Summary statistics text
    subplot(2, 3, 6);
    axis off;
    
    summary_text = {
        sprintf('Total ROIs: %d', length(stats.dfof_mean)),
        sprintf('Valid Data: %.1f%%', 100 * stats.fraction_valid),
        sprintf('Mean SNR: %.2f', stats.mean_snr),
        sprintf('Active ROIs: %d (%.1f%%)', stats.event_summary.rois_with_events, ...
            100 * stats.event_summary.fraction_active_rois),
        sprintf('Total Events: %d', stats.event_summary.total_events),
        sprintf('Events/Active ROI: %.1f', stats.event_summary.mean_events_per_active_roi),
        '',
        sprintf('Good Quality ROIs: %.1f%%', 100 * stats.fraction_good_rois)
    };
    
    text(0.1, 0.9, summary_text, 'FontSize', 11, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', 'Units', 'normalized');
    
    sgtitle(sprintf('dF/F Quality Assessment - %s', metadata.filename), ...
        'FontSize', 14, 'FontWeight', 'bold');
end

function fig = create_outlier_summary(outlier_mask, stats, time_vector, metadata, config)
    % Summarize outlier detection across time and ROIs
    
    fig = figure('Name', sprintf('Outlier Summary - %s', metadata.filename), ...
        'Position', [400, 100, 1000, 600]);
    
    % Subplot 1: Outliers over time
    subplot(2, 2, 1);
    outliers_per_frame = sum(outlier_mask, 2);
    plot(time_vector, outliers_per_frame, 'r-', 'LineWidth', 1.2);
    xlabel('Time (s)');
    ylabel('Outliers per Frame');
    title('Outlier Detection Timeline');
    grid on;
    
    % Subplot 2: Outliers per ROI
    subplot(2, 2, 2);
    outliers_per_roi = sum(outlier_mask, 1);
    histogram(outliers_per_roi, 50, 'EdgeColor', 'none', 'FaceColor', [0.9, 0.3, 0.3]);
    xlabel('Outliers per ROI');
    ylabel('Count');
    title('Outlier Distribution');
    grid on;
    
    % Subplot 3: Outlier heatmap (subsampled if too many ROIs)
    subplot(2, 2, 3);
    max_rois_to_plot = 100;
    if size(outlier_mask, 2) > max_rois_to_plot
        roi_indices = round(linspace(1, size(outlier_mask, 2), max_rois_to_plot));
        outlier_subset = outlier_mask(:, roi_indices);
    else
        outlier_subset = outlier_mask;
        roi_indices = 1:size(outlier_mask, 2);
    end
    
    imagesc(roi_indices, time_vector, double(outlier_subset));
    colormap(gca, [1 1 1; 1 0 0]);  % White = no outlier, Red = outlier
    xlabel('ROI Index');
    ylabel('Time (s)');
    title('Outlier Pattern');
    colorbar('Ticks', [0, 1], 'TickLabels', {'Normal', 'Outlier'});
    
    % Subplot 4: Summary statistics
    subplot(2, 2, 4);
    axis off;
    
    summary_text = {
        sprintf('Total Outliers: %d (%.2f%%)', sum(outlier_mask(:)), ...
            100 * stats.mean_outlier_fraction),
        sprintf('Frames with Outliers: %d/%d', sum(any(outlier_mask, 2)), length(time_vector)),
        sprintf('ROIs with Outliers: %d/%d', sum(any(outlier_mask, 1)), size(outlier_mask, 2)),
        sprintf('Peak Outliers: %d (frame %d)', max(outliers_per_frame), ...
            find(outliers_per_frame == max(outliers_per_frame), 1)),
        '',
        sprintf('Detection Threshold: %.1fσ', config.outlier_threshold_sigma),
        sprintf('Rolling Window: %.1fs (%d frames)', config.rolling_window_sec, ...
            config.rolling_window_frames)
    };
    
    text(0.1, 0.9, summary_text, 'FontSize', 11, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', 'Units', 'normalized');
    
    sgtitle(sprintf('Outlier Detection Summary - %s', metadata.filename), ...
        'FontSize', 14, 'FontWeight', 'bold');
end

function fig = create_transport_plots(raw_data, baseline, stats, time_vector, metadata, config)
    % Show transport ROI detection results
    
    transport_indices = find(stats.transport_rois);
    if isempty(transport_indices)
        fig = [];
        return;
    end
    
    fig = figure('Name', sprintf('Transport ROIs - %s', metadata.filename), ...
        'Position', [500, 100, 1000, 600]);
    
    num_transport = min(6, length(transport_indices));  % Show up to 6 transport ROIs
    
    for i = 1:num_transport
        roi_idx = transport_indices(i);
        
        subplot(2, 3, i);
        
        % Plot raw data and baseline
        plot(time_vector, raw_data(:, roi_idx), 'k-', 'LineWidth', 0.8, 'DisplayName', 'Raw');
        hold on;
        plot(time_vector, baseline(:, roi_idx), 'r-', 'LineWidth', 1.2, 'DisplayName', 'Baseline');
        
        % Add linear trend line
        p = polyfit(time_vector, baseline(:, roi_idx), 1);
        trend_line = polyval(p, time_vector);
        plot(time_vector, trend_line, 'g--', 'LineWidth', 1.0, 'DisplayName', 'Trend');
        
        xlabel('Time (s)');
        ylabel('Fluorescence');
        title(sprintf('Transport ROI %d (slope=%.3f)', roi_idx, stats.transport_slopes(roi_idx)));
        legend('Location', 'best');
        grid on;
    end
    
    sgtitle(sprintf('Transport ROI Detection - %s (%d found)', ...
        metadata.filename, stats.num_transport_rois), ...
        'FontSize', 14, 'FontWeight', 'bold');
end

function sample_indices = select_sample_rois(data, num_samples)
    % Select representative ROIs for plotting
    
    [~, numROIs] = size(data);
    
    if numROIs <= num_samples
        sample_indices = 1:numROIs;
        return;
    end
    
    % Select ROIs with different characteristics
    roi_variance = var(data, 0, 1);
    roi_mean = mean(data, 1);
    roi_max = max(data, [], 1);
    
    % Normalize metrics for selection
    variance_norm = (roi_variance - min(roi_variance)) / (max(roi_variance) - min(roi_variance));
    mean_norm = (roi_mean - min(roi_mean)) / (max(roi_mean) - min(roi_mean));
    max_norm = (roi_max - min(roi_max)) / (max(roi_max) - min(roi_max));
    
    % Combined score (favor diverse characteristics)
    combined_score = variance_norm + mean_norm + max_norm;
    
    % Select ROIs with different score ranges
    [~, sorted_indices] = sort(combined_score);
    step = floor(numROIs / num_samples);
    sample_indices = sorted_indices(1:step:end);
    sample_indices = sample_indices(1:num_samples);  % Ensure exact number
end