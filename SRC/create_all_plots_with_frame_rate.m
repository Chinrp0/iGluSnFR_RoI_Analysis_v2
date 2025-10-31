function create_all_plots_with_frame_rate(results)
    % CREATE_ALL_PLOTS_WITH_FRAME_RATE - Properly generate all visualization modules
    %
    % This function calls all plotting modules with the correct frame_rate parameter
    % to fix the "Unrecognized field name 'frame_rate'" errors
    %
    % Usage:
    %   results = integrated_batch_analysis_final(folder_path);
    %   create_all_plots_with_frame_rate(results);
    %
    % This will create:
    %   1. Inter-Event Interval comparison
    %   2. Amplitude-Frequency correlation
    %   3. Event timing raster
    %   4. Cumulative events over time
    %   5. ROI trace comparison
    
    fprintf('\n=== Creating Additional Plots with Frame Rate ===\n');
    
    % Set up options with frame_rate from results
    options = struct();
    
    if isfield(results, 'processing_info') && isfield(results.processing_info, 'frame_rate')
        options.frame_rate = results.processing_info.frame_rate;
        fprintf('Using frame_rate: %d Hz\n', options.frame_rate);
    else
        options.frame_rate = 100;  % Default to 100 Hz (10ms exposure)
        fprintf('WARNING: Using default frame_rate: %d Hz\n', options.frame_rate);
    end
    
    % Also pass recording duration if available
    if isfield(results, 'processing_info') && isfield(results.processing_info, 'recording_duration_s')
        options.recording_duration_s = results.processing_info.recording_duration_s;
    end
    
    plot_count = 0;
    
    %% 1. Inter-Event Interval Comparison
    try
        fprintf('\n1. Creating Inter-Event Interval comparison...\n');
        fig = plot_iei_comparison(results, options);
        plot_count = plot_count + 1;
        fprintf('   ✓ IEI comparison created (Figure %d)\n', fig.Number);
    catch ME
        fprintf('   ✗ Failed to create IEI plot: %s\n', ME.message);
    end
    
    %% 2. Amplitude-Frequency Correlation
    try
        fprintf('\n2. Creating Amplitude-Frequency correlation...\n');
        fig = plot_amplitude_frequency_correlation(results, options);
        plot_count = plot_count + 1;
        fprintf('   ✓ Amp-Freq correlation created (Figure %d)\n', fig.Number);
    catch ME
        fprintf('   ✗ Failed to create Amp-Freq correlation: %s\n', ME.message);
    end
    
    %% 3. Event Timing Raster
    try
        fprintf('\n3. Creating Event Timing raster...\n');
        fig = plot_event_raster(results, options);
        plot_count = plot_count + 1;
        fprintf('   ✓ Event raster created (Figure %d)\n', fig.Number);
    catch ME
        fprintf('   ✗ Failed to create event raster: %s\n', ME.message);
    end
    
    %% 4. Cumulative Events Over Time
    try
        fprintf('\n4. Creating Cumulative Events plot...\n');
        fig = plot_cumulative_events(results, options);
        plot_count = plot_count + 1;
        fprintf('   ✓ Cumulative events created (Figure %d)\n', fig.Number);
    catch ME
        fprintf('   ✗ Failed to create cumulative events: %s\n', ME.message);
    end
    
    %% 5. ROI Trace Comparison
    try
        fprintf('\n5. Creating ROI Trace comparison...\n');
        fig = plot_condition_roi_traces(results, options);
        plot_count = plot_count + 1;
        fprintf('   ✓ ROI traces created (Figure %d)\n', fig.Number);
    catch ME
        fprintf('   ✗ Failed to create ROI traces: %s\n', ME.message);
    end
    
    fprintf('\n=== Summary ===\n');
    fprintf('Successfully created %d additional plots\n', plot_count);
    fprintf('Total figures now open: %d\n', length(findall(0, 'Type', 'figure')));
    fprintf('\nTo list all figures: list_open_figures()\n');
    fprintf('To close all figures: close all\n\n');
end
