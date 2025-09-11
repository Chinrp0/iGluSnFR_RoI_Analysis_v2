%% REGENERATE_ACTIVE_PLOTS - Update plots to show active ROIs only with finer binning
% Run this to regenerate your comparison plots with the improvements:
% 1. Remove 0-event ROIs from all summary plots
% 2. Use finer frequency binning to spread out the low-frequency data

close all; % Close existing figures

fprintf('=== Regenerating Active ROI Plots ===\n');

% Check if you have results from the previous analysis
if ~exist('results', 'var') || ~isfield(results, 'wt_data')
    fprintf('Results not found. Please run the batch analysis first:\n');
    fprintf('  results = fixed_batch_condition_analysis(data_folder);\n');
    return;
end

% Extract data
wt_data = results.wt_data;
mut_data = results.mut_data;
comparison = results.comparison;

fprintf('Regenerating plots with active ROI focus...\n');
fprintf('  WT: %d total ROIs, %d active (%.1f%%)\n', ...
    wt_data.total_rois, wt_data.stats.frequency.active_rois.count, ...
    100 * wt_data.stats.frequency.active_rois.fraction);
fprintf('  R213W: %d total ROIs, %d active (%.1f%%)\n', ...
    mut_data.total_rois, mut_data.stats.frequency.active_rois.count, ...
    100 * mut_data.stats.frequency.active_rois.fraction);

%% Create improved plots
new_plots = struct();

% 1. Enhanced Frequency Comparison (active ROIs only, fine binning)
new_plots.frequency_comparison = create_frequency_comparison_plot(wt_data, mut_data, comparison);

% 2. Enhanced Activity Comparison (no inactive ROIs)
new_plots.activity_comparison = create_activity_comparison_plot(wt_data, mut_data, comparison);

% 3. Keep the amplitude comparison as is (already shows all events)
new_plots.amplitude_comparison = create_amplitude_comparison_plot(wt_data, mut_data, comparison);

% 4. Keep the file summary as is
new_plots.file_summary = create_file_summary_plot(wt_data, mut_data, comparison);

fprintf('✓ Generated %d improved plots\n', length(fieldnames(new_plots)));

%% Summary of changes
fprintf('\n=== Plot Improvements ===\n');
fprintf('1. Frequency Comparison:\n');
fprintf('   - Finer binning for better resolution near 0 Hz\n');
fprintf('   - Clear labeling as "Active ROIs Only"\n');
fprintf('   - Mean frequency lines with labels\n');

fprintf('2. Activity Comparison:\n');
fprintf('   - Completely removed 0-event ROIs\n');
fprintf('   - Updated categories: 1, 2-3, 4-5, 6-10, 11+ events\n');
fprintf('   - Statistics focus on active ROIs only\n');

fprintf('3. All plots clearly labeled as active ROI analysis\n');

% Optional: Save the new plots
save_option = input('Save improved plots to files? (y/n): ', 's');
if lower(save_option) == 'y'
    try
        timestamp = datestr(now, 'yyyymmdd_HHMMSS');
        output_dir = fullfile('active_roi_plots', timestamp);
        if ~exist(output_dir, 'dir')
            mkdir(output_dir);
        end
        
        plot_names = fieldnames(new_plots);
        for i = 1:length(plot_names)
            plot_name = plot_names{i};
            fig_handle = new_plots.(plot_name);
            
            if ishandle(fig_handle)
                output_file = fullfile(output_dir, sprintf('active_rois_%s.png', plot_name));
                saveas(fig_handle, output_file, 'png');
                fprintf('  Saved: %s\n', output_file);
            end
        end
        fprintf('All plots saved to: %s\n', output_dir);
    catch
        fprintf('Failed to save plots\n');
    end
end

fprintf('\n=== Biological Interpretation ===\n');
fprintf('Active ROI Analysis Results:\n');

% Extract active ROI data
wt_active_freqs = wt_data.all_roi_frequencies(wt_data.all_roi_frequencies > 0);
mut_active_freqs = mut_data.all_roi_frequencies(mut_data.all_roi_frequencies > 0);

fprintf('Event Frequency (among active ROIs):\n');
fprintf('  WT: %.4f ± %.4f Hz (n=%d active ROIs)\n', ...
    mean(wt_active_freqs), std(wt_active_freqs), length(wt_active_freqs));
fprintf('  R213W: %.4f ± %.4f Hz (n=%d active ROIs)\n', ...
    mean(mut_active_freqs), std(mut_active_freqs), length(mut_active_freqs));

if length(wt_active_freqs) > 1 && length(mut_active_freqs) > 1
    [p_active, ~] = ranksum(wt_active_freqs, mut_active_freqs);
    fprintf('  Mann-Whitney U test: p = %.4f', p_active);
    if p_active < 0.05
        fprintf(' *SIGNIFICANT*');
    end
    fprintf('\n');
end

fprintf('\nROI Recruitment:\n');
fprintf('  WT: %.1f%% of ROIs are active (%d/%d)\n', ...
    100 * length(wt_active_freqs) / wt_data.total_rois, ...
    length(wt_active_freqs), wt_data.total_rois);
fprintf('  R213W: %.1f%% of ROIs are active (%d/%d)\n', ...
    100 * length(mut_active_freqs) / mut_data.total_rois, ...
    length(mut_active_freqs), mut_data.total_rois);

recruitment_diff = (length(wt_active_freqs) / wt_data.total_rois) - ...
                   (length(mut_active_freqs) / mut_data.total_rois);
fprintf('  R213W recruitment deficit: %.1f percentage points\n', 100 * recruitment_diff);

fprintf('\n=== Key Findings ===\n');
if recruitment_diff > 0.05  % 5% difference threshold
    fprintf('1. REDUCED RECRUITMENT: R213W shows fewer active ROIs\n');
end

if exist('p_active', 'var') && p_active < 0.05
    if mean(wt_active_freqs) > mean(mut_active_freqs)
        fprintf('2. REDUCED ACTIVITY: Active R213W ROIs fire less frequently\n');
    else
        fprintf('2. INCREASED ACTIVITY: Active R213W ROIs fire more frequently\n');
    end
else
    fprintf('2. SIMILAR ACTIVITY: No significant difference in frequency among active ROIs\n');
end

fprintf('\nConclusion: ');
if recruitment_diff > 0.05 && exist('p_active', 'var') && p_active < 0.05
    fprintf('R213W shows BOTH reduced recruitment AND reduced activity per active ROI\n');
elseif recruitment_diff > 0.05
    fprintf('R213W primarily shows reduced recruitment (fewer active ROIs)\n');
elseif exist('p_active', 'var') && p_active < 0.05
    fprintf('R213W primarily affects activity levels among recruited ROIs\n');
else
    fprintf('R213W effects are subtle and require larger sample sizes to detect\n');
end

fprintf('\n=== Analysis Complete ===\n');