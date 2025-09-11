function comparison = compare_conditions_statistical(wt_data, mut_data)
    % COMPARE_CONDITIONS_STATISTICAL - Statistical comparison between WT and R213W
    %
    % Performs appropriate statistical tests to compare:
    %   - Event frequencies between conditions
    %   - Event amplitudes between conditions  
    %   - ROI activity fractions between conditions
    %
    % Uses non-parametric tests (Mann-Whitney U) due to likely non-normal distributions
    
    comparison = struct();
    comparison.test_type = 'Mann-Whitney U (ranksum)';
    comparison.alpha = 0.05;
    
    fprintf('  Running statistical comparisons (Mann-Whitney U tests)...\n');
    
    %% === Event Frequency Comparison (Active ROIs only) ===
    wt_freqs = wt_data.all_roi_frequencies(wt_data.all_roi_frequencies > 0);
    mut_freqs = mut_data.all_roi_frequencies(mut_data.all_roi_frequencies > 0);
    
    if ~isempty(wt_freqs) && ~isempty(mut_freqs)
        [comparison.frequency_pvalue, comparison.frequency_h] = ranksum(wt_freqs, mut_freqs);
        comparison.frequency_effect_size = calculate_effect_size(wt_freqs, mut_freqs);
        
        fprintf('    Frequency: WT=%.4f Hz, R213W=%.4f Hz, p=%.4f\n', ...
            mean(wt_freqs), mean(mut_freqs), comparison.frequency_pvalue);
    else
        comparison.frequency_pvalue = NaN;
        comparison.frequency_h = 0;
        comparison.frequency_effect_size = 0;
        fprintf('    Frequency: Insufficient data for comparison\n');
    end
    
    %% === Event Amplitude Comparison (All Individual Events) ===
    wt_amps = wt_data.all_individual_events;
    mut_amps = mut_data.all_individual_events;
    
    if ~isempty(wt_amps) && ~isempty(mut_amps)
        [comparison.amplitude_pvalue, comparison.amplitude_h] = ranksum(wt_amps, mut_amps);
        comparison.amplitude_effect_size = calculate_effect_size(wt_amps, mut_amps);
        
        fprintf('    Amplitude: WT=%.4f dF/F, R213W=%.4f dF/F, p=%.4f\n', ...
            mean(wt_amps), mean(mut_amps), comparison.amplitude_pvalue);
    else
        comparison.amplitude_pvalue = NaN;
        comparison.amplitude_h = 0;
        comparison.amplitude_effect_size = 0;
        fprintf('    Amplitude: Insufficient data for comparison\n');
    end
    
    %% === ROI Activity Fraction Comparison ===
    % Compare fraction of active ROIs per file (file-level comparison)
    wt_activity_fractions = [];
    mut_activity_fractions = [];
    
    for i = 1:length(wt_data.file_summaries)
        wt_activity_fractions(i) = wt_data.file_summaries(i).active_rois / wt_data.file_summaries(i).num_rois;
    end
    
    for i = 1:length(mut_data.file_summaries)
        mut_activity_fractions(i) = mut_data.file_summaries(i).active_rois / mut_data.file_summaries(i).num_rois;
    end
    
    if length(wt_activity_fractions) >= 2 && length(mut_activity_fractions) >= 2
        [comparison.activity_pvalue, comparison.activity_h] = ranksum(wt_activity_fractions, mut_activity_fractions);
        comparison.activity_effect_size = calculate_effect_size(wt_activity_fractions, mut_activity_fractions);
        
        fprintf('    Activity: WT=%.1f%%, R213W=%.1f%%, p=%.4f\n', ...
            100*mean(wt_activity_fractions), 100*mean(mut_activity_fractions), comparison.activity_pvalue);
    else
        comparison.activity_pvalue = NaN;
        comparison.activity_h = 0;
        comparison.activity_effect_size = 0;
        fprintf('    Activity: Insufficient files for comparison\n');
    end
    
    %% === Events Per ROI Comparison (All ROIs including inactive) ===
    wt_events_per_roi = wt_data.all_roi_event_counts;
    mut_events_per_roi = mut_data.all_roi_event_counts;
    
    if ~isempty(wt_events_per_roi) && ~isempty(mut_events_per_roi)
        [comparison.events_per_roi_pvalue, comparison.events_per_roi_h] = ranksum(wt_events_per_roi, mut_events_per_roi);
        comparison.events_per_roi_effect_size = calculate_effect_size(wt_events_per_roi, mut_events_per_roi);
        
        fprintf('    Events/ROI: WT=%.2f, R213W=%.2f, p=%.4f\n', ...
            mean(wt_events_per_roi), mean(mut_events_per_roi), comparison.events_per_roi_pvalue);
    else
        comparison.events_per_roi_pvalue = NaN;
        comparison.events_per_roi_h = 0;
        comparison.events_per_roi_effect_size = 0;
    end
    
    %% === Summary Statistics ===
    comparison.summary = struct();
    comparison.summary.wt_n_files = wt_data.num_files;
    comparison.summary.mut_n_files = mut_data.num_files;
    comparison.summary.wt_n_rois = wt_data.total_rois;
    comparison.summary.mut_n_rois = mut_data.total_rois;
    comparison.summary.wt_n_events = wt_data.total_events;
    comparison.summary.mut_n_events = mut_data.total_events;
    
    % Count significant differences
    p_values = [comparison.frequency_pvalue, comparison.amplitude_pvalue, ...
                comparison.activity_pvalue, comparison.events_per_roi_pvalue];
    comparison.summary.significant_tests = sum(p_values < 0.05, 'omitnan');
    comparison.summary.total_tests = sum(~isnan(p_values));
    
    fprintf('    Summary: %d/%d tests significant (p < 0.05)\n', ...
        comparison.summary.significant_tests, comparison.summary.total_tests);
end

function effect_size = calculate_effect_size(group1, group2)
    % Calculate Cohen's d effect size for two groups
    %
    % Effect size interpretation:
    %   Small: 0.2, Medium: 0.5, Large: 0.8
    
    if isempty(group1) || isempty(group2) || length(group1) < 2 || length(group2) < 2
        effect_size = 0;
        return;
    end
    
    % Remove NaN values
    group1 = group1(~isnan(group1));
    group2 = group2(~isnan(group2));
    
    if isempty(group1) || isempty(group2)
        effect_size = 0;
        return;
    end
    
    % Calculate means and standard deviations
    mean1 = mean(group1);
    mean2 = mean(group2);
    std1 = std(group1);
    std2 = std(group2);
    
    % Pooled standard deviation
    n1 = length(group1);
    n2 = length(group2);
    pooled_std = sqrt(((n1-1)*std1^2 + (n2-1)*std2^2) / (n1 + n2 - 2));
    
    % Cohen's d
    if pooled_std > 0
        effect_size = abs(mean1 - mean2) / pooled_std;
    else
        effect_size = 0;
    end
end

function save_condition_plots(plot_handles, folder_path)
    % Save condition comparison plots to files
    
    try
        [~, folder_name] = fileparts(folder_path);
        output_dir = fullfile('condition_plots', folder_name);
        
        if ~exist(output_dir, 'dir')
            mkdir(output_dir);
        end
        
        plot_names = fieldnames(plot_handles);
        
        for i = 1:length(plot_names)
            plot_name = plot_names{i};
            fig_handle = plot_handles.(plot_name);
            
            if ishandle(fig_handle)
                timestamp = datestr(now, 'yyyymmdd_HHMMSS');
                output_file = fullfile(output_dir, sprintf('%s_%s_%s.png', ...
                    folder_name, plot_name, timestamp));
                
                saveas(fig_handle, output_file, 'png');
                fprintf('    Saved: %s\n', output_file);
            end
        end
        
        fprintf('  All plots saved to: %s\n', output_dir);
        
    catch ME
        warning('Failed to save condition plots: %s', ME.message);
    end
end

function results_table = export_condition_summary(results, output_file)
    % Export condition comparison results to CSV table
    %
    % Creates a summary table with key metrics for both conditions
    
    try
        wt = results.wt_data;
        mut = results.mut_data;
        comp = results.comparison;
        
        % Create summary table
        metrics = {
            'Number of Files';
            'Total ROIs';
            'Total Events';
            'Active ROIs';
            'Active ROI Fraction (%)';
            'Mean Frequency (Hz, active ROIs)';
            'Std Frequency (Hz, active ROIs)';
            'Mean Amplitude (dF/F)';
            'Std Amplitude (dF/F)';
            'Mean Events per ROI'
        };
        
        wt_values = [
            wt.num_files;
            wt.total_rois;
            wt.total_events;
            wt.stats.frequency.active_rois.count;
            100 * wt.stats.frequency.active_rois.fraction;
            wt.stats.frequency.active_rois.mean;
            wt.stats.frequency.active_rois.std;
            wt.stats.amplitude.mean;
            wt.stats.amplitude.std;
            wt.stats.events_per_roi.mean
        ];
        
        mut_values = [
            mut.num_files;
            mut.total_rois;
            mut.total_events;
            mut.stats.frequency.active_rois.count;
            100 * mut.stats.frequency.active_rois.fraction;
            mut.stats.frequency.active_rois.mean;
            mut.stats.frequency.active_rois.std;
            mut.stats.amplitude.mean;
            mut.stats.amplitude.std;
            mut.stats.events_per_roi.mean
        ];
        
        p_values = [
            NaN;  % Files
            NaN;  % ROIs
            NaN;  % Events
            comp.activity_pvalue;  % Active ROIs (using activity comparison)
            comp.activity_pvalue;  % Active fraction
            comp.frequency_pvalue;
            NaN;  % Frequency std
            comp.amplitude_pvalue;
            NaN;  % Amplitude std
            comp.events_per_roi_pvalue
        ];
        
        % Create table
        results_table = table(metrics, wt_values, mut_values, p_values, ...
            'VariableNames', {'Metric', 'WT', 'R213W', 'P_Value'});
        
        % Write to file if specified
        if nargin > 1 && ~isempty(output_file)
            writetable(results_table, output_file);
            fprintf('  Exported summary table to: %s\n', output_file);
        end
        
    catch ME
        warning('Failed to export condition summary: %s', ME.message);
        results_table = [];
    end
end