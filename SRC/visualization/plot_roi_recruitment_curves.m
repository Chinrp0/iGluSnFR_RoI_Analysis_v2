function fig = plot_roi_recruitment_curves(results, options)
    [lblA, lblB, colA, colB] = spont_group_labels(results);   % group labels/colours
    % PLOT_ROI_RECRUITMENT_CURVES - Shows ROI recruitment at different frequency thresholds
    % More nuanced view of active ROI fraction
    %
    % Args:
    %   results - Output from spont_batch_analysis()
    %   options - Optional struct with plotting preferences
    
    if nargin < 2
        options = struct();
    end
    
    fprintf('Creating ROI recruitment curves...\n');
    
    %% Extract frequency data
    wt_freqs = results.wt_data.all_roi_frequencies;
    mut_freqs = results.mut_data.all_roi_frequencies;
    
    wt_n_total = length(wt_freqs);
    mut_n_total = length(mut_freqs);
    
    fprintf(['  ' lblA ': %d total ROIs\n'], wt_n_total);
    fprintf(['  ' lblB ': %d total ROIs\n'], mut_n_total);
    
    %% Define frequency thresholds
    thresholds = [0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2];
    
    %% Calculate recruitment at each threshold
    wt_recruitment_count = zeros(size(thresholds));
    mut_recruitment_count = zeros(size(thresholds));
    wt_recruitment_fraction = zeros(size(thresholds));
    mut_recruitment_fraction = zeros(size(thresholds));
    
    for i = 1:length(thresholds)
        wt_recruitment_count(i) = sum(wt_freqs >= thresholds(i));
        mut_recruitment_count(i) = sum(mut_freqs >= thresholds(i));
        wt_recruitment_fraction(i) = 100 * wt_recruitment_count(i) / wt_n_total;
        mut_recruitment_fraction(i) = 100 * mut_recruitment_count(i) / mut_n_total;
    end
    
    %% Create figure
    fig = figure('Name', ['ROI Recruitment Curves: ' lblA ' vs ' lblB], ...
        'Position', [100, 100, 1400, 900]);
    
    %% Subplot 1: Recruitment Fraction vs Threshold
    subplot(2, 3, [1, 2]);
    plot(thresholds, wt_recruitment_fraction, 'b-o', 'LineWidth', 2.5, ...
        'MarkerSize', 8, 'MarkerFaceColor', colA, 'DisplayName', lblA);
    hold on;
    plot(thresholds, mut_recruitment_fraction, 'r-s', 'LineWidth', 2.5, ...
        'MarkerSize', 8, 'MarkerFaceColor', colB, 'DisplayName', lblB);
    
    xlabel('Frequency Threshold (Hz)');
    ylabel('Active ROI Fraction (%)');
    title('ROI Recruitment by Frequency Threshold');
    legend('Location', 'best');
    grid on;
    xlim([0, max(thresholds)]);
    ylim([0, max([wt_recruitment_fraction, mut_recruitment_fraction]) * 1.1]);
    
    %% Subplot 2: Absolute Count
    subplot(2, 3, 3);
    plot(thresholds, wt_recruitment_count, 'b-o', 'LineWidth', 2.5, ...
        'MarkerSize', 8, 'MarkerFaceColor', colA, 'DisplayName', lblA);
    hold on;
    plot(thresholds, mut_recruitment_count, 'r-s', 'LineWidth', 2.5, ...
        'MarkerSize', 8, 'MarkerFaceColor', colB, 'DisplayName', lblB);
    
    xlabel('Frequency Threshold (Hz)');
    ylabel('Active ROI Count');
    title('ROI Count by Frequency Threshold');
    legend('Location', 'best');
    grid on;
    xlim([0, max(thresholds)]);
    
    %% Subplot 3: Difference in Recruitment
    subplot(2, 3, 4);
    recruitment_diff = wt_recruitment_fraction - mut_recruitment_fraction;
    
    bar(thresholds, recruitment_diff, 'FaceColor', [0.5, 0.5, 0.5], 'EdgeColor', 'black');
    hold on;
    yline(0, 'k-', 'LineWidth', 1.5);
    
    xlabel('Frequency Threshold (Hz)');
    ylabel(['Difference in Recruitment (% ' lblA ' - % ' lblB ')']);
    title(['Recruitment Difference (Positive = ' lblA ' has more)']);
    grid on;
    xlim([min(thresholds) - 0.005, max(thresholds) + 0.005]);
    
    %% Subplot 4: Log scale view
    subplot(2, 3, 5);
    semilogx(thresholds(2:end), wt_recruitment_fraction(2:end), 'b-o', 'LineWidth', 2.5, ...
        'MarkerSize', 8, 'MarkerFaceColor', colA, 'DisplayName', lblA);
    hold on;
    semilogx(thresholds(2:end), mut_recruitment_fraction(2:end), 'r-s', 'LineWidth', 2.5, ...
        'MarkerSize', 8, 'MarkerFaceColor', colB, 'DisplayName', lblB);
    
    xlabel('Frequency Threshold (Hz, log scale)');
    ylabel('Active ROI Fraction (%)');
    title('ROI Recruitment (Log Scale)');
    legend('Location', 'best');
    grid on;
    
    %% Subplot 5: Key threshold comparison
    subplot(2, 3, 6); axis off;
    
    % Define key thresholds
    key_thresholds = [0.01, 0.02, 0.05, 0.1];
    
    summary_text = {
        'RECRUITMENT AT KEY THRESHOLDS';
        '';
        ['Threshold | ' lblA ' % | ' lblB ' % | Diff'];
        '----------|------|---------|-----';
    };
    
    for i = 1:length(key_thresholds)
        thresh = key_thresholds(i);
        thresh_idx = find(thresholds == thresh);
        if ~isempty(thresh_idx)
            wt_val = wt_recruitment_fraction(thresh_idx);
            mut_val = mut_recruitment_fraction(thresh_idx);
            diff_val = wt_val - mut_val;
            
            summary_text{end+1} = sprintf('%.3f Hz  | %5.1f | %7.1f | %+5.1f', ...
                thresh, wt_val, mut_val, diff_val);
        end
    end
    
    summary_text{end+1} = '';
    summary_text{end+1} = 'Total ROIs:';
    summary_text{end+1} = sprintf(['  ' lblA ': %d ROIs'], wt_n_total);
    summary_text{end+1} = sprintf(['  ' lblB ': %d ROIs'], mut_n_total);
    summary_text{end+1} = '';
    summary_text{end+1} = 'Interpretation:';
    summary_text{end+1} = ['  Positive diff = ' lblA ' has'];
    summary_text{end+1} = '    more active ROIs';
    summary_text{end+1} = ['  Negative diff = ' lblB ' has'];
    summary_text{end+1} = '    more active ROIs';
    
    text(0.05, 0.95, summary_text, 'FontSize', 9, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', 'Units', 'normalized', 'FontName', 'FixedWidth');
    
    sgtitle(['ROI Recruitment Curves: ' lblA ' vs ' lblB], 'FontSize', 16, 'FontWeight', 'bold');
    
    fprintf('  Created ROI recruitment curves\n');
end
