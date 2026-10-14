function mask = prestim_outlier_mask(baseline_stats, nPre, numROIs)
    % PRESTIM_OUTLIER_MASK - Pre-stim outlier mask from prestim_baseline, or all-false
    %
    % Returns baseline_stats.outlier_mask when present and sized [nPre x numROIs];
    % otherwise an all-false mask (no exclusion, i.e. the plain pre-stim mean).

    mask = false(nPre, numROIs);
    if isstruct(baseline_stats) && isfield(baseline_stats, 'outlier_mask') ...
            && isequal(size(baseline_stats.outlier_mask), [nPre, numROIs])
        mask = logical(baseline_stats.outlier_mask);
    end
end
