function [base_mean, n_used] = prestim_base_mean(pre_raw, outlier_col)
    % PRESTIM_BASE_MEAN - Pre-stim dF/F baseline level excluding outlier frames
    %
    % Mean of the pre-stim dF/F frames that prestim_baseline did NOT flag as
    % upward outliers, so the amplitude / threshold zero point uses the same
    % frames as F0. Falls back to the mean of all finite frames when no mask is
    % available or fewer than 5 clean frames remain.
    %
    % Inputs:
    %   pre_raw     - [prestimFrames x 1] dF/F trace over prestim_idx (may hold NaN)
    %   outlier_col - [prestimFrames x 1] logical outlier mask for this ROI, or []
    %
    % Outputs:
    %   base_mean - scalar baseline level (dF/F units)
    %   n_used    - number of frames averaged

    finite = isfinite(pre_raw);
    keep   = finite;
    if ~isempty(outlier_col) && numel(outlier_col) == numel(pre_raw)
        keep = finite & ~outlier_col(:);
    end
    if sum(keep) < 5
        keep = finite;
    end
    n_used    = sum(keep);
    base_mean = mean(pre_raw(keep));
end
