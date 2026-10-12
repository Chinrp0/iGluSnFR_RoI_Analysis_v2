function [lo, hi] = twoap_residual_window(stim1, stim2, window_ms, frame_dur)
    % TWOAP_RESIDUAL_WINDOW - Frames used for the pulse-2 local baseline (residual)
    %
    % The residual is the mean trace over the last window_ms before stim2
    % (config.twoAP.residual_window_ms), converted to frames with the
    % recording's frame interval. At least one frame is used, and the window
    % never reaches back to stim1 itself. Shared by twoap_analyzer (per ROI) and
    % twoap_batch_analysis (averaged traces) so both use the same window.
    %
    % Returns lo:hi (frame indices); hi < lo means no usable frame (ISI of one
    % frame).
    if ~(isscalar(window_ms) && isfinite(window_ms) && window_ms > 0)
        error('twoap_residual_window:window', ...
            'config.twoAP.residual_window_ms must be a positive number of ms (got %g).', window_ms);
    end
    n  = max(1, round(window_ms / frame_dur));
    hi = stim2 - 1;
    lo = max(stim2 - n, stim1 + 1);
end
