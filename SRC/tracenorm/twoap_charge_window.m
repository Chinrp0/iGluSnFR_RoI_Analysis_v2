function [lo, hi] = twoap_charge_window(stim1, stim2, nF, post_ms, frame_dur)
    % TWOAP_CHARGE_WINDOW - Frames integrated for the 2AP total charge Q
    %
    % Q (q_pair) is the integral of dF/F above the pre-stim1 baseline over ONE
    % window spanning both pulses: stim1 .. stim2 + post_ms
    % (config.twoAP.auc_post_ms), converted to frames with the recording's
    % frame interval and clipped to the trace. Shared by twoap_analyzer (per
    % ROI), twoap_batch_analysis (averaged traces) and the figures that shade
    % the window, so the shading is exactly the integrated frames.
    %
    % Returns lo:hi (frame indices).
    if ~(isscalar(post_ms) && isfinite(post_ms) && post_ms > 0)
        error('twoap_charge_window:window', ...
            'config.twoAP.auc_post_ms must be a positive number of ms (got %g).', post_ms);
    end
    lo = stim1;
    hi = min(stim2 + max(1, round(post_ms / frame_dur)), nF);
end
