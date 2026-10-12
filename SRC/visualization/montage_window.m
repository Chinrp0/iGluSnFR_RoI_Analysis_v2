function [f0, f1] = montage_window(tax, stim_frame, def_pre, def_post, nFrames)
    % MONTAGE_WINDOW - absolute frame window [f0 f1] for a montage panel.
    %
    % Uses config.visualization.x_axis.montage_range_frames when set (an
    % absolute frame interval, clamped to the trace); otherwise falls back to
    % the default stim-centred window [stim_frame-def_pre, stim_frame+def_post].
    % Lets the user pin the montage x-range (e.g. [200 500]) from the config.

    if ~isempty(tax.montage_range_frames)
        r  = tax.montage_range_frames;
        f0 = max(round(min(r)), 1);
        f1 = min(round(max(r)), nFrames);
    else
        f0 = max(stim_frame - def_pre, 1);
        f1 = min(stim_frame + def_post, nFrames);
    end
    if f1 <= f0, f1 = min(f0 + 1, nFrames); end   % guard degenerate window
end
