function n = apply_ymin_zero(figs)
    % APPLY_YMIN_ZERO - set the lower y-limit of every axes to 0 (or its own floor)
    %
    % Used by run_pipeline (config.visualization.y_min_zero) so all figures
    % share a zero-based y-axis without each plot function having to do it.
    % The upper limit is kept (MATLAB's current auto/fixed top).
    %
    % An axes can ask for a different floor with  setappdata(ax, 'ymin', v)
    % - dF/F trace axes do this so the negative baseline noise stays visible
    % (config.visualization.ymin_montage / ymin_mean_trace).
    %
    % Skipped:
    %   - log-scale y axes (0 is not a valid limit there)
    %   - axes holding an image (imagesc rows; y is an index, not a value)
    %   - axes whose current top is at or below the floor
    %
    %   figs - figure handles (default: all open figures)
    %   n    - number of axes changed

    if nargin < 1, figs = findall(groot, 'Type', 'figure'); end
    n = 0;
    for f = figs(:)'
        axs = findall(f, 'Type', 'axes');
        for ax = axs(:)'
            if strcmpi(ax.YScale, 'log'), continue; end
            if ~isempty(findobj(ax, 'Type', 'image')), continue; end
            lo = 0;
            if isappdata(ax, 'ymin'), lo = getappdata(ax, 'ymin'); end
            yl = ax.YLim;
            if yl(2) <= lo, continue; end
            ax.YLim = [lo yl(2)];
            n = n + 1;
        end
    end
end
