function scroll_time_axis(ax, xlims)
    % SCROLL_TIME_AXIS - configure an axes to scroll horizontally through a trace.
    %
    % The FULL trace is plotted, but only xlims is shown initially; the user can
    % click-drag to pan left/right (scroll through the rest of the trace) and
    % scroll-wheel to zoom. Panning and zooming are constrained to the x
    % dimension so the vertical (dF/F) scale stays put while scrolling. This
    % replaces cramming the whole windowed slice into a tiny panel: the region of
    % interest fills the panel and the rest of the trace is a drag away.
    %
    % Inputs:
    %   ax    - target axes (default gca)
    %   xlims - initial [lo hi] x-limits in the axes' display units
    if nargin < 1 || isempty(ax), ax = gca; end
    if nargin < 2 || isempty(xlims) || numel(xlims) < 2 || xlims(2) <= xlims(1)
        return;
    end
    xlim(ax, xlims);
    ax.Toolbar.Visible = 'on';
    try
        % Configurable per-axes interactions (R2019a+). Constrain to x so drag =
        % horizontal scroll and wheel = horizontal zoom; keep datatips.
        ax.Interactions = [ ...
            panInteraction('Dimensions', 'x'), ...
            zoomInteraction('Dimensions', 'x'), ...
            dataTipInteraction];
    catch
        % Older MATLAB: default pan/zoom from the axes toolbar still work.
    end
end
