function share_ylim(axs, config)
    % SHARE_YLIM - give every axes in AXS the same y-range (the union of their
    % current limits, i.e. scaled to the largest trace plotted)
    %
    % config.visualization.link_trace_axes picks how:
    %   true  - linkaxes(axs, 'y'): zooming / panning y in one panel moves all
    %           of them (they stay on one scale while you explore)
    %   false - fixed shared limits only; each panel zooms / pans on its own
    %
    % apply_ymin_zero still sets the floor afterwards; the shared top is kept.
    %
    %   axs    - axes handles (non-axes / deleted handles are ignored)
    %   config - pipeline config (for visualization.link_trace_axes)

    axs = axs(isgraphics(axs, 'axes'));
    if numel(axs) < 2, return; end
    yl = vertcat(axs.YLim);
    yl = [min(yl(:, 1)), max(yl(:, 2))];
    link = true;
    if isfield(config, 'visualization') && isfield(config.visualization, 'link_trace_axes')
        link = config.visualization.link_trace_axes;
    end
    set(axs, 'YLim', yl);
    if link, linkaxes(axs, 'y'); end
end
