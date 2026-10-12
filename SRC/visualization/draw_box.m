function draw_box(ax, xc, vals, color, width)
    % DRAW_BOX - Lightweight box-and-whisker for one group, drawn natively.
    %
    % A dependency-free, fast replacement for Statistics Toolbox boxplot() (which
    % is slow on its first call per session and pulls in label machinery we don't
    % need). Draws, at x-centre xc, the IQR box (25-75%), the median line and
    % Tukey whiskers (to the most extreme point within 1.5*IQR). Outliers are not
    % drawn as separate marks - the caller overlays the raw points itself.
    %
    % Inputs:
    %   ax    - target axes
    %   xc    - x-centre of the box
    %   vals  - data vector (NaNs ignored)
    %   color - RGB edge/fill colour
    %   width - box width in x units (default 0.5)

    if nargin < 5 || isempty(width), width = 0.5; end
    v = vals(isfinite(vals));
    if isempty(v), return; end

    hw = width / 2;
    if isscalar(v)
        % A single value (e.g. one coverslip): draw just a median tick.
        line(ax, [xc-hw, xc+hw], [v v], 'Color', color, 'LineWidth', 2);
        return;
    end

    q  = quantile(v, [0.25 0.50 0.75]);
    q1 = q(1); med = q(2); q3 = q(3);
    iqr = q3 - q1;
    lo  = min(v(v >= q1 - 1.5*iqr));   % lower whisker cap (within-fence min)
    hi  = max(v(v <= q3 + 1.5*iqr));   % upper whisker cap (within-fence max)

    % IQR box
    patch(ax, [xc-hw xc+hw xc+hw xc-hw], [q1 q1 q3 q3], color, ...
        'FaceAlpha', 0.12, 'EdgeColor', color, 'LineWidth', 1.5);
    % Median
    line(ax, [xc-hw xc+hw], [med med], 'Color', color, 'LineWidth', 2);
    % Whiskers + caps
    line(ax, [xc xc], [q3 hi], 'Color', color, 'LineWidth', 1);
    line(ax, [xc xc], [lo q1], 'Color', color, 'LineWidth', 1);
    line(ax, [xc-hw/2 xc+hw/2], [hi hi], 'Color', color, 'LineWidth', 1);
    line(ax, [xc-hw/2 xc+hw/2], [lo lo], 'Color', color, 'LineWidth', 1);
end
