function hb = grouped_hist(ax, data, edges, colors, names)
    % GROUPED_HIST - side-by-side relative-frequency histogram of several groups
    %
    % Each bin holds one bar per group next to each other (no overlap), and each
    % group is normalised to its own number of finite values, so groups of
    % different size compare directly. Values outside EDGES are not drawn but
    % still count in the denominator.
    %
    %   ax     - target axes
    %   data   - cell array, one numeric vector per group
    %   edges  - bin edges, or a scalar bin count (edges then span all groups' data)
    %   colors - nG x 3 colour matrix (or cell of 1x3)
    %   names  - cellstr of legend names
    %
    % Returns the bar handles of the non-empty groups (empty if none).

    hb = gobjects(1, 0);
    nG = numel(data);
    data = cellfun(@(v) v(isfinite(v)), data, 'UniformOutput', false);
    has = ~cellfun(@isempty, data);
    if ~any(has), return; end
    if isscalar(edges)
        [~, edges] = histcounts([data{has}], edges);
    end
    edges = edges(:)';
    ctr = (edges(1:end-1) + edges(2:end)) / 2;
    Y = zeros(numel(ctr), nG);
    for g = find(has)
        Y(:, g) = histcounts(data{g}, edges)' / numel(data{g});
    end
    if iscell(colors), colors = vertcat(colors{:}); end
    gi = find(has);
    hb = bar(ax, ctr, Y(:, gi), 1, 'grouped', 'EdgeColor', 'none');
    for k = 1:numel(hb)
        hb(k).FaceColor   = colors(gi(k), :);
        hb(k).DisplayName = names{gi(k)};
    end
end
