function [nr, nc, per_page, n_pages] = montage_grid(nR, mc)
    % MONTAGE_GRID - tile layout for an ROI montage of nR panels
    %
    % mc = config.visualization.montage. Columns come from mc.columns. With
    % mc.rows empty (or absent) every panel goes on one figure (rows = as many
    % as needed). With mc.rows set, each figure is a fixed rows x columns grid
    % and panels beyond that continue on further figures (pages), so the
    % max_rois settings still decide how many ROIs are shown in total.
    %
    % Outputs: nr x nc grid, per_page panels per figure, n_pages figures.
    nc = max(1, mc.columns);
    if isfield(mc, 'rows') && ~isempty(mc.rows)
        nr = max(1, mc.rows);
    else
        nr = max(1, ceil(nR / nc));
    end
    per_page = nr * nc;
    n_pages  = max(1, ceil(nR / per_page));
end
