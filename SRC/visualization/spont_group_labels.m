function [lblA, lblB, colA, colB] = spont_group_labels(a, b)
    % SPONT_GROUP_LABELS - Names + colours of the two groups in a spont comparison
    %
    % The spontaneous plots compare two groups held in the slots wt_data (A)
    % and mut_data (B). The slot names are historical; the actual groups (e.g.
    % "Doc2b-WT control" vs "Doc2b-R213W control") are set by
    % spont_batch_analysis in each slot's .condition / .color.
    %
    %   spont_group_labels(results)            - results with .wt_data / .mut_data
    %   spont_group_labels(wt_data, mut_data)  - the two condition structs

    if nargin == 1
        b = a.mut_data;
        a = a.wt_data;
    end
    [lblA, colA] = label_of(a, 'WT',    [0.2 0.6 1.0]);
    [lblB, colB] = label_of(b, 'R213W', [1.0 0.4 0.2]);
end

function [lbl, col] = label_of(d, def_lbl, def_col)
    lbl = def_lbl; col = def_col;
    if isfield(d, 'condition') && ~isempty(d.condition), lbl = d.condition; end
    if isfield(d, 'color') && ~isempty(d.color), col = d.color; end
end
