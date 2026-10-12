function rid = parse_region_id(filename)
    % PARSE_REGION_ID - Region (c#) token from a name like "...Cs1-c1..."
    %
    % A coverslip (Cs#) is imaged at several regions (c1, c2, ...); each file is
    % one region. parse_coverslip_id keeps the Cs# unit; this keeps the c# region
    % label WITHIN that coverslip, so the individual regions can be overlaid on
    % the coverslip grand-average to show within-coverslip variability.
    %
    % The c# is matched only when it follows the Cs# token (so a stray 'c' in the
    % rest of the name is ignored). Falls back to the whole stem when no c# is
    % present, so an unlabelled file becomes its own region rather than merging.

    [~, stem] = fileparts(char(filename));
    tok = regexpi(stem, 'Cs[\s_\-]*0*\d+[\s_\-]*c[\s_\-]*0*(\d+)', 'tokens', 'once');
    if isempty(tok)
        rid = stem;                      % fallback: treat file as its own region
    else
        rid = sprintf('c%s', tok{1});    % normalized, e.g. 'c1'
    end
end
