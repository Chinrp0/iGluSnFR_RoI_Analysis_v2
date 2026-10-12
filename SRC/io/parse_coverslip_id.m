function csid = parse_coverslip_id(filename)
    % PARSE_COVERSLIP_ID - Coverslip identity from a name like "...Cs1-c1..."
    %
    % The Cs# token is the experimental unit; the c# region token is the
    % pseudoreplicate and is discarded, so multiple regions on one coverslip
    % collapse to the same unit before any statistics. Falls back to the whole
    % stem when no Cs# is present, so an unparseable file becomes its own unit
    % rather than silently merging with others.
    %
    % Note: coverslip numbering (Cs1, Cs2, ...) typically restarts per genotype
    % folder, so this ID is only unique WITHIN a group. Callers aggregate per
    % group, so that is the intended scope.

    [~, stem] = fileparts(char(filename));
    tok = regexpi(stem, 'Cs[\s_\-]*0*(\d+)', 'tokens', 'once');
    if isempty(tok)
        csid = stem;                     % fallback: treat file as its own unit
    else
        csid = sprintf('Cs%s', tok{1});  % normalized, e.g. 'Cs1'
    end
end
