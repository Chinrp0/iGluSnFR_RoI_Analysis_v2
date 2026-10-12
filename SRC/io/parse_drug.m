function drug = parse_drug(filename, config)
    % PARSE_DRUG - Perfusion condition from a name like "..._5AP+ttx-01_mean"
    %
    % The drug is the "+token" that follows the paradigm, lower-cased:
    %   CP_Snfr4_Doc2b-WT_10ms_Cs1-c1_5AP+egta-01  -> 'egta'
    %   CP_Snfr4_Doc2b-WT_10ms_Cs1-c1_5AP-01       -> config.drugs.control_label
    % Several "+tokens" are joined with '+' (e.g. '+egta+ttx' -> 'egta+ttx').

    [~, stem] = fileparts(char(filename));
    tok = regexp(stem, '\+([A-Za-z0-9]+)', 'tokens');
    if isempty(tok)
        drug = config.drugs.control_label;
    else
        drug = lower(strjoin(cellfun(@(t) t{1}, tok, 'UniformOutput', false), '+'));
    end
end
