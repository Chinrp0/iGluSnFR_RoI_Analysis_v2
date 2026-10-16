function id = recording_id(filename)
    % RECORDING_ID - Short recording ID: coverslip-region-repeat, e.g. 'Cs1-c6-02'
    %
    % From a name like "..._Cs1-c6_1AP-02_mean.csv". Used in figure labels and
    % table ID columns so a panel can be traced to its row and raw file.
    % Falls back to the file stem when the name has no Cs#-c#, and drops the
    % repeat part when there is no trailing "-##".

    [~, stem] = fileparts(char(filename));
    cs = parse_coverslip_id(filename);
    rg = parse_region_id(filename);
    rep = regexp(stem, '-(\d+)(_mean)?$', 'tokens', 'once');
    if strcmp(cs, stem) || strcmp(rg, stem)
        id = stem;
    elseif isempty(rep)
        id = sprintf('%s-%s', cs, rg);
    else
        id = sprintf('%s-%s-%s', cs, rg, rep{1});
    end
end
