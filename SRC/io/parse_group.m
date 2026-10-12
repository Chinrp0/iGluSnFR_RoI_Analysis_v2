function name = parse_group(filename, config)
    % PARSE_GROUP - Map a filename to a configured group label (fail-loud)
    %
    % The group (e.g. WT vs R213W) is derived by matching the filename stem
    % against the case-insensitive patterns in config.groups. This is the SINGLE
    % mechanism every path uses for group identity, so swapping proteins is a
    % one-line edit in tracenorm_config (config.groups) and nothing else.
    %
    % A name that matches no pattern - or more than one - is an error, not a
    % silent drop or a guess: mislabelling a genotype is worse than stopping.
    %
    % Inputs:
    %   filename - file name or full path (extension and folders are ignored)
    %   config   - tracenorm_config() struct (uses config.groups)
    %
    % Output:
    %   name - the matched group's .name

    if nargin < 2 || ~isfield(config, 'groups') || isempty(config.groups)
        error('parse_group:noGroups', ...
            'config.groups is missing/empty; define groups in tracenorm_config.');
    end

    [~, stem] = fileparts(char(filename));

    groups  = config.groups;
    matched = false(1, numel(groups));
    for g = 1:numel(groups)
        if ~isempty(regexpi(stem, groups(g).pattern, 'once'))
            matched(g) = true;
        end
    end

    hits = find(matched);
    if isempty(hits)
        error('parse_group:noMatch', ...
            ['File "%s" matches no configured group pattern (%s). ' ...
             'Add a pattern to config.groups or fix the filename.'], ...
            stem, strjoin({groups.name}, ', '));
    elseif numel(hits) > 1
        error('parse_group:ambiguous', ...
            'File "%s" matches multiple groups (%s); patterns must be mutually exclusive.', ...
            stem, strjoin({groups(hits).name}, ', '));
    end

    name = groups(hits).name;
end
