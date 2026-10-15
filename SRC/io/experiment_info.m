function info = experiment_info(folder, files, config)
    % EXPERIMENT_INFO - sensor and experiment date for a data folder
    %
    % Both come from the folder names, e.g.
    %   E:\...\2026-10-06_CP_Ms-hipp_Snfr4f_doc2b_10ms\PPF\PPF_2026-10-11--1\E_RoI_Raw-Mean
    %   -> date '2026-10-06', sensor 'Snfr4f'
    %   date   - nearest folder (data folder upward) whose name STARTS with
    %            yyyy-mm-dd (processing folders like PPF_2026-10-11--1 do not)
    %   sensor - nearest folder name containing a sensor token (SnFR3, Snfr4f,
    %            GluSnFR, gpiSnFR3, iglu3fast, iglusnfr4f, ...); falls back to
    %            the first data file's name
    % config.visualization.sensor / .experiment_date override when non-empty.
    % Either field is '' when it cannot be found.

    if nargin < 2, files = []; end
    if nargin < 3, config = struct(); end
    info = struct('sensor', '', 'date', '');

    parts = strsplit(strrep(char(folder), '/', '\'), '\');
    parts = fliplr(parts(~cellfun(@isempty, parts)));   % nearest first
    sensor_pat = '(gpi)?(i?glu)?snfr(\d[a-z]*)?|iglu(\d[a-z]*)?';

    for k = 1:numel(parts)
        if isempty(info.date)
            d = regexp(parts{k}, '^\d{4}-\d{2}-\d{2}', 'match', 'once');
            if ~isempty(d), info.date = d; end
        end
        if isempty(info.sensor)
            info.sensor = regexpi(parts{k}, sensor_pat, 'match', 'once');
        end
    end
    if isempty(info.sensor) && ~isempty(files)
        info.sensor = regexpi(files(1).name, sensor_pat, 'match', 'once');
    end

    if isfield(config, 'visualization')
        v = config.visualization;
        if isfield(v, 'sensor') && ~isempty(v.sensor), info.sensor = v.sensor; end
        if isfield(v, 'experiment_date') && ~isempty(v.experiment_date)
            info.date = v.experiment_date;
        end
    end
end
