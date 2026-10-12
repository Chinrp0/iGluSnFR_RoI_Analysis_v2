function isi_ms = parse_isi_ms(filename, config)
    % PARSE_ISI_MS - Paired-pulse inter-stimulus interval (ms) from the filename
    %
    % PPF recordings encode the ISI in the name (e.g. "..._PPF-30ms-01_mean"),
    % so each 2AP file carries its own interval. The pattern lives in
    % config.twoAP.isi_pattern (single source of truth), with one capture group
    % for the numeric interval. Returns the ISI in milliseconds, or errors if the
    % name has no interval - a 2AP file with no parseable ISI cannot be analysed.
    %
    % Inputs:
    %   filename - file name or full path
    %   config   - tracenorm_config() struct (uses config.twoAP.isi_pattern)
    %
    % Output:
    %   isi_ms - inter-stimulus interval in ms (numeric scalar)

    if nargin < 2 || ~isfield(config, 'twoAP') || ~isfield(config.twoAP, 'isi_pattern')
        error('parse_isi_ms:noPattern', ...
            'config.twoAP.isi_pattern is missing; define it in tracenorm_config.');
    end

    [~, stem] = fileparts(char(filename));
    tok = regexpi(stem, config.twoAP.isi_pattern, 'tokens', 'once');
    if isempty(tok)
        error('parse_isi_ms:noMatch', ...
            'File "%s" has no ISI matching pattern "%s".', stem, config.twoAP.isi_pattern);
    end
    isi_ms = str2double(tok{1});
    if ~isfinite(isi_ms) || isi_ms <= 0
        error('parse_isi_ms:badValue', 'Parsed a non-positive ISI (%g ms) from "%s".', isi_ms, stem);
    end
end
