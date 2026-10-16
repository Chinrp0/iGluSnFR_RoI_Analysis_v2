function summary = export_roi_sets(analysis_dir)
%% EXPORT_ROI_SETS - colour-coded Fiji ROI sets of the active / failed ROIs
% Takes an F_Analysis_* folder written by run_pipeline and, for every
% recording in its per-ROI table, copies that recording's Fiji ROI zip from
% <run-root>\D_RoI with each ROI's outline coloured by its analysis outcome:
%
%   <F_Analysis>\D_RoI_All\     every ROI, coloured
%   <F_Analysis>\D_RoI_active\  active ROIs only
%   <F_Analysis>\D_RoI_Fail\    failed ROIs only (red)
%   <F_Analysis>\roi_export_summary.csv   one row per ROI: class, colour, reason
%
% A run with several frame-interval sub-folders (e.g. ...\5ms\, ...\10ms\)
% gets these next to each sub-folder's table.
%
% Active rules (the same ROIs the pipeline's per-ROI statistics use):
%   1AP  n_success >= 1 (responded on at least one trial)            green
%   2AP  gate_either (QC passed and pulse 1 or pulse 2 responded):
%          both pulses green, pulse 1 only cyan, pulse 2 only yellow
%   HFS  qc_pass and not noise_excluded and is_responder             green
% Every other ROI is a fail (red).
%
% ROI k of a recording is column Mean(roi-k) of its _mean.csv (step 5
% measures the zip in Manager order) and entry roi-k.roi in its zip (step 4
% names them in that order), so the .roi entries are copied unchanged apart
% from the stroke colour and keep their roi-NNN names.
%
% Usage:
%   export_roi_sets                % folder dialog
%   export_roi_sets('E:\...\1AP_2026-10-11--1\F_Analysis_2026-10-14--1')

    if nargin < 1 || isempty(analysis_dir)
        analysis_dir = uigetdir('', 'Choose the F_Analysis folder');
        if isequal(analysis_dir, 0), summary = table(); return; end
    end
    analysis_dir = char(analysis_dir);
    if ~isfolder(analysis_dir)
        error('export_roi_sets:noFolder', 'Not a folder: %s', analysis_dir);
    end

    % The run-root holds D_RoI; F_Analysis sits directly inside it.
    run_root = fileparts(analysis_dir);
    roi_dir = fullfile(run_root, 'D_RoI');
    if ~isfolder(roi_dir)
        error('export_roi_sets:noRoiDir', ['No D_RoI folder next to %s.\n' ...
            'Pick an F_Analysis folder inside a run-root.'], analysis_dir);
    end
    zips = dir(fullfile(roi_dir, '*.zip'));
    zip_names = {zips.name};

    tables = find_roi_tables(analysis_dir);
    if isempty(tables)
        error('export_roi_sets:noTable', ['No per-ROI table (oneap_per_roi_release.csv, ' ...
            'twoap_per_roi_metrics.csv, hfs_per_roi_metrics.csv) under %s'], analysis_dir);
    end

    summary = table();
    for t = 1:numel(tables)
        tb = tables(t);
        fprintf('\n%s: %s\n', tb.paradigm, tb.path);
        T = readtable(tb.path, 'TextType', 'char', 'Delimiter', ',');
        C = classify_rois(T, tb.paradigm);
        out_dir = fileparts(tb.path);
        S = write_roi_sets(T, C, roi_dir, zip_names, out_dir);
        writetable(S, fullfile(out_dir, 'roi_export_summary.csv'));
        summary = [summary; S]; %#ok<AGROW>
        fprintf('  Wrote D_RoI_All / D_RoI_active / D_RoI_Fail + roi_export_summary.csv to:\n  %s\n', out_dir);
    end
end

%% ------------------------------------------------------------------------
function tables = find_roi_tables(analysis_dir)
    % Per-ROI tables of every paradigm, in the folder or one level down
    % (frame-interval sub-folders).
    names = {'oneap_per_roi_release.csv', '1AP'; ...
             'twoap_per_roi_metrics.csv', '2AP'; ...
             'hfs_per_roi_metrics.csv',   'HFS'};
    tables = struct('path', {}, 'paradigm', {});
    for k = 1:size(names, 1)
        hits = [dir(fullfile(analysis_dir, names{k, 1})); ...
                dir(fullfile(analysis_dir, '*', names{k, 1}))];
        for h = 1:numel(hits)
            tables(end+1) = struct('path', fullfile(hits(h).folder, hits(h).name), ...
                'paradigm', names{k, 2}); %#ok<AGROW>
        end
    end
end

%% ------------------------------------------------------------------------
function C = classify_rois(T, paradigm)
    % Per row of the per-ROI table: class label, ARGB stroke colour, reason.
    n = height(T);
    cls = repmat({'fail'}, n, 1);
    reason = repmat({''}, n, 1);
    switch paradigm
        case '1AP'
            act = T.n_success >= 1;
            cls(act) = {'active'};
            reason(~act & T.n_valid_trials == 0) = {'no valid trials'};
            reason(~act & T.n_valid_trials > 0)  = {'no successful trial'};
        case '2AP'
            ge = as_logical(T.gate_either);
            r1 = as_logical(T.is_responder1);
            r2 = as_logical(T.is_responder2_incr);
            cls(ge & r1 & r2)  = {'both'};
            cls(ge & r1 & ~r2) = {'A1_only'};
            cls(ge & ~r1 & r2) = {'A2_only'};
            qc = isfinite(T.a2_decay_sub);   % NaN only for ROIs that failed noise QC
            reason(~ge & ~qc) = {'failed noise QC'};
            reason(~ge & qc)  = {'no response to either pulse'};
        case 'HFS'
            qc = as_logical(T.qc_pass);
            nx = as_logical(T.noise_excluded);
            rs = as_logical(T.is_responder);
            act = qc & ~nx & rs;
            cls(act) = {'active'};
            reason(~qc) = {'not analysed'};
            reason(qc & nx) = {'noise excluded'};
            reason(qc & ~nx & ~rs) = {'non-responder'};
    end
    C = struct('class', {cls}, 'reason', {reason});
end

function v = as_logical(x)
    % readtable gives logical, numeric or 'true'/'false' text for a logical column.
    if iscell(x), v = strcmpi(x, 'true') | strcmp(x, '1');
    else,         v = logical(x); v(isnan(double(x))) = false;
    end
    v = v(:);
end

function argb = class_colour(cls)
    switch cls
        case {'active', 'both'}, argb = uint32(hex2dec('FF00FF00'));  % green
        case 'A1_only',          argb = uint32(hex2dec('FF00FFFF'));  % cyan
        case 'A2_only',          argb = uint32(hex2dec('FFFFFF00'));  % yellow
        otherwise,               argb = uint32(hex2dec('FFFF0000'));  % red (fail)
    end
end

function name = colour_name(cls)
    switch cls
        case {'active', 'both'}, name = 'green';
        case 'A1_only',          name = 'cyan';
        case 'A2_only',          name = 'yellow';
        otherwise,               name = 'red';
    end
end

%% ------------------------------------------------------------------------
function S = write_roi_sets(T, C, roi_dir, zip_names, out_dir)
    sets = {'D_RoI_All', 'D_RoI_active', 'D_RoI_Fail'};
    for s = 1:numel(sets)
        d = fullfile(out_dir, sets{s});
        if ~isfolder(d), mkdir(d); end
    end

    files = unique(T.filename, 'stable');
    rows = cell(numel(files), 1);
    for f = 1:numel(files)
        sel = find(strcmp(T.filename, files{f}));
        zip_name = match_roi_zip(files{f}, zip_names);
        if isempty(zip_name)
            fprintf('  WARNING: no ROI zip in D_RoI for %s - skipped\n', files{f});
            continue
        end

        tmp = tempname; mkdir(tmp);
        cleanup = onCleanup(@() rmdir(tmp, 's'));
        entries = unzip(fullfile(roi_dir, zip_name), tmp);
        [~, stems] = cellfun(@fileparts, entries, 'UniformOutput', false);
        if numel(entries) ~= numel(sel)
            fprintf(['  WARNING: %s has %d ROIs but the table has %d for %s - ' ...
                'ROI numbers may not match, skipped\n'], zip_name, numel(entries), numel(sel), files{f});
            continue
        end

        roi = T.roi(sel);
        cls = C.class(sel);
        names = arrayfun(@(k) sprintf('roi-%03d', k), roi, 'UniformOutput', false);
        [found, where] = ismember(names, stems);
        if ~all(found)
            fprintf('  WARNING: %s is missing %s - skipped\n', zip_name, strjoin(names(~found), ', '));
            continue
        end

        % Colour every entry in place, then zip each set in ROI-number order
        % so the Manager lists them as roi-001, roi-002, ...
        [~, order] = sort(roi);
        for j = 1:numel(sel)
            set_stroke_colour(entries{where(j)}, class_colour(cls{j}));
        end
        rel = cellfun(@(p) p(numel(tmp)+2:end), entries(where(order)), 'UniformOutput', false);
        is_fail = strcmp(cls(order), 'fail');
        write_zip(fullfile(out_dir, 'D_RoI_All', zip_name),    rel,            tmp);
        write_zip(fullfile(out_dir, 'D_RoI_active', zip_name), rel(~is_fail),  tmp);
        write_zip(fullfile(out_dir, 'D_RoI_Fail', zip_name),   rel(is_fail),   tmp);
        clear cleanup

        n = numel(sel);
        rows{f} = table(repmat(files(f), n, 1), repmat({zip_name}, n, 1), roi, names, cls, ...
            cellfun(@colour_name, cls, 'UniformOutput', false), C.reason(sel), ...
            'VariableNames', {'filename', 'roi_zip', 'roi', 'roi_name', 'class', 'colour', 'fail_reason'});
        fprintf('  %s: %d ROIs, %d active, %d fail\n', zip_name, n, sum(~is_fail), sum(is_fail));
    end
    S = vertcat(rows{:});
    if isempty(S)
        S = table({}, {}, [], {}, {}, {}, {}, 'VariableNames', ...
            {'filename', 'roi_zip', 'roi', 'roi_name', 'class', 'colour', 'fail_reason'});
    end
end

function zip_name = match_roi_zip(csv_name, zip_names)
    % Same rule as macro 5: the ROI core (zip name minus _RoI.zip) equals the
    % recording stem or is the stem plus a processing suffix (_SNRmax, _dF).
    stem = regexprep(csv_name, '_mean\.csv$', '');
    stem = regexprep(stem, '\.(csv|xlsx|xls)$', '');
    cores = regexprep(zip_names, '(_RoI)?\.zip$', '');
    hit = strcmp(cores, stem) | startsWith(cores, [stem '_']);
    zip_name = '';
    if any(hit)
        idx = find(hit);
        zip_name = zip_names{idx(1)};
        if numel(idx) > 1
            fprintf('  WARNING: several ROI zips match %s - using %s\n', csv_name, zip_name);
        end
    end
end

function set_stroke_colour(roi_file, argb)
    % ImageJ .roi header: "Iout" magic, stroke colour = big-endian ARGB int at
    % byte offset 40 (RoiDecoder STROKE_COLOR).
    fid = fopen(roi_file, 'r+', 'ieee-be');
    if fid < 0, error('export_roi_sets:open', 'Cannot open %s', roi_file); end
    c = onCleanup(@() fclose(fid));
    magic = fread(fid, 4, '*char')';
    if ~strcmp(magic, 'Iout')
        error('export_roi_sets:notRoi', 'Not an ImageJ ROI file: %s', roi_file);
    end
    fseek(fid, 40, 'bof');
    fwrite(fid, argb, 'uint32');
end

function write_zip(zip_path, rel_names, root)
    if exist(zip_path, 'file'), delete(zip_path); end
    if isempty(rel_names), return; end      % no empty zips (Fiji cannot open them)
    zip(zip_path, rel_names, root);
end
