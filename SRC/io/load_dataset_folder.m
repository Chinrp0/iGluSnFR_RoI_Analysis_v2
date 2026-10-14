function dataset = load_dataset_folder(folder, config)
    % LOAD_DATASET_FOLDER - Load a raw-mean folder + its Recordings.csv into memory
    %
    % Reads the metadata (Recordings.csv), cross-checks it against the data files
    % in BOTH directions, loads every matched _mean.csv once, and checks each
    % file's shape against the metadata. The returned struct is what the
    % pipeline consumes, so a run can reuse it without re-reading the CSVs (see
    % load_dataset.m / run_pipeline.m).
    %
    % Problems never stop the load; they are warned about, the file is skipped,
    % and everything is listed in dataset.validation and the printed summary:
    %   missing_data  - metadata row whose Mean_CSV is not in the folder
    %   orphan_files  - data file with no metadata row
    %   rejected      - matched file that failed to load or whose frames / ROIs
    %                   disagree with N_Frames / N_ROIs, or whose stim frames are
    %                   missing / out of range
    % Stim frames: check_stim_frames compares the refined and TTL frames with the
    % ROIs' response onset and switches a trial to TTL only when the refined
    % frame is impossible; overrides are printed, added to the file's flags and
    % listed in the summary's stim_check column.
    %
    % Inputs:
    %   folder - raw-mean folder (holds the _mean.csv files and Recordings.csv)
    %   config - tracenorm_config() struct (metadata, groups, drugs)
    %
    % Output: dataset struct with
    %   .folder, .metadata_file, .loaded_at
    %   .files(k)   - .name .path .genotype .drug .coverslip .region .raw .meta
    %   .validation - .missing_data .orphan_files .rejected (cellstr / struct)
    %   .summary    - table, one row per metadata row or orphan file

    if nargin < 2 || isempty(config)
        evalc('config = tracenorm_config();');
    end
    if ~isfolder(folder)
        error('load_dataset_folder:noFolder', 'Data folder not found: %s', folder);
    end

    fprintf('=== Loading dataset ===\nFolder: %s\n', folder);

    %% Metadata file (name matched case-insensitively: Recordings.csv / recordings.csv)
    csvs = dir(fullfile(folder, '*.csv'));
    is_meta = strcmpi({csvs.name}, config.metadata.filename);
    if ~any(is_meta)
        error('load_dataset_folder:noMetadata', ...
            ['No %s in %s. Timing (frames, frame interval, stim frames) is read ' ...
             'only from the metadata, so the folder cannot be analysed without it.'], ...
            config.metadata.filename, folder);
    end
    meta_path = fullfile(folder, csvs(find(is_meta, 1)).name);
    recs = read_recordings_metadata(meta_path);
    fprintf('Metadata: %s (%d recordings)\n', meta_path, numel(recs));

    %% Two-way cross-check: metadata rows <-> data files
    data_files = list_data_files(folder);
    data_files = data_files(~strcmpi({data_files.name}, config.metadata.filename));
    data_names = {data_files.name};
    meta_names = {recs.mean_csv};

    [~, dup_idx] = unique(lower(meta_names), 'stable');
    if numel(dup_idx) < numel(meta_names)
        dups = meta_names(setdiff(1:numel(meta_names), dup_idx));
        error('load_dataset_folder:duplicateMetadata', ...
            '%s lists the same data file more than once: %s', ...
            config.metadata.filename, strjoin(unique(dups), ', '));
    end

    in_folder = ismember(lower(meta_names), lower(data_names));
    has_meta  = ismember(lower(data_names), lower(meta_names));
    missing_data = meta_names(~in_folder);
    orphan_files = data_names(~has_meta);
    for k = 1:numel(missing_data)
        warning('load_dataset_folder:missingData', ...
            'In %s but not in the folder (skipped): %s', config.metadata.filename, missing_data{k});
    end
    for k = 1:numel(orphan_files)
        warning('load_dataset_folder:orphanFile', ...
            'Data file has no %s row (skipped): %s', config.metadata.filename, orphan_files{k});
    end

    %% Load + validate each matched file
    loader   = csv_loader_v2();
    files    = struct('name', {}, 'path', {}, 'genotype', {}, 'drug', {}, ...
        'coverslip', {}, 'region', {}, 'dataType', {}, 'raw', {}, 'meta', {});
    rejected = struct('name', {}, 'reason', {});
    status   = repmat({'loaded'}, numel(recs), 1);
    stim_note = repmat({''}, numel(recs), 1);
    status(~in_folder) = {'missing data file'};

    for i = find(in_folder)
        r = recs(i);
        fpath = fullfile(folder, data_names{find(strcmpi(data_names, r.mean_csv), 1)});
        reason = check_metadata_row(r);
        raw = [];
        if isempty(reason)
            try
                raw = loader.loadSingleFile(fpath, struct('dataType', 'single'));
                reason = check_shape(raw, r);
            catch ME
                reason = ME.message;
            end
        end
        if isempty(reason)
            try
                genotype = parse_group(r.mean_csv, config);
            catch ME
                reason = ME.message;
            end
        end
        if ~isempty(reason)
            warning('load_dataset_folder:rejected', '%s skipped: %s', r.mean_csv, reason);
            rejected(end+1) = struct('name', r.mean_csv, 'reason', reason); %#ok<AGROW>
            status{i} = ['rejected: ' reason];
            continue;
        end

        % Refined vs TTL stim frames, decided from the ROIs' own response onset
        % (the refined frame comes from the full field of view).
        if ~strcmp(experiment_to_datatype(r.experiment), 'Spont')
            [r_chk, stim_note{i}] = check_stim_frames(raw, r, config);
            if isempty(check_metadata_row(r_chk))
                r = r_chk;
            else
                stim_note{i} = [stim_note{i} ' (TTL frames invalid, refined kept)'];
            end
            if ~isempty(stim_note{i})
                fprintf('  Stim check %s: %s\n', r.mean_csv, stim_note{i});
            end
        end

        f = struct();
        f.name      = r.mean_csv;
        f.path      = fpath;
        f.genotype  = genotype;
        f.drug      = parse_drug(r.mean_csv, config);
        f.coverslip = parse_coverslip_id(r.mean_csv);
        f.region    = parse_region_id(r.mean_csv);
        f.dataType  = experiment_to_datatype(r.experiment);
        f.raw       = raw;
        f.meta      = r;
        files(end+1) = f; %#ok<AGROW>
    end

    %% Summary table (metadata rows first, then orphan data files)
    nM = numel(recs);
    file_col   = [meta_names(:); orphan_files(:)];
    status_col = [status; repmat({'no metadata row'}, numel(orphan_files), 1)];
    frames_col = [[recs.n_frames]'; nan(numel(orphan_files), 1)];
    rois_col   = [[recs.n_rois]';   nan(numel(orphan_files), 1)];
    fi_col     = [[recs.frame_interval_ms]'; nan(numel(orphan_files), 1)];
    trials_col = [[recs.n_trials]'; nan(numel(orphan_files), 1)];
    exp_col    = [{recs.experiment}'; repmat({''}, numel(orphan_files), 1)];
    nflag_col  = [cellfun(@numel, {recs.flags})'; nan(numel(orphan_files), 1)];
    stim_col   = [stim_note; repmat({''}, numel(orphan_files), 1)];
    summary = table(file_col, exp_col, frames_col, rois_col, fi_col, trials_col, nflag_col, status_col, stim_col, ...
        'VariableNames', {'file', 'experiment', 'n_frames', 'n_rois', 'frame_interval_ms', ...
        'n_trials', 'n_flags', 'status', 'stim_check'});

    dataset = struct();
    dataset.folder        = folder;
    dataset.metadata_file = meta_path;
    dataset.loaded_at     = datetime('now');
    dataset.files         = files;
    % Paradigm(s) present, from the Experiment column (run_pipeline requires one).
    if isempty(files), dataset.dataType = {}; else, dataset.dataType = unique({files.dataType}, 'stable'); end
    dataset.recordings    = recs;
    dataset.validation    = struct('missing_data', {missing_data}, ...
        'orphan_files', {orphan_files}, 'rejected', rejected);
    dataset.summary       = summary;

    fprintf('Loaded %d/%d metadata recordings | %d missing data file | %d orphan file | %d rejected\n', ...
        numel(files), nM, numel(missing_data), numel(orphan_files), numel(rejected));
    if numel(files) < nM || ~isempty(orphan_files)
        disp(summary(~strcmp(summary.status, 'loaded'), :));
    end
    print_group_counts(files);
end

function reason = check_metadata_row(r)
    % Timing fields the pipeline relies on must be present and self-consistent.
    % Stim frames are required for evoked paradigms only (spontaneous has none).
    reason = '';
    dt = experiment_to_datatype(r.experiment);
    if isempty(dt)
        reason = sprintf('Experiment "%s" has no analysis path', r.experiment);
    elseif ~(isfinite(r.n_frames) && r.n_frames >= 1)
        reason = 'N_Frames is NA/invalid';
    elseif ~(isfinite(r.frame_interval_ms) && r.frame_interval_ms > 0)
        reason = 'Frame_Interval_ms is NA/invalid';
    elseif strcmp(dt, 'Spont')
        return;
    elseif isempty(r.stim_frames) || any(~isfinite(r.stim_frames))
        reason = 'Stim_Frames_Refined is NA/invalid';
    elseif any(r.stim_frames < 1 | r.stim_frames > r.n_frames)
        reason = sprintf('Stim_Frames_Refined outside 1..%d', r.n_frames);
    elseif any(diff(r.stim_frames) <= 0)
        reason = 'Stim_Frames_Refined not increasing';
    end
end

function reason = check_shape(raw, r)
    reason = '';
    [nF, nR] = size(raw);
    if nF ~= r.n_frames
        reason = sprintf('%d frames in file, N_Frames = %d', nF, r.n_frames);
    elseif isfinite(r.n_rois) && nR ~= r.n_rois
        reason = sprintf('%d ROI columns in file, N_ROIs = %d', nR, r.n_rois);
    end
end

function print_group_counts(files)
    if isempty(files), return; end
    key = strcat({files.genotype}, {' | '}, {files.drug});
    [u, ~, j] = unique(key, 'stable');
    for k = 1:numel(u)
        fprintf('  %-28s %d file(s)\n', u{k}, sum(j == k));
    end
end
