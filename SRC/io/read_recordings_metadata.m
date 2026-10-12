function recs = read_recordings_metadata(filepath)
    % READ_RECORDINGS_METADATA - Parse the Fiji macro's Recordings.csv
    %
    % One row per recording (written by 5_GluSnFR_icyROI_to_mean_batch). Frames
    % are 1-based and equal the row number in the matching _mean.csv. Columns
    % that do not apply to an experiment type hold the literal "NA" (e.g. Train_Hz
    % for a 1AP file); those parse to NaN / empty here, never to an error.
    % Multi-valued columns are ';'-separated (one value per trial or pulse).
    %
    % Output: struct array, one element per row, with fields
    %   source_filename   - .czi name (Source_Filename)
    %   mean_csv          - data file name (Mean_CSV)
    %   n_rois, n_frames  - expected data-file shape
    %   experiment        - e.g. '1AP', 'PPF', 'HFS'
    %   frame_interval_ms - ms per frame
    %   interval_source, stim_source
    %   n_trials
    %   ttl_stim_frames   - [1 x n] raw TTL stim frames
    %   stim_frames       - [1 x n] refined stim frames (the stim frames used)
    %   refine_shift_frames
    %   ppf_interval_ms, ppf_mode, train_pulses, train_hz
    %   flags             - cellstr of Step1_Flags tokens ({} when NA)
    %   trial_flags       - {1 x n_trials} cellstr of the flags for each trial
    %                       (tokens tagged _T<k>); untagged flags stay in .flags only
    %   row               - row number in Recordings.csv (for messages)

    if ~exist(filepath, 'file')
        error('read_recordings_metadata:missing', 'Metadata file not found: %s', filepath);
    end

    opts = detectImportOptions(filepath, 'Delimiter', ',', 'FileType', 'text');
    opts = setvartype(opts, opts.VariableNames, 'char');   % parse everything ourselves
    opts.VariableNamingRule = 'preserve';
    T = readtable(filepath, opts);

    required = {'Source_Filename', 'Mean_CSV', 'N_ROIs', 'N_Frames', 'Experiment', ...
        'Frame_Interval_ms', 'N_Trials', 'Stim_Frames_Refined'};
    names = T.Properties.VariableNames;
    missing = required(~ismember(lower(required), lower(names)));
    if ~isempty(missing)
        error('read_recordings_metadata:columns', ...
            '%s is missing required column(s): %s', filepath, strjoin(missing, ', '));
    end
    % Every column as an [n x 1] cellstr; absent optional columns read as all NA.
    C = struct();
    allcols = {'Source_Filename', 'Mean_CSV', 'N_ROIs', 'N_Frames', 'Experiment', ...
        'Frame_Interval_ms', 'Interval_Source', 'Stim_Source', 'N_Trials', ...
        'TTL_Stim_Frames', 'Stim_Frames_Refined', 'Refine_Shift_Frames', ...
        'PPF_Interval_ms', 'PPF_Mode', 'Train_Pulses', 'Train_Hz', 'Step1_Flags'};
    for c = 1:numel(allcols)
        C.(allcols{c}) = column_text(T, names, allcols{c});
    end

    n = height(T);
    recs = repmat(empty_rec(), n, 1);
    for i = 1:n
        r = empty_rec();
        r.row               = i + 1;   % +1 for the header line
        r.source_filename   = strtrim(C.Source_Filename{i});
        r.mean_csv          = na_text(C.Mean_CSV{i});
        r.n_rois            = num1(C.N_ROIs{i});
        r.n_frames          = num1(C.N_Frames{i});
        r.experiment        = na_text(C.Experiment{i});
        r.frame_interval_ms = num1(C.Frame_Interval_ms{i});
        r.n_trials          = num1(C.N_Trials{i});
        r.stim_frames       = numlist(C.Stim_Frames_Refined{i});

        r.interval_source     = na_text(C.Interval_Source{i});
        r.stim_source         = na_text(C.Stim_Source{i});
        r.ttl_stim_frames     = numlist(C.TTL_Stim_Frames{i});
        r.refine_shift_frames = numlist(C.Refine_Shift_Frames{i});
        r.ppf_interval_ms     = num1(C.PPF_Interval_ms{i});
        r.ppf_mode            = na_text(C.PPF_Mode{i});
        r.train_pulses        = num1(C.Train_Pulses{i});
        r.train_hz            = num1(C.Train_Hz{i});

        flagtxt = na_text(C.Step1_Flags{i});
        if isempty(flagtxt)
            r.flags = {};
        else
            r.flags = strtrim(strsplit(flagtxt, ';'));
            r.flags = r.flags(~cellfun(@isempty, r.flags));
        end
        nT = max(r.n_trials, 0);
        if ~isfinite(nT), nT = 0; end
        r.trial_flags = cell(1, nT);
        for k = 1:nT
            tagged = ~cellfun(@isempty, regexp(r.flags, sprintf('_T%d(P\\d+)?(\\(|$)', k), 'once'));
            r.trial_flags{k} = r.flags(tagged);
        end

        if isempty(r.mean_csv)   % fall back to the macro's naming convention
            [~, stem] = fileparts(r.source_filename);
            r.mean_csv = [stem '_mean.csv'];
        end
        recs(i) = r;
    end
end

function c = column_text(T, names, name)
    idx = find(strcmpi(names, name), 1);
    if isempty(idx)
        c = repmat({'NA'}, height(T), 1);
    else
        c = cellstr(T.(names{idx}));
    end
end

function r = empty_rec()
    r = struct('row', NaN, 'source_filename', '', 'mean_csv', '', 'n_rois', NaN, ...
        'n_frames', NaN, 'experiment', '', 'frame_interval_ms', NaN, ...
        'interval_source', '', 'stim_source', '', 'n_trials', NaN, ...
        'ttl_stim_frames', [], 'stim_frames', [], 'refine_shift_frames', [], ...
        'ppf_interval_ms', NaN, 'ppf_mode', '', 'train_pulses', NaN, 'train_hz', NaN, ...
        'flags', {{}}, 'trial_flags', {{}});
end

function s = na_text(s)
    % '' for NA / blank, otherwise the trimmed text.
    if isempty(s), s = ''; return; end
    s = strtrim(char(s));
    if strcmpi(s, 'NA') || strcmpi(s, 'NaN'), s = ''; end
end

function v = num1(s)
    % Scalar number, NaN for NA / blank.
    s = na_text(s);
    if isempty(s), v = NaN; else, v = str2double(s); end
end

function v = numlist(s)
    % ';'-separated numbers -> row vector, [] for NA / blank.
    s = na_text(s);
    if isempty(s), v = []; return; end
    v = str2double(strsplit(s, ';'));
    v = v(:)';
end
