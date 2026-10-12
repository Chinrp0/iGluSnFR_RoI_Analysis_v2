function fileResult = twoap_pipeline(fileEntry, options, config)
    % TWOAP_PIPELINE - Process ONE 2AP (paired-pulse / PPF) recording
    %
    %   epoch around stim1 -> pre-stim baseline -> dF/F -> twoap_analyzer
    %
    % Timing comes only from the recording's Recordings.csv row (fileEntry.meta):
    %   stim1 = first refined stim frame (Stim_Frames_Refined)
    %   ISI   = PPF_Interval_ms
    %   stim2 = stim1 + round(ISI / Frame_Interval_ms)
    % so the two stims are exactly the reported ISI apart. If the refined stim2,
    % or the ISI in the filename, disagrees with that, the file is flagged in
    % fileResult.isi_check (printed by the batch driver and written to the
    % metrics table) and the ISI-derived stim2 is used.
    %
    % The trace is cut to [stim1 - epoch_pre_ms, stim1 + epoch_post_ms] (NaN
    % padded where the recording is shorter), so stim1 sits at the same index
    % in every file and coverslip x ISI averages line up.
    %
    % Inputs:
    %   fileEntry - one element of dataset.files (from load_dataset_folder)
    %   options   - struct, optional. Fields: verbose (logical)
    %   config    - tracenorm_config() struct, optional
    %
    % Output:
    %   fileResult - struct with filename, group, drug, ISI, dF/F epoch, ISI
    %                check and the twoap_analyzer result.

    if nargin < 2, options = struct(); end
    if ~isfield(options, 'verbose'), options.verbose = false; end
    if nargin < 3 || isempty(config)
        config = struct();
        evalc('config = tracenorm_config();');   % suppress banner
    end
    if ~isstruct(fileEntry) || ~isfield(fileEntry, 'meta')
        error('twoap_pipeline:input', ...
            'twoap_pipeline takes a dataset.files entry (run load_dataset first).');
    end
    config.dataType = '2AP';
    config.verbose  = options.verbose;

    r   = fileEntry.meta;
    raw = fileEntry.raw;
    fi  = r.frame_interval_ms;
    [numFrames, numROIs] = size(raw);

    if numel(r.stim_frames) ~= 2
        error('twoap_pipeline:stims', ...
            '%s: expected 2 refined stim frames (one pair), got %d. Multi-trial PPF is not supported yet.', ...
            fileEntry.name, numel(r.stim_frames));
    end

    % --- ISI and stim2 from the metadata, cross-checked ---------------------
    isi_ms = r.ppf_interval_ms;
    checks = {};
    name_isi = NaN;
    try
        name_isi = parse_isi_ms(fileEntry.name, config);
    catch
    end
    if ~(isfinite(isi_ms) && isi_ms > 0)
        if ~isfinite(name_isi)
            error('twoap_pipeline:isi', '%s: PPF_Interval_ms is NA and no ISI in the filename.', fileEntry.name);
        end
        isi_ms = name_isi;
        checks{end+1} = sprintf('PPF_Interval_ms NA, ISI %g ms from filename', name_isi);
    elseif isfinite(name_isi) && name_isi ~= isi_ms
        checks{end+1} = sprintf('filename ISI %g ms != PPF_Interval_ms %g ms', name_isi, isi_ms);
    end
    isi_frames = round(isi_ms / fi);
    gap = diff(r.stim_frames);
    if gap ~= isi_frames
        checks{end+1} = sprintf('refined stims %d frames apart, ISI %g ms / %.4g ms = %d frames; stim2 set from ISI', ...
            gap, isi_ms, fi, isi_frames);
    end
    isi_check = strjoin(checks, '; ');

    % --- Epoch around stim1 -------------------------------------------------
    s1     = r.stim_frames(1);
    pre_f  = round(config.twoAP.epoch_pre_ms  / fi);
    post_f = round(config.twoAP.epoch_post_ms / fi);
    epLen  = pre_f + post_f + 1;
    ep_s1  = pre_f + 1;
    lo = max(s1 - pre_f, 1); hi = min(s1 + post_f, numFrames);
    data = nan(epLen, numROIs, 'single');
    data((lo:hi) - s1 + ep_s1, :) = raw(lo:hi, :);

    config = apply_recording_timing(config, fi, epLen);
    config.twoAP.stim1_frame = ep_s1;
    config.twoAP.stim2_frame = ep_s1 + isi_frames;
    config.twoAP.isi_ms      = isi_ms;
    if s1 + isi_frames >= numFrames
        error('twoap_pipeline:stim2OutOfRange', ...
            '%s: stim2 (frame %d) is past the end of the recording (%d frames).', ...
            fileEntry.name, s1 + isi_frames, numFrames);
    end

    % --- Pre-stim baseline -> dF/F ---
    [baseline, baseline_stats] = prestim_baseline(data, config);
    [dfof_data, dfof_stats]    = dfof_calculator(data, baseline, config);

    % --- Per-ROI 2AP analysis ---
    twoap = twoap_analyzer(dfof_data, baseline_stats, config);

    % --- Assemble file result ---
    fileResult = struct();
    fileResult.filename       = fileEntry.name;
    fileResult.filepath       = fileEntry.path;
    fileResult.group          = fileEntry.genotype;
    fileResult.genotype       = fileEntry.genotype;
    fileResult.drug           = fileEntry.drug;
    fileResult.meta           = r;
    fileResult.isi_ms         = isi_ms;
    fileResult.isi_check      = isi_check;           % '' when stims match the ISI
    fileResult.ppf_mode       = r.ppf_mode;          % merged | separate (reported only)
    fileResult.numFrames      = epLen;               % frames in the epoch
    fileResult.recFrames      = numFrames;           % frames in the recording
    fileResult.epoch_bounds   = [lo hi];
    fileResult.numROIs        = numROIs;
    fileResult.dfof_data      = dfof_data;
    fileResult.baseline_stats = baseline_stats;
    fileResult.dfof_stats     = dfof_stats;
    fileResult.twoap          = twoap;
    fileResult.config         = config;
end
