function fileResult = hfs_pipeline(fileEntry, options, config)
    % HFS_PIPELINE - Process ONE high-frequency train (HFS) recording
    %
    %   epoch around stim1 -> pre-stim baseline -> dF/F -> hfs_analyzer
    %
    % Timing comes only from the recording's Recordings.csv row (fileEntry.meta):
    %   stim1  = first refined stim frame (Stim_Frames_Refined)
    %   pulses = Train_Pulses, rate = Train_Hz (apply_train_timing places them)
    % The trace is cut to [stim1 - epoch_pre_ms, train end + post_train_ms]
    % (NaN padded where the recording is shorter), so stim1 sits at the same
    % index in every file and coverslip averages line up.
    %
    % fileResult.train_check lists anything worth a look before interpreting
    % the file: a filename train ("50AP-20Hz") that disagrees with the metadata,
    % and phasic peaks that sit on the stim frames (median latency below
    % config.hfs.min_peak_latency_ms on the file's mean trace), which usually
    % means the refined stim frame is late.
    %
    % Inputs:
    %   fileEntry - one element of dataset.files (from load_dataset_folder)
    %   options   - struct, optional. Fields: verbose (logical)
    %   config    - tracenorm_config() struct, optional
    %
    % Output:
    %   fileResult - struct with filename, genotype, drug, dF/F epoch, checks
    %                and the hfs_analyzer result.

    if nargin < 2, options = struct(); end
    if ~isfield(options, 'verbose'), options.verbose = false; end
    if nargin < 3 || isempty(config)
        evalc('config = tracenorm_config();');
    end
    if ~isstruct(fileEntry) || ~isfield(fileEntry, 'meta')
        error('hfs_pipeline:input', 'hfs_pipeline takes a dataset.files entry (run load_dataset first).');
    end
    config.dataType = 'HFS';
    config.verbose  = options.verbose;

    r   = fileEntry.meta;
    raw = fileEntry.raw;
    fi  = r.frame_interval_ms;
    [numFrames, numROIs] = size(raw);
    nP = r.train_pulses; hz = r.train_hz;
    if ~(isfinite(nP) && isfinite(hz))
        error('hfs_pipeline:train', '%s: Train_Pulses / Train_Hz are NA in %s.', ...
            fileEntry.name, config.metadata.filename);
    end
    checks = {};
    tok = regexpi(fileEntry.name, '(\d+)\s*AP[\s_\-]*(\d+(?:\.\d+)?)\s*Hz', 'tokens', 'once');
    if ~isempty(tok) && (str2double(tok{1}) ~= nP || str2double(tok{2}) ~= hz)
        checks{end+1} = sprintf('filename says %sAP-%sHz, metadata %gAP-%gHz (metadata used)', ...
            tok{1}, tok{2}, nP, hz);
    end

    % --- Epoch around stim1 -------------------------------------------------
    s1     = r.stim_frames(1);
    pre_f  = round(config.hfs.epoch_pre_ms / fi);
    post_f = round((nP * 1000 / hz + config.hfs.post_train_ms) / fi);
    epLen  = pre_f + post_f + 1;
    ep_s1  = pre_f + 1;
    lo = max(s1 - pre_f, 1); hi = min(s1 + post_f, numFrames);
    data = nan(epLen, numROIs, 'single');
    data((lo:hi) - s1 + ep_s1, :) = raw(lo:hi, :);

    config = apply_recording_timing(config, fi, epLen);
    config = apply_train_timing(config, ep_s1, nP, hz);
    if s1 + config.hfs.train_end_frame - ep_s1 > numFrames
        checks{end+1} = sprintf('train runs past the end of the recording (%d frames)', numFrames);
    end

    % --- Pre-stim baseline -> dF/F -> per-ROI metrics -----------------------
    [baseline, baseline_stats] = prestim_baseline(data, config);
    [dfof_data, dfof_stats]    = dfof_calculator(data, baseline, config);
    hfs = hfs_analyzer(dfof_data, baseline_stats, config);

    % File mean trace (QC-passed ROIs) for the latency check.
    mean_trace = mean(dfof_data(:, hfs.analyzed), 2, 'omitnan');
    % Peaks pinned to either edge of the per-pulse window mean the stim frames
    % are off by a few frames: early edge = stim late (peak on the stim frame),
    % late edge = the window is catching the NEXT pulse's peak.
    mt = hfs_train_metrics(mean_trace, [], config, false);
    win_end_ms = min(round(config.hfs.peak_search_ms / fi), round(config.hfs.isi_ms / fi) - 1) * fi;
    tol = config.hfs.min_peak_latency_ms;
    if mt.resolvable && (mt.median_lat_ms < tol || mt.median_lat_ms > win_end_ms - tol)
        checks{end+1} = sprintf(['phasic peaks %.1f ms after the stim frames (median; window 0-%.0f ms): ' ...
            'stim frames probably misaligned, per-pulse values unreliable (Refine_Shift_Frames = %s)'], ...
            mt.median_lat_ms, win_end_ms, num2str(r.refine_shift_frames));
    end

    fileResult = struct();
    fileResult.filename       = fileEntry.name;
    fileResult.filepath       = fileEntry.path;
    fileResult.genotype       = fileEntry.genotype;
    fileResult.drug           = fileEntry.drug;
    fileResult.coverslip      = fileEntry.coverslip;
    fileResult.region         = fileEntry.region;
    fileResult.meta           = r;
    fileResult.train_check    = strjoin(checks, '; ');
    fileResult.median_lat_ms  = mt.median_lat_ms;
    fileResult.numFrames      = epLen;
    fileResult.recFrames      = numFrames;
    fileResult.epoch_bounds   = [lo hi];
    fileResult.numROIs        = numROIs;
    fileResult.dfof_data      = dfof_data;
    fileResult.baseline_stats = baseline_stats;
    fileResult.dfof_stats     = dfof_stats;
    fileResult.hfs            = hfs;
    fileResult.config         = config;
end
