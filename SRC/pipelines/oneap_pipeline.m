function fileResult = oneap_pipeline(fileEntry, options, config)
    % ONEAP_PIPELINE - Process ONE 1AP recording (one or more trials) per ROI
    %
    % Each trial is analysed on its own:
    %   cut epoch at the trial's refined stim frame -> pre-stim baseline ->
    %   dF/F -> oneap_analyzer
    % so every (ROI, trial) pair gets its own success/failure call. The analyzer
    % outputs are then laid side by side as "trial units": column u of every
    % per-ROI array is ROI unit_roi(u) on trial unit_trial(u). Downstream code
    % that treats columns as traces (averages, montages, coverslip means) works
    % on these units unchanged; per-ROI release statistics across trials are in
    % fileResult.roi_summary.
    %
    % Timing comes only from the recording's Recordings.csv row (fileEntry.meta):
    % the frame interval, N_Frames and the refined stim frames.
    %
    % Inputs:
    %   fileEntry - one element of dataset.files (from load_dataset_folder)
    %   options   - struct, optional. Fields: verbose (logical)
    %   config    - tracenorm_config() struct, optional
    %
    % Output:
    %   fileResult - struct with filename, genotype, drug, trial-unit dF/F epochs,
    %                merged analyzer result (.oneap), per-ROI summary and flags.

    if nargin < 2, options = struct(); end
    if ~isfield(options, 'verbose'), options.verbose = false; end
    if nargin < 3 || isempty(config)
        config = struct();
        evalc('config = tracenorm_config();');   % suppress banner
    end
    if ~isstruct(fileEntry) || ~isfield(fileEntry, 'meta')
        error('oneap_pipeline:input', ...
            'oneap_pipeline takes a dataset.files entry (run load_dataset first).');
    end
    config.dataType = '1AP';
    config.verbose  = options.verbose;

    r     = fileEntry.meta;
    raw   = fileEntry.raw;
    stims = r.stim_frames;
    [numFrames, numROIs] = size(raw);
    nT = numel(stims);
    if isfinite(r.n_trials) && r.n_trials ~= nT
        error('oneap_pipeline:trials', ...
            '%s: N_Trials = %d but %d refined stim frames listed.', fileEntry.name, r.n_trials, nT);
    end

    % --- Trial epochs (fixed length, stim at the same index in every epoch) ---
    [pre_f, post_f] = epoch_frames(config, r.frame_interval_ms);
    epLen     = pre_f + post_f + 1;
    ep_stim   = pre_f + 1;
    config    = apply_recording_timing(config, r.frame_interval_ms, epLen);
    config.oneAP.stim_frame = ep_stim;

    dfof_units = nan(epLen, numROIs * nT, 'single');
    trial_res  = cell(1, nT);
    ep_bounds  = zeros(nT, 2);
    for k = 1:nT
        s = stims(k);
        % Never reach into a neighbouring trial: start after the previous trial's
        % post window, stop before the next stim.
        lo = s - pre_f;  hi = s + post_f;
        if k > 1,  lo = max(lo, stims(k-1) + post_f + 1); end
        if k < nT, hi = min(hi, stims(k+1) - 1);          end
        lo = max(lo, 1); hi = min(hi, numFrames);
        ep_bounds(k, :) = [lo hi];

        ep = nan(epLen, numROIs, 'single');
        ep((lo:hi) - s + ep_stim, :) = raw(lo:hi, :);

        [baseline, baseline_stats] = prestim_baseline(ep, config);
        dfof = dfof_calculator(ep, baseline, config);
        trial_res{k} = oneap_analyzer(dfof, baseline_stats, config);
        dfof_units(:, (k-1)*numROIs + (1:numROIs)) = dfof;
    end

    oneap = merge_trials(trial_res, numROIs);
    oneap.unit_roi   = repmat(1:numROIs, 1, nT);
    oneap.unit_trial = repelem(1:nT, numROIs);
    oneap.unit_flags = repelem(cellfun(@(c) strjoin(c, ';'), pad_flags(r, nT), ...
        'UniformOutput', false), numROIs);

    % --- Per-ROI release statistics across trials --------------------------
    % valid trial = the ROI could be scored on it (qc codes 1-3: responder or
    % subthreshold; not all-NaN / too noisy). Failures count as 0 amplitude for
    % efficacy and are left out of potency.
    amp   = reshape(oneap.peak_amplitude, numROIs, nT);
    succ  = reshape(oneap.is_responder,   numROIs, nT);
    valid = reshape(oneap.qc_category >= 1 & oneap.qc_category <= 3, numROIs, nT);
    n_valid   = sum(valid, 2)';
    n_success = sum(succ & valid, 2)';
    amp_s = amp; amp_s(~(succ & valid)) = NaN;
    amp_e = amp; amp_e(valid & ~succ) = 0; amp_e(~valid) = NaN;
    rs = struct();
    rs.n_valid    = n_valid;
    rs.n_success  = n_success;
    rs.p_response = n_success ./ n_valid;            % NaN when no valid trial
    rs.potency    = mean(amp_s, 2, 'omitnan')';      % mean amplitude of successes
    rs.efficacy   = mean(amp_e, 2, 'omitnan')';      % failures counted as 0
    rs.efficacy(n_valid == 0) = NaN;

    % --- Assemble file result ---
    fileResult = struct();
    fileResult.filename       = fileEntry.name;
    fileResult.filepath       = fileEntry.path;
    fileResult.genotype       = fileEntry.genotype;
    fileResult.drug           = fileEntry.drug;
    fileResult.coverslip      = fileEntry.coverslip;
    fileResult.meta           = r;
    fileResult.numFrames      = epLen;              % frames per trial epoch
    fileResult.recFrames      = numFrames;          % frames in the recording
    fileResult.numROIs        = numROIs;            % physical ROIs
    fileResult.numTrials      = nT;
    fileResult.numUnits       = numROIs * nT;       % (ROI, trial) columns
    fileResult.epoch_bounds   = ep_bounds;          % [start stop] recording frames
    fileResult.epoch_stim     = ep_stim;
    fileResult.dfof_data      = dfof_units;
    fileResult.trials         = trial_res;
    fileResult.oneap          = oneap;
    fileResult.roi_summary    = rs;
    fileResult.config         = config;
end

function [pre_f, post_f] = epoch_frames(config, frame_interval_ms)
    pre_f  = round(config.oneAP.epoch_pre_ms  / frame_interval_ms);
    post_f = round(config.oneAP.epoch_post_ms / frame_interval_ms);
end

function c = pad_flags(r, nT)
    c = r.trial_flags;
    if numel(c) < nT, c(end+1:nT) = {{}}; end
    c = c(1:nT);
end

function m = merge_trials(trial_res, numROIs)
    % Concatenate the per-ROI arrays of every trial's analyzer result along the
    % ROI dimension (trial-major), keeping scalar / shared fields from trial 1.
    m = trial_res{1};
    rowFields = {'peak_amplitude', 'peak_frame', 'time_to_peak_ms', 'is_responder', ...
        'prestim_sd', 'upper_threshold', 'lower_threshold', 'auc', 'auc_sync', ...
        'auc_async', 'tau_ms', 'tau_r2', 'tau_A', 'tau_offset', 'tau_peak_frame', ...
        'prestim_base', 'sync_count', 'async_count', 'n_spont', 'total_events', ...
        'event_peak_times_ms', 'release_events', 'crossed_upper', 'noise_excluded', 'qc_category'};
    for f = 1:numel(rowFields)
        name = rowFields{f};
        m.(name) = cell2mat_or_cat(cellfun(@(t) t.(name), trial_res, 'UniformOutput', false));
    end
    % Event bin lengths: from the first trial that ran the detector
    for k = 1:numel(trial_res)
        if isfinite(trial_res{k}.events_ms.async), m.events_ms = trial_res{k}.events_ms; break; end
    end
    m.temporal_hist = cell2mat(cellfun(@(t) t.temporal_hist, trial_res, 'UniformOutput', false));
    m.numROIs = numROIs * numel(trial_res);

    c = m.qc_category;
    m.n_responders = sum(m.is_responder);
    m.responder_fraction = m.n_responders / max(numel(c), 1);
    m.responder_fraction_analyzable = m.n_responders / max(sum(c >= 1 & c <= 3), 1);
    m.n_sync_responders  = sum(c == 1);
    m.n_async_responders = sum(c == 2);
    m.total_sync_events  = sum(m.sync_count, 'omitnan');
    m.total_async_events = sum(m.async_count, 'omitnan');
    m.qc_counts = struct('excluded', sum(c == 0), 'responder_sync', sum(c == 1), ...
        'responder_async', sum(c == 2), 'subthreshold', sum(c == 3), ...
        'noise_excluded', sum(c == 4));
end

function out = cell2mat_or_cat(parts)
    if iscell(parts{1})
        out = [parts{:}];
    else
        out = cat(2, parts{:});
    end
end
