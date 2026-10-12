function config = tracenorm_config()
    % TRACENORM_CONFIG - SINGLE SOURCE OF TRUTH for the iGluSnFR analysis pipeline
    %
    % RECORDING TIMING COMES FROM THE METADATA, NOT FROM HERE. Every data folder
    % carries a Recordings.csv (written by the Fiji macro) giving, per file, the
    % number of frames, the frame interval (ms) and the refined stimulus (TTL)
    % frames. load_dataset_folder() reads it and apply_recording_timing() writes
    % those values into the config per file/run, so the frame rate, frame count
    % and stim frame are never hardcoded for 1AP. All analysis windows below are
    % in MILLISECONDS and converted to frames with the metadata frame interval.
    %
    % >>> EDIT PER RUN <<<
    %   dataType    : '1AP' | '2AP' | 'Spont'
    % The data folder is NOT set here: load_dataset / run_pipeline ask for
    % it with a folder dialog that opens at data_root (select_data_folder).
    %
    % Groups are genotype x drug. Genotype is parsed from the filename with the
    % patterns in config.groups (parse_group); drug from the "+drug" token after
    % the paradigm (parse_drug), e.g. "..._5AP+ttx-01" -> ttx, none -> control.

    config = struct();

    %% ======================================================================
    %% ACQUISITION - paradigm selector (EDIT PER RUN) + folder-dialog start
    %% ======================================================================
    config.dataType    = '1AP';  % '1AP' | '2AP' | 'Spont'

    % Where the folder dialog opens (falls back to the current folder if this
    % drive/path is not available on this PC).
    config.data_root   = 'E:\Data\GluSnFR\Ms\';

    % LEGACY exposure, used ONLY by the 2AP and Spont paths, which have not been
    % moved onto Recordings.csv yet. The 1AP path ignores it: apply_recording_timing
    % overwrites exposure_ms / frame_rate / frame_dur_ms from the metadata.
    config.exposure_ms = 5;      % ms per frame (2AP = 5 ms; spont = 25 ms)

    % --- Derived from exposure_ms (replaced per run by apply_recording_timing) ---
    config.frame_rate     = 1000 / config.exposure_ms;  % Hz, derived
    config.frame_dur_ms   = config.exposure_ms;         % ms/frame, alias for readability

    %% ======================================================================
    %% METADATA (Recordings.csv, see io/read_recordings_metadata.m)
    %% ======================================================================
    config.metadata = struct();
    config.metadata.filename = 'Recordings.csv';   % matched case-insensitively
    % Frames per file must equal N_Frames, and ROI columns must equal N_ROIs.
    % A file that disagrees is skipped (with a warning) rather than analysed on
    % a timing it doesn't have.
    config.metadata.interval_tol_frac = 0.001;     % Files pooled in one run must
                                                   % share Frame_Interval_ms to
                                                   % within this fraction.

    %% ======================================================================
    %% GROUPS - config-driven, fail-loud parsing (see io/parse_group.m)
    %% ======================================================================
    % name    : genotype label used in legends/titles/aggregation
    % pattern : case-insensitive regexp matched against the filename stem
    % color   : RGB for all figures for this genotype (drugs are shades of it)
    % Patterns must be mutually exclusive: "Doc2b-WT" must NOT also match plain
    % "WT", so the plain genotypes only match as a whole "_WT_" token.
    %   CP_Snfr4_WT_Cs2-c1_5AP-02            -> WT
    %   CP_Snfr4_Doc2b-WT_10ms_Cs1-c1_5AP-01 -> Doc2b-WT
    % A filename that matches no pattern (or more than one) makes parse_group
    % error rather than mislabel.
    config.groups = struct( ...
        'name',    {'WT',              'Doc2b-WT',      'Doc2b-R213W',      'R213W'}, ...
        'pattern', {'(^|_)WT(_|$)',    'Doc2b[-_]WT',   'Doc2b[-_]R213W',   '(^|_)R213W(_|$)'}, ...
        'color',   {[0.3 0.3 0.3],     [0.2 0.6 1.0],   [1.0 0.4 0.2],      [0.8 0.2 0.6]});

    % Drug conditions (perfusion). parse_drug returns the "+token" after the
    % paradigm in lower case, or control_label when there is none. drug_order
    % sets the plotting order; drugs not listed here are appended after it.
    config.drugs = struct();
    config.drugs.control_label = 'control';
    config.drugs.order         = {'control', 'egta', 'ttx'};

    %% ======================================================================
    %% 1AP (single evoked action potential, one or more trials per file)
    %% ======================================================================
    % All windows are in MILLISECONDS relative to stimulus onset and converted to
    % frames with the metadata frame interval. The stim frames themselves come
    % from Recordings.csv (Stim_Frames_Refined, one per trial). Each trial is cut
    % into its own epoch [stim - epoch_pre_ms, stim + epoch_post_ms] (clipped so
    % it never reaches into the neighbouring trial) and analysed on its own, so a
    % failed trial at a low-release-probability site stays a failure instead of
    % being averaged away. oneAP.stim_frame is set at run time to the stim index
    % WITHIN the epoch; it is not a property of the recording.
    config.oneAP = struct();
    config.oneAP.epoch_pre_ms       = 1000;     % Pre-stim baseline kept per trial epoch
    config.oneAP.epoch_post_ms      = 1500;     % Post-stim trace kept per trial epoch
    config.oneAP.sync_window_ms     = [0, 15];  % Synchronous release window (timing label)
    config.oneAP.async_window_ms    = [15, 250];% Asynchronous release window (timing label)
                                                % Evoked iGlu peak lands ~15-35 ms post-stim
                                                % (median ~25 ms), so the responder call and
                                                % amplitude use the LARGEST peak, not the sync window.
    config.oneAP.bin_size_ms        = 5;        % Temporal (PSTH) histogram bin
    config.oneAP.auc_window_ms      = [0, 250]; % Area-under-curve window (baseline-subtracted)
    config.oneAP.tau_peak_window_ms = 50;       % Decay-onset peak searched in [0, this] post-stim
    config.oneAP.tau_min_r2         = 0.5;      % Min fit R^2 to keep a decay tau
    config.oneAP.baseline_guard_ms  = 0;        % Time before stim dropped from baseline (settling)
    config.oneAP.upper_threshold_sigma  = 3.0;  % Event / responder threshold (peak height)
    config.oneAP.lower_threshold_sigma  = 1.5;  % Lower threshold (reference / decay end)
    config.oneAP.min_peak_distance_ms   = 10;   % Min spacing between detected peaks (findpeaks)
    config.oneAP.min_frames_above_lower = 2;    % Schmitt validation: after a peak crosses the UPPER
                                                % (3σ) threshold it must stay above the LOWER (1.5σ)
                                                % threshold for at least this many consecutive frames
                                                % (i.e. show a real decay). Rejects single-frame noise
                                                % spikes WITHOUT discarding fast real transients that
                                                % clear 3σ for only one frame. 1 = height-only.
    config.oneAP.max_noise_sd           = 0.02; % Drop ROIs whose pre-stim robust SD exceeds this.

    %% ======================================================================
    %% 2AP (paired-pulse facilitation, PPF)
    %% ======================================================================
    % The inter-pulse interval (ISI) is encoded in the filename (e.g. PPF-30ms),
    % parsed by parse_isi_ms() using isi_pattern. stim2_frame is derived per file:
    %   stim2_frame = stim1_frame + round(ISI_ms / exposure_ms)
    config.twoAP = struct();
    config.twoAP.expected_frames    = 600;      % Frames per PPF trace (this dataset)
    config.twoAP.stim1_frame        = 266;      % First pulse frame (shared with 1AP, data-confirmed)
    config.twoAP.isi_pattern        = 'PPF[\s_\-]*0*(\d+)\s*ms';  % -> ISI in ms from filename
    config.twoAP.peak_search_ms     = 40;       % Search [0, this] ms after each stim for that pulse's peak
                                                % (capped at the ISI for pulse 1 so it can't reach pulse 2).
    config.twoAP.auc_post_ms        = 250;      % AUC = TOTAL CHARGE over ONE window that encompasses BOTH
                                                % pulses: [stim1, stim2 + this]. A single pre-stim1 baseline
                                                % is used throughout (no per-pulse residual), so the metric is
                                                % the total glutamate charge of the paired response (dF/F*ms),
                                                % NOT a ratio. 250 ms matches the 1AP auc window and captures
                                                % the full iGluSnFR3fast decay after pulse 2 while staying
                                                % inside the 600-frame trace at the longest ISI (stim2~=466 ->
                                                % +50 frames = 516). NB the window length = ISI + this, so it
                                                % grows with ISI by construction: raw charge is ISI-dependent
                                                % (interpret cross-ISI accordingly; no normalization applied).
    config.twoAP.montage_pre_ms     = 250;      % Montage x-window: ms of trace shown BEFORE stim1 (to judge
    config.twoAP.montage_post_ms    = 300;      % baseline/noise) and AFTER stim2 (to see the full decay).
    config.twoAP.sync_latency_ms    = 10;       % "Sync time-lock" PPR: read BOTH pulses at this fixed
                                                % latency post-stim (the synchronous-release phase, before
                                                % the peak). See the three PPR methods in twoap_analyzer.
    config.twoAP.baseline_guard_ms  = 0;        % Time before stim1 dropped from baseline
    config.twoAP.tau_peak_window_ms = 50;       % Decay-onset peak window (post-stim2)
    config.twoAP.tau_min_r2         = 0.5;
    config.twoAP.upper_threshold_sigma  = 3.5;  % Responder threshold (peak must clear this)
    config.twoAP.lower_threshold_sigma  = 1.5;  % Schmitt decay threshold (hysteresis)
    config.twoAP.min_frames_above_lower = 2;    % Schmitt validation (SAME as 1AP): after the peak
                                                % clears the UPPER threshold, the run of frames >= the
                                                % LOWER threshold spanning it must be >= this. Keeps fast
                                                % transients that clear the upper for one frame, rejects
                                                % single-frame noise spikes.
    config.twoAP.max_noise_sd           = 0.02;
    config.twoAP.min_peak_distance_ms   = 10;   % findpeaks min spacing (SAME as 1AP). A pulse peak must
                                                % be a genuine LOCAL MAXIMUM >= upper (via findpeaks),
                                                % not just the window max, so a monotonic noise ramp no
                                                % longer counts as a response.
    config.twoAP.baseline_excursion_pct = 99;   % Noise-ceiling QC: the pulse-1 peak must exceed the
    config.twoAP.peak_over_baseline_k   = 1.0;  % k * (this percentile of pre-stim1 excursions above
                                                % baseline). Rejects ROIs whose "response" is no bigger
                                                % than their own largest baseline swing (uniform-noise
                                                % ROIs that clear 3-sigma by chance). k=1 => peak must beat
                                                % its 99th-pct baseline excursion. Set k=0 to disable.

    %% ======================================================================
    %% Mode-dependent frame count / duration (DERIVED, never a hidden constant)
    %% ======================================================================
    % 1AP: NaN until apply_recording_timing sets it from Recordings.csv N_Frames.
    switch upper(config.dataType)
        case '1AP',   config.expected_frames = NaN;
        case '2AP',   config.expected_frames = config.twoAP.expected_frames;
        case 'SPONT', config.expected_frames = 1200;   % 30 s at 25 ms/frame (this dataset)
        otherwise,    error('tracenorm_config:dataType', ...
                          'Unknown dataType "%s" (expected 1AP | 2AP | Spont).', config.dataType);
    end
    config.recording_duration_s = config.expected_frames / config.frame_rate;
    config.min_valid_frames     = 400;          % Minimum frames required for baseline calc

    %% ======================================================================
    %% Baseline (spontaneous rolling-median path)
    %% ======================================================================
    config.baseline_method       = 'iterative_rolling_median';
    config.rolling_window_sec    = 0.50;        % Rolling window in seconds
    config.rolling_window_frames = round(config.rolling_window_sec * config.frame_rate);
    config.outlier_threshold_sigma = 2.0;       % SDs for outlier detection
    config.max_iterations        = 3;           % Refinement iterations

    config.verbose               = false;       % Suppress detailed console output
    config.smooth_baseline       = false;       % Don't smooth final baseline
    config.smooth_method         = 'movmedian';
    config.smooth_window         = 5;           % frames
    config.transport_slope_threshold = 0.01;    % Baseline drift flag, in F units per SECOND
                                                % (frame-rate independent; see baseline_detector).
    config.min_baseline_frames   = 400;

    %% ======================================================================
    %% dF/F
    %% ======================================================================
    config.dfof_method       = 'divide';        % (F - F0) / F0
    config.min_baseline_value = 0.01;           % Avoid division by zero

    %% ======================================================================
    %% Corrected Schmitt trigger (spontaneous event detection)
    %% ======================================================================
    config.corrected_schmitt = struct();
    config.corrected_schmitt.upper_threshold_sigma = 3.5;  % Event start
    config.corrected_schmitt.lower_threshold_sigma = 1.5;  % Event end
    config.corrected_schmitt.min_event_duration    = 3;    % frames (validation)
    config.corrected_schmitt.merge_gap_frames      = 2;    % Merge events <= this many frames apart
    config.corrected_schmitt.noise_exclusion_window = 7;   % Exclude sustained events >= this (frames)
    config.corrected_schmitt.sustained_percentile  = 85;   % Sustained-event percentile threshold
    config.corrected_schmitt.includes_outliers_in_noise = true;
    config.corrected_schmitt.excludes_sustained_events  = true;

    %% === Legacy Schmitt trigger (kept for comparison / pure detector) ===
    config.event_detection = struct();
    config.event_detection.method = 'schmitt_trigger';
    config.event_detection.upper_threshold_sigma = 3.0;
    config.event_detection.lower_threshold_sigma = 1.5;
    config.event_detection.min_event_duration    = 3;      % frames
    config.event_detection.merge_gap_frames      = 2;      % frames

    %% ======================================================================
    %% Quality metrics
    %% ======================================================================
    config.quality = struct();
    config.quality.min_baseline_stability = 0.8;
    config.quality.max_noise_level        = 0.2;
    config.quality.min_snr                = 2.0;

    %% ======================================================================
    %% Temporal (all DERIVED from exposure_ms; no independent time base stored)
    %% ======================================================================
    config.temporal = struct();
    config.temporal.frame_duration_ms    = config.frame_dur_ms;          % = exposure_ms
    config.temporal.recording_duration_s = config.recording_duration_s;
    config.temporal.max_iei_s            = 10;    % Max IEI to analyze (s)
    config.temporal.bin_size_ms          = 50;    % IEI histogram bins (ms)

    %% ======================================================================
    %% Visualization
    %% ======================================================================
    config.visualization = struct();
    config.visualization.plot_baseline = true;
    config.visualization.plot_events   = true;
    config.visualization.max_traces_to_plot = 100;
    % Legacy convenience aliases; group colors are the source of truth in
    % config.groups (use color_for/parse_group). Default to the WT/R213W groups,
    % falling back to the first two configured groups so older plotters keep
    % working even if the group names change.
    wtIdx  = find(strcmpi({config.groups.name}, 'WT'),    1);
    mutIdx = find(strcmpi({config.groups.name}, 'R213W'), 1);
    if isempty(wtIdx),  wtIdx  = 1; end
    if isempty(mutIdx), mutIdx = min(2, numel(config.groups)); end
    config.visualization.colors = struct( ...
        'wt',  config.groups(wtIdx).color, ...
        'mut', config.groups(mutIdx).color);

    % --- Figure x-axis display units (time-based figures) ---------------------
    % Every time-from-stimulus axis (1AP/2AP average traces, PSTH, ROI
    % montages) and the QC peak-time histogram read their unit and tick spacing
    % from here. Traces are computed in frames/ms internally and only CONVERTED
    % for display, so detection, stats and tables are unaffected by this choice.
    % See visualization/make_time_axis.m + format_time_axis.m.
    config.visualization.x_axis = struct();
    config.visualization.x_axis.unit = 'ms';   % 'frames' | 'ms' | 's'  (axis LABELS only)
    % Plotted x-RANGE as an absolute frame interval [start end] ([] = auto per
    % figure). Frames are a property of the recording (like stim_frame), so the
    % range is given in frames regardless of the display unit above; it is
    % converted to the chosen unit for display. Applies to the dF/F trace
    % figures (1AP/2AP average trace + ROI montages). PSTH / QC histograms keep
    % their own natural domain (event / peak time).
    %   e.g. [0 600] = whole trace   |   [200 500] = zoom around the stim
    config.visualization.x_axis.range_frames         = [];   % avg-trace (wide time figure)
    % For 1AP these are frames WITHIN a trial epoch (stim at epoch_pre_ms).
    config.visualization.x_axis.montage_range_frames = [];   % per-panel ROI montage window

    %% ======================================================================
    %% Output
    %% ======================================================================
    config.output = struct();
    config.output.save_figures  = false;
    config.output.figure_format = 'png';
    config.output.figure_dpi    = 300;
    config.output.save_data     = false;

    %% ======================================================================
    %% Summary
    %% ======================================================================
    fprintf('\n=== Configuration Loaded ===\n');
    fprintf('Analysis Mode: %s\n', config.dataType);
    if strcmpi(config.dataType, '1AP')
        fprintf('Timing: per file from %s (frames, frame interval, stim frames)\n', ...
            config.metadata.filename);
    else
        fprintf('Exposure: %.3g ms/frame  ->  Frame Rate: %.4g Hz (legacy config value)\n', ...
            config.exposure_ms, config.frame_rate);
        fprintf('Recording: %d frames = %.2f s\n', config.expected_frames, config.recording_duration_s);
    end
    fprintf('Genotypes: %s | Drugs: %s\n', strjoin({config.groups.name}, ', '), ...
        strjoin(config.drugs.order, ', '));
    if strcmpi(config.dataType, '1AP')
        fprintf('1AP epochs: -%g..+%g ms per trial, Sync: %d-%d ms, Async: %d-%d ms, Bin: %d ms\n', ...
            config.oneAP.epoch_pre_ms, config.oneAP.epoch_post_ms, ...
            config.oneAP.sync_window_ms(1), config.oneAP.sync_window_ms(2), ...
            config.oneAP.async_window_ms(1), config.oneAP.async_window_ms(2), ...
            config.oneAP.bin_size_ms);
    elseif strcmpi(config.dataType, '2AP')
        fprintf('2AP Stim1 Frame: %d (%.0f ms); ISI parsed per file; peak search %.0f ms\n', ...
            config.twoAP.stim1_frame, (config.twoAP.stim1_frame - 1) * config.frame_dur_ms, ...
            config.twoAP.peak_search_ms);
    else
        fprintf('Baseline Window: %.0f ms (%d frames); Corrected Schmitt %.1f/%.1f sigma\n', ...
            config.rolling_window_sec * 1000, config.rolling_window_frames, ...
            config.corrected_schmitt.upper_threshold_sigma, ...
            config.corrected_schmitt.lower_threshold_sigma);
    end
    fprintf('============================\n\n');
end
