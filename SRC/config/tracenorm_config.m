function config = tracenorm_config()
    % TRACENORM_CONFIG - SINGLE SOURCE OF TRUTH for the iGluSnFR analysis pipeline
    %
    % Everything time-related is derived from ONE acquisition knob, exposure_ms
    % (milliseconds per frame). Frame rate, all ms<->frame conversions, rolling
    % windows and figure time axes follow from it. Nothing stores Hz or a frame
    % count independently, so there is no hidden constant that can shadow this.
    %
    % >>> EDIT THESE TWO VALUES PER RUN (a folder = one paradigm) <<<
    %   exposure_ms : 1AP / 2AP (PPF) = 5 ms   |   spontaneous = 25 ms
    %   dataType    : '1AP' | '2AP' | 'Spont'
    % These are set by hand for each run; the pipeline never auto-detects them.
    %
    % Groups (WT vs R213W, or future proteins) are parsed from the filename using
    % the patterns in config.groups via parse_group(). Swapping proteins is a
    % one-line edit here and nowhere else.

    config = struct();

    %% ======================================================================
    %% ACQUISITION - the single time knob + paradigm selector (EDIT PER RUN)
    %% ======================================================================
    config.exposure_ms = 5;      % ms per frame. 1AP/2AP = 5 ms; spont = 25 ms.
    config.dataType    = '2AP';  % '1AP' | '2AP' | 'Spont'

    % Data folder for THIS run (the paradigm's raw-mean folder). Edit per run,
    % together with dataType/exposure above. run_pipeline_fixed reads this.
    %config.folder_path = 'E:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglu3fast_NGR\1AP\GPU_Processed_Images_1AP\5_raw_mean'; 
    config.folder_path = 'E:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglu3fast_NGR\PPF\PPF_2026-06-30--12\E_raw_mean';
    %config.folder_path = 'E:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglusnfr4f_NGR\1AP\1AP_2026-07-07--1\E_RoI_Raw-Mean';
    % 1AP:   E:\...\iglu3fast_NGR\1AP\GPU_Processed_Images_1AP\5_raw_mean
    % Spont: E:\...\iglusnfr4f_NGR\Spont\<run>\E_raw_mean

    % --- Everything below is DERIVED from exposure_ms (do not hardcode Hz) ---
    config.frame_rate     = 1000 / config.exposure_ms;  % Hz, derived
    config.frame_dur_ms   = config.exposure_ms;         % ms/frame, alias for readability

    %% ======================================================================
    %% GROUPS - config-driven, fail-loud parsing (see io/parse_group.m)
    %% ======================================================================
    % name    : label used in legends/titles/aggregation
    % pattern : case-insensitive regexp matched against the filename stem
    % color   : RGB for all figures for this group
    % To add/swap a protein, add/edit one row here. A filename that matches no
    % pattern (or more than one) makes parse_group error rather than mislabel.
    config.groups = struct( ...
        'name',    {'WT',            'R213W'}, ...
        'pattern', {'WT',            'R213W'}, ...
        'color',   {[0.2 0.6 1.0],   [1.0 0.4 0.2]});

    %% ======================================================================
    %% 1AP (single evoked action potential)
    %% ======================================================================
    % All windows are in MILLISECONDS relative to stimulus onset and converted to
    % frames at run time via exposure_ms. Only stim_frame / expected_frames are in
    % frames because they are properties of the specific recording, not the rate.
    config.oneAP = struct();
    config.oneAP.expected_frames    = 600;      % Frames per 1AP trace (this dataset)
    config.oneAP.stim_frame         = 266;      % Frame of stimulus onset (data-confirmed)
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
    config.twoAP.auc_window_ms      = 40;       % AUC PPR: integrate dF/F over [0, this] ms after each stim.
                                                % Both pulses use the SAME window = min(this, ISI) so A1/A2
                                                % are comparable. A1 above pre-stim1 baseline, A2 above the
                                                % local pre-stim2 residual (same baselines as the peak methods).
    config.twoAP.montage_pre_ms     = 250;      % Montage x-window: ms of trace shown BEFORE stim1 (to judge
    config.twoAP.montage_post_ms    = 300;      % baseline/noise) and AFTER stim2 (to see the full decay).
    config.twoAP.sync_latency_ms    = 10;       % "Sync time-lock" PPR: read BOTH pulses at this fixed
                                                % latency post-stim (the synchronous-release phase, before
                                                % the peak). See the three PPR methods in twoap_analyzer.
    config.twoAP.baseline_guard_ms  = 0;        % Time before stim1 dropped from baseline
    config.twoAP.tau_peak_window_ms = 50;       % Decay-onset peak window (post-stim2)
    config.twoAP.tau_min_r2         = 0.5;
    config.twoAP.upper_threshold_sigma  = 3.0;  % Responder threshold (peak must clear this)
    config.twoAP.lower_threshold_sigma  = 1.5;  % Schmitt decay threshold (hysteresis)
    config.twoAP.min_frames_above_lower = 2;    % Schmitt validation (SAME as 1AP): after the peak
                                                % clears the UPPER threshold, the run of frames >= the
                                                % LOWER threshold spanning it must be >= this. Keeps fast
                                                % transients that clear the upper for one frame, rejects
                                                % single-frame noise spikes.
    config.twoAP.max_noise_sd           = 0.02;

    %% ======================================================================
    %% Mode-dependent frame count / duration (DERIVED, never a hidden constant)
    %% ======================================================================
    switch upper(config.dataType)
        case '1AP',   config.expected_frames = config.oneAP.expected_frames;
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
    config.visualization.x_axis.montage_range_frames = [100, 500];   % per-panel ROI montage window

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
    fprintf('Exposure: %.3g ms/frame  ->  Frame Rate: %.4g Hz\n', ...
        config.exposure_ms, config.frame_rate);
    fprintf('Recording: %d frames = %.2f s\n', config.expected_frames, config.recording_duration_s);
    fprintf('Groups: %s\n', strjoin({config.groups.name}, ', '));
    if strcmpi(config.dataType, '1AP')
        fprintf('1AP Stim Frame: %d (%.0f ms), Sync: %d-%d ms, Async: %d-%d ms, Bin: %d ms\n', ...
            config.oneAP.stim_frame, (config.oneAP.stim_frame - 1) * config.frame_dur_ms, ...
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
