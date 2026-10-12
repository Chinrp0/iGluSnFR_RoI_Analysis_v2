function fileResult = twoap_pipeline(filepath, options, config)
    % TWOAP_PIPELINE - Process ONE 2AP (paired-pulse / PPF) raw-mean file
    %
    % Mirrors oneap_pipeline but for paired-pulse data:
    %   load raw mean -> parse ISI from filename -> derive stim2 frame ->
    %   pre-stim baseline -> dF/F -> twoap_analyzer (PPR by two methods)
    %
    % The inter-stimulus interval is read from the filename (e.g. PPF-30ms) via
    % parse_isi_ms, and the second stimulus frame is derived from it and the
    % exposure - so nothing about the ISI is hardcoded.
    %
    % Inputs:
    %   filepath - path to a single raw-mean CSV/XLSX (frames x ROIs)
    %   options  - struct, optional. Fields: verbose (logical)
    %   config   - tracenorm_config() struct, optional (pass from the batch driver
    %              to avoid reloading per file)
    %
    % Output:
    %   fileResult - struct with filename, group, ISI, dF/F data and the
    %                twoap_analyzer result.

    if nargin < 2, options = struct(); end
    if ~isfield(options, 'verbose'), options.verbose = false; end

    if nargin < 3 || isempty(config)
        config = struct();
        evalc('config = tracenorm_config();');   % suppress banner
    end
    config.dataType = '2AP';
    config.verbose  = options.verbose;

    % --- Load raw mean fluorescence ---
    loader = csv_loader_v2();
    [data, metadata] = loader.loadSingleFile(filepath, struct('dataType', 'single'));
    [numFrames, numROIs] = size(data);

    % --- ISI (from filename) -> stim2 frame ---
    isi_ms = parse_isi_ms(filepath, config);
    stim1  = config.twoAP.stim1_frame;
    stim2  = stim1 + round(isi_ms / config.exposure_ms);
    if stim2 >= numFrames
        error('twoap_pipeline:stim2OutOfRange', ...
            '%s: derived stim2 frame %d >= numFrames %d (ISI %g ms).', ...
            metadata.filename, stim2, numFrames, isi_ms);
    end
    config.twoAP.isi_ms      = isi_ms;
    config.twoAP.stim2_frame = stim2;

    if numFrames ~= config.twoAP.expected_frames && options.verbose
        fprintf('    NOTE: %s has %d frames (expected %d for 2AP)\n', ...
            metadata.filename, numFrames, config.twoAP.expected_frames);
    end

    % --- Pre-stim baseline -> dF/F ---
    [baseline, baseline_stats] = prestim_baseline(data, config);
    [dfof_data, dfof_stats]    = dfof_calculator(data, baseline, config);

    % --- Per-ROI 2AP analysis ---
    twoap = twoap_analyzer(dfof_data, baseline_stats, config);

    % --- Assemble file result ---
    fileResult = struct();
    fileResult.filename       = metadata.filename;
    fileResult.filepath       = filepath;
    fileResult.group          = parse_group(metadata.filename, config);
    fileResult.isi_ms         = isi_ms;
    fileResult.numFrames      = numFrames;
    fileResult.numROIs        = numROIs;
    fileResult.dfof_data      = dfof_data;
    fileResult.baseline_stats = baseline_stats;
    fileResult.dfof_stats     = dfof_stats;
    fileResult.twoap          = twoap;
    fileResult.config         = config;
end
