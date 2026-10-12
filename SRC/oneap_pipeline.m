function fileResult = oneap_pipeline(filepath, options, config)
    % ONEAP_PIPELINE - Process ONE 1AP raw-mean CSV into per-ROI evoked metrics
    %
    % Mirrors the single-file path of main_pipeline but for evoked (1AP) data:
    %   load raw mean -> pre-stim baseline -> dF/F -> oneap_analyzer
    %
    % Inputs:
    %   filepath - path to a single raw-mean CSV (frames x ROIs)
    %   options  - struct, optional. Fields: verbose (logical)
    %   config   - tracenorm_config() struct, optional. Pass it in (e.g. from
    %              oneap_batch_analysis) to avoid reloading config per file.
    %
    % Output:
    %   fileResult - struct with filename, genotype, dF/F data, baseline stats
    %                and the oneap_analyzer result.

    if nargin < 2, options = struct(); end
    if ~isfield(options, 'verbose'), options.verbose = false; end

    if nargin < 3 || isempty(config)
        config = struct();
        evalc('config = tracenorm_config();');   % suppress banner
    end
    config.dataType = '1AP';
    config.verbose = options.verbose;

    % --- Load raw mean fluorescence ---
    % The loader is layout-aware (skips header rows, strips the frame-index
    % column only when present), so it handles both the macro CSVs and the older
    % 1AP .xlsx exports without any per-type flags.
    loader = csv_loader_v2();
    [data, metadata] = loader.loadSingleFile(filepath, struct('dataType', 'single'));

    [numFrames, numROIs] = size(data);
    if numFrames ~= config.oneAP.expected_frames && options.verbose
        fprintf('    NOTE: %s has %d frames (expected %d for 1AP)\n', ...
            metadata.filename, numFrames, config.oneAP.expected_frames);
    end

    % --- Pre-stim baseline -> dF/F ---
    [baseline, baseline_stats] = prestim_baseline(data, config);
    [dfof_data, dfof_stats]    = dfof_calculator(data, baseline, config);

    % --- Per-ROI 1AP analysis ---
    oneap = oneap_analyzer(dfof_data, baseline_stats, config);

    % --- Assemble file result ---
    fileResult = struct();
    fileResult.filename       = metadata.filename;
    fileResult.filepath       = filepath;
    fileResult.genotype       = parse_group(metadata.filename, config);
    fileResult.numFrames      = numFrames;
    fileResult.numROIs        = numROIs;
    fileResult.dfof_data      = dfof_data;
    fileResult.baseline_stats = baseline_stats;
    fileResult.dfof_stats     = dfof_stats;
    fileResult.oneap          = oneap;
    fileResult.config         = config;
end
