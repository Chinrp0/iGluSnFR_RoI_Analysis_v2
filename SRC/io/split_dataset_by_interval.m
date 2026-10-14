function parts = split_dataset_by_interval(dataset, config)
    % SPLIT_DATASET_BY_INTERVAL - One sub-dataset per frame interval
    %
    % Groups are pooled onto one time axis, so files recorded at different frame
    % intervals (e.g. 5 ms and 10 ms PPF in one folder) are analysed as separate
    % sets rather than resampled. Intervals within config.metadata.interval_tol_frac
    % of each other count as the same (10.0075 vs 10.0080 ms).
    %
    % Output: struct array of datasets (same fields as the input), main set
    % (most files, then shortest interval) LAST, each with
    %   .files              - the files at that interval
    %   .frame_interval_ms  - the interval (median of its files)
    %   .label              - e.g. '10ms', used for the output sub-folder

    fi  = arrayfun(@(f) f.meta.frame_interval_ms, dataset.files);
    tol = config.metadata.interval_tol_frac;
    [~, order] = sort(fi);
    cluster = zeros(size(fi));
    c = 0; ref = NaN;
    for k = order(:)'
        if isnan(ref) || abs(fi(k) - ref) > tol * ref
            c = c + 1; ref = fi(k);
        end
        cluster(k) = c;
    end

    parts = repmat(dataset, 1, c);
    for k = 1:c
        m = cluster == k;
        parts(k).files = dataset.files(m);
        parts(k).frame_interval_ms = median(fi(m));
        parts(k).label = sprintf('%gms', round(median(fi(m)), 2));
    end
    % Run order: the MAIN set last, because run_pipeline keeps only the last
    % set's figures open. Main = most files; ties go to the shorter interval
    % (e.g. a folder of 5 ms files plus one 25 ms file shows the 5 ms set).
    nf = arrayfun(@(p) numel(p.files), parts);
    [~, order] = sortrows([nf(:), -[parts.frame_interval_ms]']);
    parts = parts(order);
    if c > 1
        fprintf('Frame intervals differ: analysing %d sets separately (%s)\n', c, ...
            strjoin(arrayfun(@(p) sprintf('%s: %d files', p.label, numel(p.files)), ...
            parts, 'UniformOutput', false), ', '));
    end
end
