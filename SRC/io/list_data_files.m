function files = list_data_files(folder)
    % LIST_DATA_FILES - List tabular data files (.csv, .xlsx, .xls) in a folder
    %
    % Returns a dir-style struct array (fields .name, .folder, .isdir, ...)
    % so it is a drop-in replacement for dir(fullfile(folder, '*.csv')).
    %
    % Excel temporary lock files ("~$name.xlsx", created while a workbook is
    % open) and directories are excluded, and results are sorted by name for
    % deterministic ordering across CSV/XLSX mixes.

    files = dir(fullfile(folder, '*.csv'));
    files = [files; dir(fullfile(folder, '*.xlsx'))];
    files = [files; dir(fullfile(folder, '*.xls'))];

    if isempty(files)
        return;
    end

    names = {files.name};
    keep = ~startsWith(names, '~$') & ~[files.isdir];
    files = files(keep);

    if ~isempty(files)
        [~, order] = sort({files.name});
        files = files(order);
    end
end
