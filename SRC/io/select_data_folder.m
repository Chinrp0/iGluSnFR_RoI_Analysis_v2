function folder = select_data_folder(config)
    % SELECT_DATA_FOLDER - Ask for the raw-mean data folder with a folder dialog
    %
    % The dialog opens at config.data_root (or the current folder when that path
    % is not available on this PC). Cancelling stops the run with an error rather
    % than continuing on some default folder.
    %
    % Without a display (e.g. matlab -batch) there is no dialog: set the folder
    % in a variable first instead,  data_folder = 'E:\...'; load_dataset

    start = pwd;
    if isfield(config, 'data_root') && ~isempty(config.data_root) && isfolder(config.data_root)
        start = config.data_root;
    end

    try
        folder = uigetdir(start, 'Select the raw-mean data folder (contains Recordings.csv)');
    catch ME
        error('select_data_folder:noDisplay', ...
            ['Cannot open the folder dialog (%s)\nSet the folder first, e.g.\n' ...
             '  data_folder = ''E:\\Data\\GluSnFR\\Ms\\...''; load_dataset'], ME.message);
    end
    if isequal(folder, 0)
        error('select_data_folder:cancelled', 'No data folder selected; run cancelled.');
    end
    fprintf('Data folder: %s\n', folder);
end
