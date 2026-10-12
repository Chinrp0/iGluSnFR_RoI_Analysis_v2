function output_dir = make_figures_dir(folder_path)
    % MAKE_FIGURES_DIR - Dated, auto-incremented output directory for figures.
    %
    % Creates "<data_name>_figures_<yyyymmdd>_<run>" one level ABOVE the data
    % folder, incrementing <run> so repeated runs on the same day never collide.
    % Shared by run_pipeline and the batch analyzers so every entry point
    % uses the same output-folder convention.
    %
    % Input:  folder_path - the data folder that was analysed
    % Output: output_dir  - absolute path of the created directory
    [parent, data_name] = fileparts(folder_path);
    [grandparent, ~]     = fileparts(parent);

    date_str = datestr(now, 'yyyymmdd');
    base     = sprintf('%s_figures_%s', data_name, date_str);

    run_num = 1;
    while exist(fullfile(grandparent, sprintf('%s_%d', base, run_num)), 'dir')
        run_num = run_num + 1;
    end

    output_dir = fullfile(grandparent, sprintf('%s_%d', base, run_num));
    mkdir(output_dir);
end
