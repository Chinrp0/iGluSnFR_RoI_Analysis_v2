function output_dir = make_figures_dir(folder_path)
    % MAKE_FIGURES_DIR - Dated, auto-incremented output directory for a run.
    %
    % Creates "F_Analysis_<yyyy-mm-dd>--<run>" one level ABOVE the data folder,
    % i.e. next to the folder holding the _mean.csv files:
    %   ...\5AP_2026-10-11--1\E_RoI_Raw-Mean\         <- data (folder_path)
    %   ...\5AP_2026-10-11--1\F_Analysis_2026-10-12--1 <- output
    % <run> increments so repeated runs on the same day never collide.
    % Shared by run_pipeline and the batch analyzers so every entry point
    % uses the same output-folder convention.
    %
    % Input:  folder_path - the data folder that was analysed
    % Output: output_dir  - absolute path of the created directory
    parent = fileparts(strip(char(folder_path), 'right', filesep));

    base = sprintf('F_Analysis_%s', datestr(now, 'yyyy-mm-dd'));

    run_num = 1;
    while exist(fullfile(parent, sprintf('%s--%d', base, run_num)), 'dir')
        run_num = run_num + 1;
    end

    output_dir = fullfile(parent, sprintf('%s--%d', base, run_num));
    mkdir(output_dir);
end
