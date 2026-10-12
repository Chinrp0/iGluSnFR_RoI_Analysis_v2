function saved = save_open_figures(output_dir, figs)
    % SAVE_OPEN_FIGURES - Save figures as both PNG (viewing) and FIG (editing).
    %
    % Filenames are "<number>_<Name>" with the figure Name sanitised to word
    % characters and hyphens (falls back to "Figure_<n>" when Name is empty).
    % Shared by run_pipeline and the batch analyzers.
    %
    % Inputs:
    %   output_dir - destination directory (must exist)
    %   figs       - optional array of figure handles; defaults to ALL open
    %                figures (findall). Non-figure handles are ignored.
    % Output:
    %   saved      - number of figures written
    if nargin < 2 || isempty(figs)
        figs = findall(0, 'Type', 'figure');
    end
    figs = figs(isgraphics(figs, 'figure'));

    saved = 0;
    for i = 1:numel(figs)
        fig = figs(i);
        if ~isempty(fig.Name)
            name = fig.Name;
        else
            name = sprintf('Figure_%d', fig.Number);
        end
        name = strrep(name, ' ', '_');
        name = regexprep(name, '[^\w\-]', '');
        stem = sprintf('%d_%s', fig.Number, name);

        try
            saveas(fig, fullfile(output_dir, [stem '.png']));
            saveas(fig, fullfile(output_dir, [stem '.fig']));
            saved = saved + 1;
        catch ME
            fprintf('  Could not save %s: %s\n', stem, ME.message);
        end
    end
end
