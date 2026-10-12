function saved = save_open_figures(output_dir, config, figs)
    % SAVE_OPEN_FIGURES - Save figures in the configured image format (+ optional FIG).
    %
    % Format, resolution and whether to also write an editable .fig come from
    % config.output (tracenorm_config): figure_format, figure_dpi, save_fig.
    % Filenames are "<number>_<Name>" with the figure Name sanitised to word
    % characters and hyphens (falls back to "Figure_<n>" when Name is empty).
    % Shared by run_pipeline and the batch analyzers.
    %
    % Inputs:
    %   output_dir - destination directory (must exist)
    %   config     - tracenorm_config struct (uses config.output)
    %   figs       - optional array of figure handles; defaults to ALL open
    %                figures (findall). Non-figure handles are ignored.
    % Output:
    %   saved      - number of figures written
    if nargin < 3 || isempty(figs)
        figs = findall(0, 'Type', 'figure');
    end
    figs = figs(isgraphics(figs, 'figure'));

    out = config.output;
    fmt = lower(regexprep(out.figure_format, '^\.', ''));
    vector_formats = {'pdf', 'svg', 'eps'};
    raster_formats = {'png', 'tiff', 'tif', 'jpg', 'jpeg'};
    if ~any(strcmp(fmt, [vector_formats, raster_formats]))
        error('save_open_figures:format', ...
            'Unsupported config.output.figure_format "%s" (use %s).', ...
            out.figure_format, strjoin([raster_formats, vector_formats], ' | '));
    end

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
            img = fullfile(output_dir, [stem '.' fmt]);
            if any(strcmp(fmt, vector_formats))
                exportgraphics(fig, img, 'ContentType', 'vector');
            else
                exportgraphics(fig, img, 'Resolution', out.figure_dpi);
            end
            if out.save_fig
                savefig(fig, fullfile(output_dir, [stem '.fig']));
            end
            saved = saved + 1;
        catch ME
            fprintf('  Could not save %s: %s\n', stem, ME.message);
        end
    end
end
