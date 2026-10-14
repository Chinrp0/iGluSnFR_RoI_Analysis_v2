function n = annotate_frame_time(frame_ms, figs)
    % ANNOTATE_FRAME_TIME - print the frame time on figures that involve time
    %
    % Adds "Frame time: <ms> ms (<Hz> Hz)" in the bottom-left corner of every
    % figure with a time axis or a time-derived quantity: any x/y label, axes
    % title, tiled-layout title or figure name mentioning ms, (s), time, tau,
    % interval / IEI, events/s or frames. Used by run_pipeline
    % (config.visualization.show_frame_time). Re-running replaces the label.
    %
    %   frame_ms - frame interval in ms (from Recordings.csv)
    %   figs     - figure handles (default: all open figures)
    %   n        - number of figures annotated

    if nargin < 2, figs = findall(groot, 'Type', 'figure'); end
    n = 0;
    if isempty(frame_ms) || ~isfinite(frame_ms) || frame_ms <= 0, return; end
    pat = '\<ms\>|\(s\)|time|\\tau|interval|\<IEI\>|events/s|frames';
    txt = sprintf('Frame time: %.4g ms (%.4g Hz)', frame_ms, 1000 / frame_ms);

    for f = figs(:)'
        if ~is_time_figure(f, pat), continue; end
        delete(findall(f, 'Tag', 'frame_time_note'));
        annotation(f, 'textbox', [0.005 0.002 0.4 0.03], 'String', txt, ...
            'Tag', 'frame_time_note', 'EdgeColor', 'none', 'FontSize', 8, ...
            'Color', [0.3 0.3 0.3], 'Interpreter', 'none', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
        n = n + 1;
    end
end

function tf = is_time_figure(f, pat)
    strs = {f.Name};
    for ax = findall(f, 'Type', 'axes')'
        strs = [strs, label_text(ax.XLabel), label_text(ax.YLabel), label_text(ax.Title)]; %#ok<AGROW>
    end
    for tl = findall(f, 'Type', 'tiledlayout')'
        strs = [strs, label_text(tl.Title), label_text(tl.XLabel)]; %#ok<AGROW>
    end
    tf = any(~cellfun(@isempty, regexpi(strs, pat, 'once')));
end

function s = label_text(h)
    % Label String as a single char row (multi-line labels joined).
    s = h.String;
    if iscell(s), s = strjoin(s, ' '); end
    if isstring(s), s = char(strjoin(s, ' ')); end
    if size(s, 1) > 1, s = strjoin(cellstr(s), ' '); end
    s = {char(s)};
end
