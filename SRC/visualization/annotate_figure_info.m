function n = annotate_figure_info(info, frame_ms, figs)
    % ANNOTATE_FIGURE_INFO - stamp sensor, experiment date and frame time on figures
    %
    % Adds e.g. "Snfr4f | 2026-10-06 | Frame time: 10 ms (100 Hz)" in the
    % bottom-left corner. Sensor and date go on every figure; the frame time
    % only on figures with a time axis or a time-derived quantity: any x/y
    % label, axes title, tiled-layout title or figure name mentioning ms, (s),
    % time, tau, interval / IEI, events/s or frames. Used by run_pipeline
    % (config.visualization.show_experiment_info / show_frame_time).
    % Re-running replaces the label.
    %
    %   info     - struct with .sensor and .date (experiment_info); fields that
    %              are '' are left out. Pass [] for none.
    %   frame_ms - frame interval in ms (from Recordings.csv); [] / NaN for none
    %   figs     - figure handles (default: all open figures)
    %   n        - number of figures annotated

    if nargin < 3, figs = findall(groot, 'Type', 'figure'); end
    n = 0;
    head = {};
    if isstruct(info)
        head = {info.sensor, info.date};
        head = head(~cellfun(@isempty, head));
    end
    time_txt = '';
    if ~isempty(frame_ms) && isfinite(frame_ms) && frame_ms > 0
        time_txt = sprintf('Frame time: %.4g ms (%.4g Hz)', frame_ms, 1000 / frame_ms);
    end
    pat = '\<ms\>|\(s\)|time|\\tau|interval|\<IEI\>|events/s|frames';

    for f = figs(:)'
        items = head;
        if ~isempty(time_txt) && is_time_figure(f, pat), items{end+1} = time_txt; end %#ok<AGROW>
        if isempty(items), continue; end
        txt = strjoin(items, ' | ');
        delete(findall(f, 'Tag', 'frame_time_note'));
        reserve_bottom_strip(f);
        annotation(f, 'textbox', [0.005 0.002 0.6 0.03], 'String', txt, ...
            'Tag', 'frame_time_note', 'EdgeColor', 'none', 'FontSize', 8, ...
            'Color', [0.3 0.3 0.3], 'Interpreter', 'none', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
        n = n + 1;
    end
end

function reserve_bottom_strip(f)
    % Tiled layouts with tight padding run to the figure's bottom edge, so the
    % stamp would sit on the bottom x tick labels. Lift such layouts by a
    % ~20 px strip (subplot figures already have a bottom margin).
    hpx = getpixelposition(f); hpx = hpx(4);
    strip = min(0.08, 20 / max(hpx, 1));
    for tl = findall(f, 'Type', 'tiledlayout', 'Parent', f)'
        op = tl.OuterPosition;   % normalized units
        if op(2) < strip
            tl.OuterPosition = [op(1), strip, op(3), op(4) - (strip - op(2))];
        end
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
