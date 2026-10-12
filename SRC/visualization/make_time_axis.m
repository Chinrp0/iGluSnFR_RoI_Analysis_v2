function tax = make_time_axis(config)
    % MAKE_TIME_AXIS - display-unit helper for figure time axes.
    %
    % Single place that turns the user's choice in
    %   config.visualization.x_axis.unit  ('frames' | 'ms' | 's')
    % into the conversions every plot needs. Traces are computed in
    % frames / ms internally; only the DISPLAY is converted, so nothing
    % downstream (detection, stats, tables) changes with the unit.
    %
    % Returns a struct:
    %   .unit                  - unit label for axis text ('frames' | 'ms' | 's')
    %   .scale                 - multiply an ms value to get the display value
    %   .frame_dur_ms          - ms per frame (from config)
    %   .from_ms(ms)           - ms            -> display units
    %   .from_frames_rel(nf)   - frames-from-a-reference -> display units
    %   .range_frames          - absolute [start end] frame range or []
    %   .montage_range_frames  - absolute [start end] frame range for montages or []
    %
    % See format_time_axis() for labelling an axis and applying the range.

    frame_dur = config.frame_dur_ms;   % ms per frame (single time knob)

    unit = 'ms'; rng = []; mrng = [];
    if isfield(config, 'visualization') && isfield(config.visualization, 'x_axis')
        xa = config.visualization.x_axis;
        if isfield(xa, 'unit') && ~isempty(xa.unit),   unit = xa.unit;  end
        if isfield(xa, 'range_frames'),                 rng  = xa.range_frames;  end
        if isfield(xa, 'montage_range_frames'),          mrng = xa.montage_range_frames; end
    end

    switch lower(strtrim(unit))
        case {'frame', 'frames'}
            scale = 1 / frame_dur; ulabel = 'frames';
        case {'ms', 'millisecond', 'milliseconds'}
            scale = 1;             ulabel = 'ms';
        case {'s', 'sec', 'secs', 'second', 'seconds'}
            scale = 1e-3;          ulabel = 's';
        otherwise
            error('make_time_axis:unit', ...
                'Unknown config.visualization.x_axis.unit "%s" (use frames | ms | s).', unit);
    end

    tax = struct();
    tax.unit                 = ulabel;
    tax.scale                = scale;
    tax.frame_dur_ms         = frame_dur;
    tax.from_ms              = @(ms) ms .* scale;                 % ms -> display
    tax.from_frames_rel      = @(nf) nf .* frame_dur .* scale;    % frames-from-ref -> display
    tax.range_frames         = rng;
    tax.montage_range_frames = mrng;
end
