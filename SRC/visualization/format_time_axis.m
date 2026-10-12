function format_time_axis(ax, tax, refLabel, stim_frame)
    % FORMAT_TIME_AXIS - label a time x-axis in the configured unit and, when a
    % frame range is configured, clamp the axis to it.
    %
    %   ax         - target axes
    %   tax        - struct from make_time_axis()
    %   refLabel   - the quantity on the x-axis, e.g. 'Time from stimulus'. The
    %                unit is appended as ' (<unit>)'. Pass '' to skip the label.
    %   stim_frame - (optional) stimulus frame that x = 0 corresponds to. Pass
    %                it to apply tax.range_frames (absolute frames) as xlim;
    %                omit or pass [] to leave the current limits untouched.
    %
    % The plotted data must ALREADY be in display units (convert with
    % tax.from_ms / tax.from_frames_rel before plotting).

    if ~isempty(refLabel)
        xlabel(ax, sprintf('%s (%s)', refLabel, tax.unit));
    end

    if nargin >= 4 && ~isempty(stim_frame) && ~isempty(tax.range_frames)
        % Absolute frame range -> display units, relative to the stim frame.
        xlim(ax, sort(tax.from_frames_rel(tax.range_frames(:).' - stim_frame)));
    end
end
