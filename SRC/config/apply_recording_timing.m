function config = apply_recording_timing(config, frame_interval_ms, n_frames)
    % APPLY_RECORDING_TIMING - Write a recording's metadata timing into the config
    %
    % Recordings.csv (not tracenorm_config) is the source of truth for the frame
    % interval and frame count. This sets every config field that is derived from
    % them, so code that reads config.frame_rate / frame_dur_ms / exposure_ms /
    % expected_frames sees the recording's real values and nothing falls back to
    % a hardcoded rate.
    %
    % Inputs:
    %   config            - tracenorm_config() struct
    %   frame_interval_ms - ms per frame (Recordings.csv Frame_Interval_ms)
    %   n_frames          - frames in the analysed trace (N_Frames, or the epoch
    %                       length for per-trial 1AP epochs)

    if ~(isscalar(frame_interval_ms) && isfinite(frame_interval_ms) && frame_interval_ms > 0)
        error('apply_recording_timing:interval', ...
            'Frame interval must be a positive finite number of ms (got %g).', frame_interval_ms);
    end
    if ~(isscalar(n_frames) && isfinite(n_frames) && n_frames >= 1)
        error('apply_recording_timing:frames', ...
            'Frame count must be a positive finite number (got %g).', n_frames);
    end

    config.exposure_ms          = frame_interval_ms;
    config.frame_dur_ms         = frame_interval_ms;
    config.frame_rate           = 1000 / frame_interval_ms;
    config.expected_frames      = n_frames;
    config.recording_duration_s = n_frames * frame_interval_ms / 1000;
    config.rolling_window_frames = round(config.rolling_window_sec * config.frame_rate);

    config.temporal.frame_duration_ms    = frame_interval_ms;
    config.temporal.recording_duration_s = config.recording_duration_s;
end
