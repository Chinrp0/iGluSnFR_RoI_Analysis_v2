function config = apply_train_timing(config, stim1_frame, n_pulses, train_hz)
    % APPLY_TRAIN_TIMING - Write an HFS train's pulse frames into config.hfs
    %
    % Pulse n sits at stim1 + round((n-1) * isi_ms / frame_dur_ms), rounded from
    % the start of the train so the pulses do not drift (49 x 50 ms at 5.0086 ms
    % frames = 489.1 frames, not 49 x 10). Call after apply_recording_timing (it
    % needs config.frame_dur_ms).
    %
    % Inputs:
    %   config      - tracenorm_config() struct with the recording timing applied
    %   stim1_frame - first pulse, as an index into the analysed trace / epoch
    %   n_pulses    - Recordings.csv Train_Pulses
    %   train_hz    - Recordings.csv Train_Hz
    %
    % Sets config.hfs.stim1_frame, .n_pulses, .train_hz, .isi_ms, .pulse_frames
    % and .train_end_frame (one inter-pulse interval after the last pulse).

    if ~(isfinite(n_pulses) && n_pulses >= 1 && n_pulses == round(n_pulses))
        error('apply_train_timing:pulses', 'Train_Pulses must be a positive integer (got %g).', n_pulses);
    end
    if ~(isfinite(train_hz) && train_hz > 0)
        error('apply_train_timing:hz', 'Train_Hz must be a positive number (got %g).', train_hz);
    end
    fd     = config.frame_dur_ms;
    isi_ms = 1000 / train_hz;
    config.hfs.stim1_frame     = stim1_frame;
    config.hfs.n_pulses        = n_pulses;
    config.hfs.train_hz        = train_hz;
    config.hfs.isi_ms          = isi_ms;
    config.hfs.pulse_frames    = stim1_frame + round((0:n_pulses - 1) * isi_ms / fd);
    config.hfs.train_end_frame = stim1_frame + round(n_pulses * isi_ms / fd);
end
