function tf = hfs_pulses_resolvable(config)
    % HFS_PULSES_RESOLVABLE - Are single pulses of the train resolvable?
    %
    % True when the inter-pulse interval spans at least
    % config.hfs.min_frames_per_pulse frames at this recording's frame interval
    % (20 Hz: 10 frames at 5 ms -> yes; 2 frames at 25 ms -> no). When false,
    % only the train-level metrics (charge, plateau, decay) are computed.
    tf = (config.hfs.isi_ms / config.frame_dur_ms) >= config.hfs.min_frames_per_pulse;
end
