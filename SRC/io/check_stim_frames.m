function [r, report] = check_stim_frames(raw, r, config)
    % CHECK_STIM_FRAMES - Choose refined vs TTL stim frames from the ROI response
    %
    % Recordings.csv gives two stim frames per pulse: the raw TTL frame and the
    % refined frame, which the macro derives from the FULL field of view. The
    % refinement normally corrects +/-1 frame of TTL jitter, but the full-FoV
    % signal can also place it on the wrong frame (e.g. +3 frames, after the ROI
    % response has already started). This check looks at the analysed ROIs
    % themselves and only overrides the refined frame when it is impossible:
    %
    %   onset = first of min_frames_above consecutive frames where the ROI-mean
    %           dF/F exceeds the level just before the search window by
    %           onset_sigma x its noise SD (dF/F per ROI from that same window)
    %   a candidate stim s is CONSISTENT when  s <= onset  and
    %           onset - s <= max_onset_latency_ms  (the response cannot start
    %           before the stimulus, and iGluSnFR responds within ~1-2 frames)
    %
    %   refined consistent                    -> refined (the default)
    %   refined inconsistent, TTL consistent  -> TTL (flag STIM_TTL_USED_T<k>)
    %   neither / no onset found              -> refined, noted in the report
    %
    % The decision is made per TRIAL on its first pulse and applied to all of
    % that trial's pulses (the macro shifts every pulse of a trial by the same
    % amount; a 2AP pair's second pulse rides on the first's decay, so it is not
    % a clean onset). Trials = N_Trials; pulses per trial = stims / N_Trials.
    %
    % Inputs:
    %   raw    - [frames x ROIs] raw fluorescence of the file
    %   r      - one read_recordings_metadata row
    %   config - tracenorm_config() (uses config.metadata.stim_check)
    %
    % Outputs:
    %   r      - row with .stim_frames set to the chosen frames, plus
    %            .stim_frames_refined (the original refined frames) and
    %            .stim_source ({1 x nTrials} 'refined' | 'ttl')
    %   report - '' when nothing to note, else a short text for the summary

    c = config.metadata.stim_check;
    r.stim_frames_refined = r.stim_frames;
    nS = numel(r.stim_frames);
    nT = r.n_trials;
    if ~(isfinite(nT) && nT >= 1 && mod(nS, nT) == 0), nT = 1; end
    per = nS / nT;
    r.stim_source = repmat({'refined'}, 1, nT);
    report = '';
    if ~c.enabled || isempty(r.stim_frames) || numel(r.ttl_stim_frames) ~= nS ...
            || isequal(r.ttl_stim_frames, r.stim_frames)
        return;
    end

    fi    = r.frame_interval_ms;
    nF    = size(raw, 1);
    nref  = max(3, round(c.ref_ms / fi));
    nback = max(2, round(c.search_back_ms / fi));
    nlat  = max(1, round(c.max_onset_latency_ms / fi));
    notes = {};
    for k = 1:nT
        i  = (k - 1) * per + 1;                  % first pulse of trial k
        rf = r.stim_frames(i); tt = r.ttl_stim_frames(i);
        if rf == tt, continue; end
        lo = min(rf, tt) - nback; hi = min(max(rf, tt) + nlat, nF);
        if lo - nref < 1
            notes{end+1} = sprintf('T%d: too close to the start to check', k); %#ok<AGROW>
            continue;
        end
        R  = double(raw(lo - nref:hi, :));
        F0 = mean(R(1:nref, :), 1, 'omitnan');
        ok = isfinite(F0) & F0 > 0;
        m  = mean((R(:, ok) - F0(ok)) ./ F0(ok), 2, 'omitnan');
        sd = mad(m(1:nref), 1) * 1.4826;
        if ~(sd > 0), sd = std(m(1:nref)); end
        above = m(nref + 1:end) > c.onset_sigma * sd;
        run = conv(double(above), ones(c.min_frames_above, 1), 'valid') == c.min_frames_above;
        o = find(run, 1);
        if isempty(o)
            notes{end+1} = sprintf('T%d: no ROI onset found, refined kept', k); %#ok<AGROW>
            continue;
        end
        onset = lo + o - 1;
        fits = @(s) s <= onset && onset - s <= nlat;
        if fits(rf), continue; end
        if fits(tt)
            idx = i:i + per - 1;
            r.stim_frames(idx) = r.ttl_stim_frames(idx);
            r.stim_source{k} = 'ttl';
            tag = sprintf('STIM_TTL_USED_T%d(onset %d, refined %d, TTL %d)', k, onset, rf, tt);
            r.flags{end+1} = tag;
            notes{end+1} = tag; %#ok<AGROW>
        else
            notes{end+1} = sprintf('T%d: ROI onset %d fits neither refined %d nor TTL %d, refined kept', ...
                k, onset, rf, tt); %#ok<AGROW>
        end
    end
    report = strjoin(notes, '; ');
end
