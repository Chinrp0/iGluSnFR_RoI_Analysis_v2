function dt = experiment_to_datatype(experiment)
    % EXPERIMENT_TO_DATATYPE - Pipeline paradigm for a Recordings.csv Experiment
    %
    %   '1AP' (incl. multi-trial 5AP files)    -> '1AP'
    %   '2AP' / 'PPF'                          -> '2AP'
    %   'spontaneous' / 'spont'                -> 'Spont'
    % Anything else (e.g. HFS trains) has no analysis path yet and returns ''.

    e = lower(strtrim(char(experiment)));
    if any(strcmp(e, {'1ap'}))
        dt = '1AP';
    elseif any(strcmp(e, {'2ap', 'ppf'}))
        dt = '2AP';
    elseif startsWith(e, 'spont')
        dt = 'Spont';
    else
        dt = '';
    end
end
