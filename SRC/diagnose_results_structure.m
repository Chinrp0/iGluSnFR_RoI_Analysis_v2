function diagnose_results_structure(results)
    % DIAGNOSE_RESULTS_STRUCTURE - Quality control script to inspect data organization
    % 
    % Usage:
    %   diagnose_results_structure(results)
    %
    % This will help us understand:
    %   1. Why file counting shows 1 instead of actual counts (10 WT, 7 R213W)
    %   2. The exact structure of file_results for proper extraction
    %   3. Data types and dimensions for fixing IEI concatenation errors
    
    fprintf('\n=== RESULTS STRUCTURE DIAGNOSTIC ===\n\n');
    
    %% 1. TOP-LEVEL STRUCTURE
    fprintf('1. TOP-LEVEL FIELDS:\n');
    fprintf('-------------------\n');
    top_fields = fieldnames(results);
    for i = 1:length(top_fields)
        fprintf('  • %s\n', top_fields{i});
    end
    
    %% 2. FILE RESULTS STRUCTURE
    fprintf('\n2. FILE RESULTS STRUCTURE:\n');
    fprintf('-------------------------\n');
    
    if isfield(results, 'file_results')
        fr = results.file_results;
        fr_fields = fieldnames(fr);
        fprintf('  Fields in file_results: %s\n', strjoin(fr_fields, ', '));
        
        % Check WT results
        if isfield(fr, 'wt_results')
            fprintf('\n  WT RESULTS:\n');
            fprintf('    Type: %s\n', class(fr.wt_results));
            fprintf('    Size: %s\n', mat2str(size(fr.wt_results)));
            fprintf('    Length: %d\n', length(fr.wt_results));
            
            if iscell(fr.wt_results)
                fprintf('    → Cell array with %d elements\n', numel(fr.wt_results));
                
                % Count non-empty cells
                non_empty = 0;
                for i = 1:length(fr.wt_results)
                    if ~isempty(fr.wt_results{i})
                        non_empty = non_empty + 1;
                    end
                end
                fprintf('    → Non-empty cells: %d\n', non_empty);
                
            elseif isstruct(fr.wt_results)
                fprintf('    → Struct array with %d elements\n', numel(fr.wt_results));
            end
        end
        
        % Check R213W results
        if isfield(fr, 'mut_results')
            fprintf('\n  R213W RESULTS:\n');
            fprintf('    Type: %s\n', class(fr.mut_results));
            fprintf('    Size: %s\n', mat2str(size(fr.mut_results)));
            fprintf('    Length: %d\n', length(fr.mut_results));
            
            if iscell(fr.mut_results)
                fprintf('    → Cell array with %d elements\n', numel(fr.mut_results));
                
                % Count non-empty cells
                non_empty = 0;
                for i = 1:length(fr.mut_results)
                    if ~isempty(fr.mut_results{i})
                        non_empty = non_empty + 1;
                    end
                end
                fprintf('    → Non-empty cells: %d\n', non_empty);
                
            elseif isstruct(fr.mut_results)
                fprintf('    → Struct array with %d elements\n', numel(fr.mut_results));
            end
        end
    else
        fprintf('  ⚠ file_results field NOT FOUND!\n');
        return;
    end
    
    %% 3. PROCESSING INFO STRUCTURE
    fprintf('\n3. PROCESSING INFO:\n');
    fprintf('------------------\n');
    
    if isfield(results, 'processing_info')
        pi = results.processing_info;
        pi_fields = fieldnames(pi);
        
        for i = 1:length(pi_fields)
            field = pi_fields{i};
            value = pi.(field);
            
            if ischar(value) || isstring(value)
                fprintf('  %s: %s\n', field, value);
            elseif isnumeric(value) && isscalar(value)
                fprintf('  %s: %d\n', field, value);
            elseif iscell(value)
                fprintf('  %s: cell array [%s] with %d elements\n', ...
                    field, mat2str(size(value)), numel(value));
            else
                fprintf('  %s: %s\n', field, class(value));
            end
        end
        
        % Show file paths if available
        if isfield(pi, 'wt_files')
            fprintf('\n  WT FILES (%d total):\n', length(pi.wt_files));
            for i = 1:min(3, length(pi.wt_files))
                [~, fname, ext] = fileparts(pi.wt_files{i});
                fprintf('    [%d] %s%s\n', i, fname, ext);
            end
            if length(pi.wt_files) > 3
                fprintf('    ... and %d more\n', length(pi.wt_files) - 3);
            end
        end
        
        if isfield(pi, 'mut_files')
            fprintf('\n  R213W FILES (%d total):\n', length(pi.mut_files));
            for i = 1:min(3, length(pi.mut_files))
                [~, fname, ext] = fileparts(pi.mut_files{i});
                fprintf('    [%d] %s%s\n', i, fname, ext);
            end
            if length(pi.mut_files) > 3
                fprintf('    ... and %d more\n', length(pi.mut_files) - 3);
            end
        end
    else
        fprintf('  ⚠ processing_info field NOT FOUND!\n');
    end
    
    %% 4. SAMPLE FILE DATA STRUCTURE
    fprintf('\n4. SAMPLE FILE DATA STRUCTURE:\n');
    fprintf('-----------------------------\n');
    
    % Get first non-empty WT file
    sample_file = [];
    if isfield(results.file_results, 'wt_results')
        for i = 1:length(results.file_results.wt_results)
            if ~isempty(results.file_results.wt_results{i})
                sample_file = results.file_results.wt_results{i};
                fprintf('  Examining WT file #%d:\n', i);
                break;
            end
        end
    end
    
    if ~isempty(sample_file)
        sample_fields = fieldnames(sample_file);
        fprintf('  Fields in individual file result:\n');
        for i = 1:length(sample_fields)
            field = sample_fields{i};
            value = sample_file.(field);
            
            if isnumeric(value) && isscalar(value)
                fprintf('    • %s: %d\n', field, value);
            elseif isnumeric(value) && ~isscalar(value)
                fprintf('    • %s: [%s] %s\n', field, mat2str(size(value)), class(value));
            elseif islogical(value)
                fprintf('    • %s: [%s] logical\n', field, mat2str(size(value)));
            elseif isstruct(value)
                fprintf('    • %s: struct with %d fields\n', field, length(fieldnames(value)));
            else
                fprintf('    • %s: %s\n', field, class(value));
            end
        end
        
        % Check event_mask specifically (critical for IEI analysis)
        if isfield(sample_file, 'event_mask')
            fprintf('\n  EVENT_MASK DETAILS:\n');
            em = sample_file.event_mask;
            fprintf('    Size: [%d frames × %d ROIs]\n', size(em, 1), size(em, 2));
            fprintf('    Type: %s\n', class(em));
            fprintf('    True values: %d (%.2f%%)\n', sum(em(:)), 100*sum(em(:))/numel(em));
        end
        
        % Check event_stats specifically
        if isfield(sample_file, 'event_stats')
            fprintf('\n  EVENT_STATS DETAILS:\n');
            es = sample_file.event_stats;
            es_fields = fieldnames(es);
            for i = 1:length(es_fields)
                field = es_fields{i};
                value = es.(field);
                
                if isnumeric(value) && isscalar(value)
                    fprintf('    • %s: %d\n', field, value);
                elseif isnumeric(value) && isvector(value)
                    fprintf('    • %s: vector [%d elements]\n', field, length(value));
                elseif isnumeric(value) && ~isscalar(value)
                    fprintf('    • %s: [%s] %s\n', field, mat2str(size(value)), class(value));
                else
                    fprintf('    • %s: %s\n', field, class(value));
                end
            end
        end
    else
        fprintf('  ⚠ No valid file data found in wt_results!\n');
    end
    
    %% 5. IEI DATA DIMENSIONS TEST
    fprintf('\n5. IEI EXTRACTION TEST:\n');
    fprintf('----------------------\n');
    
    if isfield(results, 'file_results')
        fprintf('  Testing extract_all_ieis() function...\n\n');
        
        % Test WT
        if isfield(results.file_results, 'wt_results')
            [wt_ieis, wt_stats] = extract_all_ieis(results.file_results.wt_results, 100);
            fprintf('  WT IEIs:\n');
            fprintf('    Total IEIs: %d\n', length(wt_ieis));
            fprintf('    IEI array size: %s\n', mat2str(size(wt_ieis)));
            fprintf('    IEI array type: %s\n', class(wt_ieis));
            fprintf('    Is column vector: %s\n', mat2str(iscolumn(wt_ieis)));
            fprintf('    Is row vector: %s\n', mat2str(isrow(wt_ieis)));
            fprintf('    ROIs with ≥2 events: %d\n', wt_stats.num_rois_with_multiple_events);
            fprintf('    CV array size: %s\n', mat2str(size(wt_stats.cv_per_roi)));
            
            if ~isempty(wt_ieis)
                fprintf('    Sample IEIs (first 5): %s\n', mat2str(wt_ieis(1:min(5,end))'));
            end
        end
        
        fprintf('\n');
        
        % Test R213W
        if isfield(results.file_results, 'mut_results')
            [mut_ieis, mut_stats] = extract_all_ieis(results.file_results.mut_results, 100);
            fprintf('  R213W IEIs:\n');
            fprintf('    Total IEIs: %d\n', length(mut_ieis));
            fprintf('    IEI array size: %s\n', mat2str(size(mut_ieis)));
            fprintf('    IEI array type: %s\n', class(mut_ieis));
            fprintf('    Is column vector: %s\n', mat2str(iscolumn(mut_ieis)));
            fprintf('    Is row vector: %s\n', mat2str(isrow(mut_ieis)));
            fprintf('    ROIs with ≥2 events: %d\n', mut_stats.num_rois_with_multiple_events);
            fprintf('    CV array size: %s\n', mat2str(size(mut_stats.cv_per_roi)));
            
            if ~isempty(mut_ieis)
                fprintf('    Sample IEIs (first 5): %s\n', mat2str(mut_ieis(1:min(5,end))'));
            end
        end
    end
    
    %% 6. BIOLOGICAL VARIABILITY METRICS TEST
    fprintf('\n6. BIOLOGICAL VARIABILITY EXTRACTION TEST:\n');
    fprintf('------------------------------------------\n');
    
    if isfield(results, 'file_results') && isfield(results, 'processing_info')
        fprintf('  Testing extract_file_metrics() function...\n\n');
        
        % Test WT
        if isfield(results.file_results, 'wt_results') && isfield(results.processing_info, 'wt_files')
            [wt_metrics, wt_filenames] = extract_file_metrics(results.file_results.wt_results, ...
                results.processing_info.wt_files);
            
            fprintf('  WT Metrics:\n');
            fprintf('    Files processed: %d\n', length(wt_metrics.mean_frequency));
            fprintf('    Non-zero frequencies: %d\n', sum(wt_metrics.mean_frequency > 0));
            fprintf('    Mean frequency array size: %s\n', mat2str(size(wt_metrics.mean_frequency)));
            fprintf('    Total ROIs array size: %s\n', mat2str(size(wt_metrics.total_rois)));
            
            if ~isempty(wt_filenames)
                fprintf('    Sample filenames:\n');
                for i = 1:min(3, length(wt_filenames))
                    fprintf('      [%d] %s (freq=%.4f Hz, %d events)\n', ...
                        i, wt_filenames{i}, wt_metrics.mean_frequency(i), wt_metrics.total_events(i));
                end
            end
        end
        
        fprintf('\n');
        
        % Test R213W
        if isfield(results.file_results, 'mut_results') && isfield(results.processing_info, 'mut_files')
            [mut_metrics, mut_filenames] = extract_file_metrics(results.file_results.mut_results, ...
                results.processing_info.mut_files);
            
            fprintf('  R213W Metrics:\n');
            fprintf('    Files processed: %d\n', length(mut_metrics.mean_frequency));
            fprintf('    Non-zero frequencies: %d\n', sum(mut_metrics.mean_frequency > 0));
            fprintf('    Mean frequency array size: %s\n', mat2str(size(mut_metrics.mean_frequency)));
            fprintf('    Total ROIs array size: %s\n', mat2str(size(mut_metrics.total_rois)));
            
            if ~isempty(mut_filenames)
                fprintf('    Sample filenames:\n');
                for i = 1:min(3, length(mut_filenames))
                    fprintf('      [%d] %s (freq=%.4f Hz, %d events)\n', ...
                        i, mut_filenames{i}, mut_metrics.mean_frequency(i), mut_metrics.total_events(i));
                end
            end
        end
    end
    
    %% 7. SUMMARY OF ISSUES
    fprintf('\n7. DIAGNOSTIC SUMMARY:\n');
    fprintf('---------------------\n');
    
    issues_found = false;
    
    % Check file counting issue
    if isfield(results, 'processing_info')
        expected_wt = length(results.processing_info.wt_files);
        expected_mut = length(results.processing_info.mut_files);
        
        if isfield(results, 'file_results')
            actual_wt = length(results.file_results.wt_results);
            actual_mut = length(results.file_results.mut_results);
            
            if expected_wt ~= actual_wt || expected_mut ~= actual_mut
                fprintf('  ⚠ FILE COUNT MISMATCH:\n');
                fprintf('    Expected: %d WT, %d R213W\n', expected_wt, expected_mut);
                fprintf('    Actual: %d WT, %d R213W\n', actual_wt, actual_mut);
                issues_found = true;
            end
        end
    end
    
    % Check IEI dimension issue
    if exist('wt_ieis', 'var') && exist('mut_ieis', 'var')
        if ~isempty(wt_ieis) && ~isempty(mut_ieis)
            if size(wt_ieis, 2) ~= size(mut_ieis, 2) || size(wt_ieis, 1) ~= size(mut_ieis, 1)
                fprintf('  ⚠ IEI DIMENSION MISMATCH:\n');
                fprintf('    WT IEIs: %s\n', mat2str(size(wt_ieis)));
                fprintf('    R213W IEIs: %s\n', mat2str(size(mut_ieis)));
                fprintf('    → This will cause concatenation errors in boxplot\n');
                issues_found = true;
            end
        end
    end
    
    if ~issues_found
        fprintf('  ✓ No obvious structural issues detected\n');
    end
    
    fprintf('\n=== END OF DIAGNOSTIC ===\n\n');
end

function [all_ieis, stats] = extract_all_ieis(file_results, frame_rate)
    % Extract all inter-event intervals from file results
    % (Copied from plot_iei_comparison.m for testing)
    
    all_ieis = [];
    cv_per_roi = [];
    num_rois_with_multiple_events = 0;
    
    for file_idx = 1:length(file_results)
        if isempty(file_results{file_idx})
            continue;
        end
        
        file_result = file_results{file_idx};
        
        if ~isfield(file_result, 'event_mask') || isempty(file_result.event_mask)
            continue;
        end
        
        event_mask = file_result.event_mask;
        [num_frames, num_rois] = size(event_mask);
        
        % Process each ROI
        for roi = 1:num_rois
            roi_events = event_mask(:, roi);
            
            % Find event start times (transitions from false to true)
            event_starts = find(diff([false; roi_events]) == 1);
            
            if length(event_starts) >= 2
                % Calculate IEIs in seconds
                ieis_frames = diff(event_starts);
                ieis_seconds = ieis_frames / frame_rate;
                
                all_ieis = [all_ieis; ieis_seconds];
                
                % Calculate CV for this ROI
                if length(ieis_seconds) >= 2
                    roi_cv = std(ieis_seconds) / mean(ieis_seconds);
                    cv_per_roi = [cv_per_roi; roi_cv];
                    num_rois_with_multiple_events = num_rois_with_multiple_events + 1;
                end
            end
        end
    end
    
    stats = struct();
    stats.num_rois_with_multiple_events = num_rois_with_multiple_events;
    stats.cv_per_roi = cv_per_roi;
end

function [metrics, filenames] = extract_file_metrics(file_results, file_paths)
    % Extract per-file summary metrics
    % (Copied from plot_biological_variability.m for testing)
    
    n_files = length(file_results);
    
    metrics = struct();
    metrics.mean_frequency = zeros(n_files, 1);
    metrics.mean_amplitude = zeros(n_files, 1);
    metrics.active_roi_fraction = zeros(n_files, 1);
    metrics.total_events = zeros(n_files, 1);
    metrics.active_rois = zeros(n_files, 1);
    metrics.total_rois = zeros(n_files, 1);
    
    filenames = cell(n_files, 1);
    
    for i = 1:n_files
        if isempty(file_results{i})
            continue;
        end
        
        file_result = file_results{i};
        
        % Extract filename
        [~, fname, ~] = fileparts(file_paths{i});
        filenames{i} = fname;
        
        % Event statistics
        if isfield(file_result, 'event_stats')
            es = file_result.event_stats;
            
            % Mean frequency (only active ROIs)
            if isfield(es, 'frequency_hz')
                active_freqs = es.frequency_hz(es.frequency_hz > 0);
                if ~isempty(active_freqs)
                    metrics.mean_frequency(i) = mean(active_freqs);
                end
            end
            
            % Mean amplitude
            if isfield(es, 'all_event_peaks') && ~isempty(es.all_event_peaks)
                metrics.mean_amplitude(i) = mean(es.all_event_peaks);
            end
            
            % Total events
            if isfield(es, 'total_events')
                metrics.total_events(i) = es.total_events;
            end
            
            % Active ROIs
            if isfield(es, 'rois_with_events')
                metrics.active_rois(i) = es.rois_with_events;
            end
        end
        
        % Total ROIs
        if isfield(file_result, 'numROIs')
            metrics.total_rois(i) = file_result.numROIs;
            metrics.active_roi_fraction(i) = metrics.active_rois(i) / metrics.total_rois(i);
        end
    end
end
