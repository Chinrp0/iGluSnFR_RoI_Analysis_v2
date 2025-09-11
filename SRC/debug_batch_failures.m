%% DEBUG_BATCH_FAILURES - Diagnose why R213W files failed processing
% This script helps identify why certain files failed during batch processing

clear; close all; clc;

%% === CONFIGURATION ===
data_folder = 'D:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglusnfr4f_NGR\Spont\GPU_SNR_Processed\5_raw_mean';

fprintf('=== Debugging Batch Processing Failures ===\n');
fprintf('Data folder: %s\n', data_folder);

%% === STEP 1: Re-classify files and test individual processing ===

% Get file classification
[wt_files, mut_files, all_files] = classify_files_by_condition(data_folder);

fprintf('\nFile Classification:\n');
fprintf('WT files (%d):\n', length(wt_files));
for i = 1:length(wt_files)
    [~, name, ~] = fileparts(wt_files{i});
    fprintf('  %d. %s\n', i, name);
end

fprintf('R213W files (%d):\n', length(mut_files));
for i = 1:length(mut_files)
    [~, name, ~] = fileparts(mut_files{i});
    fprintf('  %d. %s\n', i, name);
end

%% === STEP 2: Test individual file processing ===

% Test first WT file
if ~isempty(wt_files)
    fprintf('\n=== Testing WT File Processing ===\n');
    wt_file = wt_files{1};
    [~, wt_name, ~] = fileparts(wt_file);
    fprintf('Testing: %s\n', wt_name);
    
    try
        % Test with verbose output
        options = struct();
        options.verbose = true;
        options.createPlots = false;
        options.continueOnError = false;
        
        tic;
        wt_result = main_pipeline(wt_file, options);
        wt_time = toc;
        
        fprintf('SUCCESS: WT file processed in %.3f s\n', wt_time);
        
        % Check result structure
        fprintf('WT result structure:\n');
        wt_fields = fieldnames(wt_result);
        for i = 1:length(wt_fields)
            field_name = wt_fields{i};
            field_value = wt_result.(field_name);
            if isstruct(field_value)
                fprintf('  %s: [struct with %d fields]\n', field_name, length(fieldnames(field_value)));
            elseif iscell(field_value)
                fprintf('  %s: [cell array %s]\n', field_name, mat2str(size(field_value)));
            else
                fprintf('  %s: [%s %s]\n', field_name, class(field_value), mat2str(size(field_value)));
            end
        end
        
        % Check fileResults specifically
        if isfield(wt_result, 'fileResults')
            fprintf('\nfileResults analysis:\n');
            if ~isempty(wt_result.fileResults)
                fprintf('  Length: %d\n', length(wt_result.fileResults));
                if ~isempty(wt_result.fileResults{1})
                    file_res_fields = fieldnames(wt_result.fileResults{1});
                    fprintf('  Fields: %s\n', strjoin(file_res_fields, ', '));
                    
                    % Check event_stats
                    if isfield(wt_result.fileResults{1}, 'event_stats')
                        event_fields = fieldnames(wt_result.fileResults{1}.event_stats);
                        fprintf('  event_stats fields: %s\n', strjoin(event_fields, ', '));
                    end
                else
                    fprintf('  fileResults{1} is EMPTY\n');
                end
            else
                fprintf('  fileResults is EMPTY\n');
            end
        else
            fprintf('  NO fileResults field!\n');
        end
        
    catch ME
        fprintf('FAILED: WT file processing failed\n');
        fprintf('Error: %s\n', ME.message);
        if ~isempty(ME.stack)
            fprintf('  In: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
        end
    end
end

%% === STEP 3: Test R213W file processing ===

if ~isempty(mut_files)
    fprintf('\n=== Testing R213W File Processing ===\n');
    mut_file = mut_files{1};
    [~, mut_name, ~] = fileparts(mut_file);
    fprintf('Testing: %s\n', mut_name);
    
    try
        % Test with verbose output
        options = struct();
        options.verbose = true;
        options.createPlots = false;
        options.continueOnError = false;
        
        tic;
        mut_result = main_pipeline(mut_file, options);
        mut_time = toc;
        
        fprintf('SUCCESS: R213W file processed in %.3f s\n', mut_time);
        
        % Check result structure
        fprintf('R213W result structure:\n');
        mut_fields = fieldnames(mut_result);
        for i = 1:length(mut_fields)
            field_name = mut_fields{i};
            field_value = mut_result.(field_name);
            if isstruct(field_value)
                fprintf('  %s: [struct with %d fields]\n', field_name, length(fieldnames(field_value)));
            elseif iscell(field_value)
                fprintf('  %s: [cell array %s]\n', field_name, mat2str(size(field_value)));
            else
                fprintf('  %s: [%s %s]\n', field_name, class(field_value), mat2str(size(field_value)));
            end
        end
        
        % Check fileResults specifically
        if isfield(mut_result, 'fileResults')
            fprintf('\nfileResults analysis:\n');
            if ~isempty(mut_result.fileResults)
                fprintf('  Length: %d\n', length(mut_result.fileResults));
                if ~isempty(mut_result.fileResults{1})
                    file_res_fields = fieldnames(mut_result.fileResults{1});
                    fprintf('  Fields: %s\n', strjoin(file_res_fields, ', '));
                    
                    % Check event_stats
                    if isfield(mut_result.fileResults{1}, 'event_stats')
                        event_fields = fieldnames(mut_result.fileResults{1}.event_stats);
                        fprintf('  event_stats fields: %s\n', strjoin(event_fields, ', '));
                    end
                else
                    fprintf('  fileResults{1} is EMPTY\n');
                end
            else
                fprintf('  fileResults is EMPTY\n');
            end
        else
            fprintf('  NO fileResults field!\n');
        end
        
    catch ME
        fprintf('FAILED: R213W file processing failed\n');
        fprintf('Error: %s\n', ME.message);
        if ~isempty(ME.stack)
            fprintf('  In: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
            
            % Print more detailed stack trace
            fprintf('  Full stack trace:\n');
            for j = 1:min(3, length(ME.stack))
                fprintf('    %d. %s (line %d)\n', j, ME.stack(j).name, ME.stack(j).line);
            end
        end
        
        % Try to identify the specific issue
        fprintf('\n=== Diagnosing R213W File Issue ===\n');
        
        % Check if file exists and is readable
        if exist(mut_file, 'file')
            fprintf('File exists: YES\n');
            file_info = dir(mut_file);
            fprintf('File size: %.2f MB\n', file_info.bytes / (1024^2));
        else
            fprintf('File exists: NO - This is the problem!\n');
        end
        
        % Try loading just the CSV data
        try
            fprintf('Testing CSV loading...\n');
            data = readmatrix(mut_file);
            fprintf('CSV load: SUCCESS - %d x %d matrix\n', size(data, 1), size(data, 2));
            
            % Check for obvious data issues
            if any(isnan(data(:)))
                nan_fraction = sum(isnan(data(:))) / numel(data);
                fprintf('Data contains %.2f%% NaN values\n', 100 * nan_fraction);
            else
                fprintf('Data contains no NaN values\n');
            end
            
            if any(isinf(data(:)))
                inf_count = sum(isinf(data(:)));
                fprintf('Data contains %d infinite values\n', inf_count);
            else
                fprintf('Data contains no infinite values\n');
            end
            
            fprintf('Data range: %.6f to %.6f\n', min(data(:)), max(data(:)));
            
        catch CSV_ME
            fprintf('CSV load: FAILED - %s\n', CSV_ME.message);
        end
    end
end

%% === STEP 4: Test batch processing with more detailed error reporting ===

fprintf('\n=== Testing Batch Processing with Detailed Errors ===\n');

% Modified version with better error reporting
all_test_files = [wt_files; mut_files];
file_results = cell(length(all_test_files), 1);

for i = 1:length(all_test_files)
    file_path = all_test_files{i};
    [~, filename, ~] = fileparts(file_path);
    
    fprintf('\nProcessing %d/%d: %s\n', i, length(all_test_files), filename);
    
    try
        % Use existing pipeline with minimal options
        pipeline_options = struct();
        pipeline_options.createPlots = false;
        pipeline_options.verbose = false;
        pipeline_options.continueOnError = false;  % Don't continue on error for debugging
        
        result = main_pipeline(file_path, pipeline_options);
        
        % Check the result structure
        if isfield(result, 'fileResults') && ~isempty(result.fileResults)
            file_results{i} = result.fileResults{1};
            fprintf('  SUCCESS: %d ROIs, %d events\n', ...
                file_results{i}.numROIs, file_results{i}.event_stats.total_events);
        else
            file_results{i} = [];
            fprintf('  FAILED: No fileResults returned\n');
        end
        
    catch ME
        file_results{i} = [];
        fprintf('  FAILED: %s\n', ME.message);
        
        % More detailed error analysis
        if contains(ME.message, 'baseline') || contains(ME.message, 'Baseline')
            fprintf('    → Baseline calculation issue\n');
        elseif contains(ME.message, 'event') || contains(ME.message, 'Event')
            fprintf('    → Event detection issue\n');
        elseif contains(ME.message, 'dF/F') || contains(ME.message, 'dfof')
            fprintf('    → dF/F calculation issue\n');
        elseif contains(ME.message, 'size') || contains(ME.message, 'dimension')
            fprintf('    → Data size/dimension issue\n');
        elseif contains(ME.message, 'NaN') || contains(ME.message, 'nan')
            fprintf('    → NaN data issue\n');
        else
            fprintf('    → Unknown issue type\n');
        end
    end
end

%% === STEP 5: Summary ===

fprintf('\n=== SUMMARY ===\n');

successful_files = sum(~cellfun(@isempty, file_results));
fprintf('Successfully processed: %d/%d files\n', successful_files, length(all_test_files));

% Separate by condition
wt_success = sum(~cellfun(@isempty, file_results(1:length(wt_files))));
mut_success = sum(~cellfun(@isempty, file_results(length(wt_files)+1:end)));

fprintf('WT success rate: %d/%d (%.1f%%)\n', wt_success, length(wt_files), 100*wt_success/length(wt_files));
fprintf('R213W success rate: %d/%d (%.1f%%)\n', mut_success, length(mut_files), 100*mut_success/length(mut_files));

if mut_success == 0
    fprintf('\n*** R213W FILES ARE FAILING ***\n');
    fprintf('Common solutions:\n');
    fprintf('1. Check file paths and permissions\n');
    fprintf('2. Verify CSV file format and contents\n');
    fprintf('3. Check for data quality issues (NaN, Inf, wrong dimensions)\n');
    fprintf('4. Run individual R213W file with verbose=true to see detailed error\n');
else
    fprintf('\nBoth conditions processed successfully - original batch issue was different\n');
end

%% === Helper function (inline) ===
function [wt_files, mut_files, all_files] = classify_files_by_condition(folder_path)
    % Inline version of classification function for debugging
    
    csv_files = dir(fullfile(folder_path, '*.csv'));
    
    if isempty(csv_files)
        error('No CSV files found in folder: %s', folder_path);
    end
    
    all_files = {};
    wt_files = {};
    mut_files = {};
    
    for i = 1:length(csv_files)
        file_path = fullfile(csv_files(i).folder, csv_files(i).name);
        filename = csv_files(i).name;
        all_files{end+1} = file_path;
        
        % Classify based on filename patterns
        if contains(filename, 'Doc2b-WT', 'IgnoreCase', true) || ...
           contains(filename, 'WT', 'IgnoreCase', true)
            wt_files{end+1} = file_path;
            
        elseif contains(filename, 'Doc2b-R213W', 'IgnoreCase', true) || ...
               contains(filename, 'R213W', 'IgnoreCase', true)
            mut_files{end+1} = file_path;
        end
    end
end