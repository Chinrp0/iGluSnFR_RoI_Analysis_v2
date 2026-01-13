%% TEST_CORRECTED_BATCH - Test the corrected batch analysis with detailed diagnostics
clear; close all; clc;

fprintf('=== Testing Corrected Batch Analysis ===\n');

% Your data folder
data_folder = 'D:\Data\GluSnFR\Ms\2025-06-17_Ms-Hipp_DIV13_Doc2b_pilot_resave\iglusnfr4f_NGR\Spont\GPU_SNR_Processed\5_raw_mean';

% Options with extra debugging
options = struct();
options.verbose = true;
options.createPlots = true;
options.savePlots = false;
options.useParallel = false;
options.continueOnError = true;

% Run corrected batch analysis
try
    fprintf('Running batch analysis with improved error handling...\n');
    tic;
    results = batch_condition_analysis(data_folder, options);
    batch_time = toc;
    
    fprintf('\n=== CORRECTED BATCH RESULTS ===\n');
    fprintf('Analysis time: %.3f s\n', batch_time);
    
    % Check results
    if isfield(results, 'wt_data') && isfield(results, 'mut_data')
        wt = results.wt_data;
        mut = results.mut_data;
        
        fprintf('\nFile Processing Results:\n');
        fprintf('  WT: %d files processed successfully\n', wt.num_files);
        fprintf('  R213W: %d files processed successfully\n', mut.num_files);
        
        if mut.num_files > 0
            fprintf('\n🎉 SUCCESS: R213W files are now processing!\n');
            
            fprintf('\nCondition Comparison:\n');
            fprintf('  WT: %d ROIs, %d events (%.1f%% active)\n', ...
                wt.total_rois, wt.total_events, 100 * wt.stats.frequency.active_rois.fraction);
            fprintf('  R213W: %d ROIs, %d events (%.1f%% active)\n', ...
                mut.total_rois, mut.total_events, 100 * mut.stats.frequency.active_rois.fraction);
            
            % Show key statistical comparisons
            if isfield(results, 'comparison')
                comp = results.comparison;
                fprintf('\nStatistical Comparisons:\n');
                
                if isfield(comp, 'frequency_pvalue') && ~isnan(comp.frequency_pvalue)
                    fprintf('  Event Frequency: p = %.4f', comp.frequency_pvalue);
                    if comp.frequency_pvalue < 0.05
                        fprintf(' *SIGNIFICANT*');
                    end
                    fprintf('\n');
                end
                
                if isfield(comp, 'amplitude_pvalue') && ~isnan(comp.amplitude_pvalue)
                    fprintf('  Event Amplitude: p = %.4f', comp.amplitude_pvalue);
                    if comp.amplitude_pvalue < 0.05
                        fprintf(' *SIGNIFICANT*');
                    end
                    fprintf('\n');
                end
            end
            
        else
            fprintf('\n❌ STILL FAILING: R213W files are not processing in batch mode\n');
            
            % Show what files were attempted
            fprintf('\nFile Processing Details:\n');
            if isfield(results, 'processing_info')
                fprintf('  WT files found: %d\n', length(results.processing_info.wt_files));
                fprintf('  R213W files found: %d\n', length(results.processing_info.mut_files));
                
                fprintf('\nR213W files that should have been processed:\n');
                for i = 1:length(results.processing_info.mut_files)
                    [~, name, ~] = fileparts(results.processing_info.mut_files{i});
                    fprintf('    %d. %s\n', i, name);
                end
            end
        end
        
    else
        fprintf('❌ ERROR: Invalid results structure returned\n');
        fprintf('Available fields: %s\n', strjoin(fieldnames(results), ', '));
    end
    
catch ME
    fprintf('\n❌ BATCH ANALYSIS FAILED\n');
    fprintf('Error: %s\n', ME.message);
    if ~isempty(ME.stack)
        fprintf('In: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
    end
end

%% === Alternative: Manual File Processing Test ===

fprintf('\n=== Manual File Processing Test ===\n');
fprintf('Testing each file type individually...\n');

try
    % Get file lists
    csv_files = dir(fullfile(data_folder, '*.csv'));
    
    wt_files = {};
    mut_files = {};
    
    for i = 1:length(csv_files)
        file_path = fullfile(csv_files(i).folder, csv_files(i).name);
        filename = csv_files(i).name;
        
        if contains(filename, 'WT', 'IgnoreCase', true)
            wt_files{end+1} = file_path;
        elseif contains(filename, 'R213W', 'IgnoreCase', true)
            mut_files{end+1} = file_path;
        end
    end
    
    fprintf('Found: %d WT files, %d R213W files\n', length(wt_files), length(mut_files));
    
    % Test one file from each condition
    if ~isempty(wt_files)
        fprintf('\nTesting WT file...\n');
        [~, wt_name, ~] = fileparts(wt_files{1});
        fprintf('  File: %s\n', wt_name);
        
        wt_options = struct('verbose', false, 'createPlots', false);
        wt_result = main_pipeline(wt_files{1}, wt_options);
        
        if isfield(wt_result, 'fileResults') && ~isempty(wt_result.fileResults)
            wt_file_result = wt_result.fileResults{1};
            fprintf('  ✓ SUCCESS: %d ROIs, %d events\n', ...
                wt_file_result.numROIs, wt_file_result.event_stats.total_events);
        else
            fprintf('  ✗ FAILED: No file results\n');
        end
    end
    
    if ~isempty(mut_files)
        fprintf('\nTesting R213W file...\n');
        [~, mut_name, ~] = fileparts(mut_files{1});
        fprintf('  File: %s\n', mut_name);
        
        mut_options = struct('verbose', false, 'createPlots', false);
        mut_result = main_pipeline(mut_files{1}, mut_options);
        
        if isfield(mut_result, 'fileResults') && ~isempty(mut_result.fileResults)
            mut_file_result = mut_result.fileResults{1};
            fprintf('  ✓ SUCCESS: %d ROIs, %d events\n', ...
                mut_file_result.numROIs, mut_file_result.event_stats.total_events);
        else
            fprintf('  ✗ FAILED: No file results\n');
        end
    end
    
catch ME2
    fprintf('Manual test failed: %s\n', ME2.message);
end

fprintf('\n=== Test Complete ===\n');