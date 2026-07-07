function setup_pipeline()
    % SETUP_PIPELINE - Initialize the fluorescent imaging analysis environment
    % UPDATED: Now verifies corrected Schmitt trigger detector and comparison tools
    % Run this once per MATLAB session to configure paths and dependencies
    
    fprintf('=== Setting up Fluorescent Imaging Pipeline with Corrected Schmitt Trigger ===\n');
    
    % Get the current directory (should be SRC folder)
    srcPath = pwd;
    fprintf('Source directory: %s\n', srcPath);
    
    % Add all subfolders to MATLAB path recursively
    fprintf('Adding all subfolders to MATLAB path...\n');
    addpath(genpath(srcPath));
    
    % Verify key modules are accessible
    fprintf('Verifying module accessibility:\n');
    
    modules = {
        'csv_loader_v2', 'io/csv_loader_v2.m';
        'main_pipeline', 'main_pipeline.m';
        'tracenorm_config', 'config/tracenorm_config.m';
        'baseline_detector', 'tracenorm/baseline_detector.m';
        'dfof_calculator', 'tracenorm/dfof_calculator.m';
        'corrected_schmitt_detector', 'tracenorm/corrected_schmitt_detector.m';  % NEW: Corrected detector
        'pure_schmitt_trigger_detector', 'tracenorm/pure_schmitt_trigger_detector.m';  % Legacy
        'quality_assessor', 'tracenorm/quality_assessor.m';
        'baseline_plotter', 'visualization/baseline_plotter_clean.m';
        'prestim_baseline', 'tracenorm/prestim_baseline.m';            % 1AP: pre-stim baseline
        'oneap_analyzer', 'tracenorm/oneap_analyzer.m';                % 1AP: per-ROI metrics
        'fit_decay_double', 'tracenorm/fit_decay_double.m';            % 1AP: averaged double-exp decay
        'oneap_pipeline', 'oneap_pipeline.m';                          % 1AP: single-file path
        'oneap_batch_analysis', 'oneap_batch_analysis.m';             % 1AP: batch driver
        'twoap_analyzer', 'tracenorm/twoap_analyzer.m';                % 2AP: per-ROI PPR metrics
        'twoap_pipeline', 'twoap_pipeline.m';                          % 2AP: single-file path
        'twoap_batch_analysis', 'twoap_batch_analysis.m';             % 2AP: batch driver (PPR vs ISI)
        'parse_group', 'io/parse_group.m';                             % shared: config-driven group parse
        'parse_coverslip_id', 'io/parse_coverslip_id.m';               % shared: coverslip id
        'parse_isi_ms', 'io/parse_isi_ms.m';                           % shared: PPF ISI from filename
        'draw_box', 'visualization/draw_box.m';                        % shared: native box/whisker
        'make_time_axis', 'visualization/make_time_axis.m';            % shared: figure x-axis unit
        'format_time_axis', 'visualization/format_time_axis.m';        % shared: x-axis label + range
        'montage_window', 'visualization/montage_window.m'             % shared: montage frame window
    };
    
    allModulesFound = true;
    for i = 1:size(modules, 1)
        moduleName = modules{i, 1};
        expectedPath = modules{i, 2};
        
        if exist(moduleName, 'file')
            fprintf('  ✓ %s\n', moduleName);
        else
            fprintf('  ✗ %s (expected at %s)\n', moduleName, expectedPath);
            allModulesFound = false;
        end
    end
    
    % Test key functionality
    if allModulesFound
        fprintf('\nTesting core functionality:\n');
        
        try
            % Test CSV loader
            loader = csv_loader_v2();
            fprintf('  ✓ CSV loader initialized\n');
            
            % Test configuration (now includes corrected Schmitt parameters)
            config = tracenorm_config();
            fprintf('  ✓ Configuration loaded:\n');
            fprintf('    - Corrected Schmitt: %.1f/%.1fσ (includes outliers in noise)\n', ...
                config.corrected_schmitt.upper_threshold_sigma, ...
                config.corrected_schmitt.lower_threshold_sigma);
            fprintf('    - Legacy Schmitt: %.1f/%.1fσ (excludes outliers from noise)\n', ...
                config.event_detection.upper_threshold_sigma, ...
                config.event_detection.lower_threshold_sigma);
            fprintf('    - Noise exclusion window: %d frames\n', ...
                config.corrected_schmitt.noise_exclusion_window);
            fprintf('    - Minimum event duration: %d frames\n', ...
                config.corrected_schmitt.min_event_duration);
            
            % Test main pipeline function (without running it)
            if exist('main_pipeline', 'file')
                fprintf('  ✓ Main pipeline accessible\n');
            end
            
            % Test corrected Schmitt detector
            if exist('corrected_schmitt_detector', 'file')
                fprintf('  ✓ Corrected Schmitt trigger detector available\n');
            end
            
            % Test legacy detector for comparison
            if exist('pure_schmitt_trigger_detector', 'file')
                fprintf('  ✓ Legacy Schmitt trigger detector available (for comparison)\n');
            end
            
            % Test quality assessor
            if exist('quality_assessor', 'file')
                fprintf('  ✓ Quality assessor module available\n');
            end
            
            fprintf('\n=== Setup Complete ===\n');
            fprintf('Ready to run:\n');
            fprintf('  results = main_pipeline(folder_path);  %% Uses corrected Schmitt by default\n');
            fprintf('  results = main_pipeline(folder_path, struct(''useCorrectedSchmitt'', false));  %% Uses legacy\n');
            fprintf('\nNew Features:\n');
            fprintf('  • Corrected Schmitt trigger (includes outliers in noise calculation)\n');
            fprintf('  • Proper sustained event exclusion (≥7 frames)\n');
            fprintf('  • 3-frame validation period for biological events\n');
            fprintf('  • Backward compatibility with legacy detector\n');
            fprintf('  • Enhanced noise estimation and thresholding\n');
            
        catch ME
            fprintf('  ✗ Error testing functionality: %s\n', ME.message);
            allModulesFound = false;
        end
    end
    
    if ~allModulesFound
        fprintf('\n=== Setup Issues Found ===\n');
        fprintf('Please check that you are in the SRC directory and all modules exist.\n');
        fprintf('Expected folder structure:\n');
        fprintf('  SRC/\n');
        fprintf('    ├── main_pipeline.m\n');
        fprintf('    ├── config/\n');
        fprintf('    │   └── tracenorm_config.m\n');
        fprintf('    ├── tracenorm/\n');
        fprintf('    │   ├── baseline_detector.m\n');
        fprintf('    │   ├── dfof_calculator.m\n');
        fprintf('    │   ├── corrected_schmitt_detector.m     ← NEW: Proper noise estimation\n');
        fprintf('    │   ├── pure_schmitt_trigger_detector.m  ← Legacy for comparison\n');
        fprintf('    │   └── quality_assessor.m\n');
        fprintf('    ├── visualization/\n');
        fprintf('    │   └── baseline_plotter_clean.m\n');
        fprintf('    ├── tests/\n');
        fprintf('    │   └── compare_schmitt_implementations.m  ← NEW: Compare methods\n');
        fprintf('    └── io/\n');
        fprintf('        └── csv_loader_v2.m\n');
    else
        % Display algorithm summary
        fprintf('\n=== Algorithm Summary ===\n');
        fprintf('CORRECTED Schmitt Trigger Logic:\n');
        fprintf('  1. Include outliers in noise calculation (they ARE noise)\n');
        fprintf('  2. Exclude sustained biological events (≥7 frames) from noise\n');
        fprintf('  3. Calculate thresholds: 3.5σ upper, 1.5σ lower\n');
        fprintf('  4. Event starts when crossing above upper threshold\n');
        fprintf('  5. Validate: must stay above lower threshold for ≥3 frames\n');
        fprintf('  6. Event ends when dropping below lower threshold\n');
        fprintf('  7. Merge nearby events (≤2 frame gap)\n');
        fprintf('\nKey Difference from Legacy:\n');
        fprintf('  Legacy: Excludes outliers from noise → UNDERESTIMATES noise\n');
        fprintf('  Corrected: Includes outliers in noise → REALISTIC noise estimation\n');
    end
end