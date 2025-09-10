%% CHECK_SETUP - Verify baseline calculation setup
% Run this script to check if all required files are in place

clear; clc;

fprintf('=== Baseline Calculation Setup Check ===\n');

% Check current directory
current_dir = pwd;
fprintf('Current directory: %s\n', current_dir);

% Expected file structure
expected_files = {
    'config/tracenorm_config.m',
    'tracenorm/baseline_detector.m', 
    'tracenorm/dfof_calculator.m',
    'visualization/baseline_plotter.m',
    'io/csv_loader_v2.m',
    'main_pipeline.m'
};

fprintf('\nChecking required files:\n');
all_files_found = true;

for i = 1:length(expected_files)
    file_path = expected_files{i};
    if exist(file_path, 'file')
        fprintf('  ✓ %s\n', file_path);
    else
        fprintf('  ✗ %s (MISSING)\n', file_path);
        all_files_found = false;
    end
end

% Add paths and check function accessibility
if all_files_found
    fprintf('\nAdding paths...\n');
    addpath(genpath(current_dir));
    
    % Check function accessibility
    fprintf('\nChecking function accessibility:\n');
    functions_to_check = {
        'tracenorm_config',
        'baseline_detector', 
        'dfof_calculator',
        'baseline_plotter',
        'csv_loader_v2'
    };
    
    all_functions_accessible = true;
    for i = 1:length(functions_to_check)
        func_name = functions_to_check{i};
        if exist(func_name, 'file')
            fprintf('  ✓ %s accessible\n', func_name);
        else
            fprintf('  ✗ %s NOT accessible\n', func_name);
            all_functions_accessible = false;
        end
    end
    
    % Test configuration loading
    if all_functions_accessible
        fprintf('\nTesting configuration:\n');
        try
            config = tracenorm_config();
            fprintf('  ✓ Configuration loaded successfully\n');
            fprintf('    Frame rate: %d Hz\n', config.frame_rate);
            fprintf('    Rolling window: %.1f s (%d frames)\n', ...
                config.rolling_window_sec, config.rolling_window_frames);
            fprintf('    Outlier threshold: %.1f σ\n', config.outlier_threshold_sigma);
        catch ME
            fprintf('  ✗ Configuration loading failed: %s\n', ME.message);
            all_functions_accessible = false;
        end
    end
    
    % Final status
    if all_functions_accessible
        fprintf('\n✅ SETUP COMPLETE - Ready to test baseline calculation!\n');
        fprintf('\nNext step: Update test_baseline_single_file.m with your CSV file path\n');
    else
        fprintf('\n❌ SETUP INCOMPLETE - Function accessibility issues\n');
    end
    
else
    fprintf('\n❌ SETUP INCOMPLETE - Missing required files\n');
    fprintf('\nExpected directory structure:\n');
    fprintf('  SRC/\n');
    fprintf('    ├── main_pipeline.m\n');
    fprintf('    ├── config/\n');
    fprintf('    │   └── tracenorm_config.m\n');
    fprintf('    ├── tracenorm/\n');
    fprintf('    │   ├── baseline_detector.m\n');
    fprintf('    │   └── dfof_calculator.m\n');
    fprintf('    ├── visualization/\n');
    fprintf('    │   └── baseline_plotter.m\n');
    fprintf('    └── io/\n');
    fprintf('        └── csv_loader_v2.m\n');
    
    fprintf('\nPlease ensure all files are in place and run from the SRC directory.\n');
end

fprintf('\n=== Setup Check Complete ===\n');