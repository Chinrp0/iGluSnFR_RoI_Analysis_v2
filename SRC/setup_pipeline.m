function setup_pipeline()
    % SETUP_PIPELINE - Initialize the fluorescent imaging analysis environment
    % UPDATED: Now verifies new modular event detection components
    % Run this once per MATLAB session to configure paths and dependencies
    
    fprintf('=== Setting up Fluorescent Imaging Pipeline ===\n');
    
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
        'schmitt_event_detector', 'tracenorm/schmitt_event_detector.m';    % NEW
        'quality_assessor', 'tracenorm/quality_assessor.m';                % NEW
        'baseline_plotter', 'visualization/baseline_plotter.m'
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
            
            % Test configuration (now includes Schmitt trigger parameters)
            config = tracenorm_config();
            fprintf('  ✓ Configuration loaded (Schmitt trigger: %.1f/%.1fσ)\n', ...
                config.event_detection.upper_threshold_sigma, ...
                config.event_detection.lower_threshold_sigma);
            
            % Test main pipeline function (without running it)
            if exist('main_pipeline', 'file')
                fprintf('  ✓ Main pipeline accessible\n');
            end
            
            % Test new modular components
            if exist('schmitt_event_detector', 'file')
                fprintf('  ✓ Schmitt trigger event detector available\n');
            end
            
            if exist('quality_assessor', 'file')
                fprintf('  ✓ Quality assessor module available\n');
            end
            
            fprintf('\n=== Setup Complete ===\n');
            fprintf('Ready to run: results = main_pipeline(folder_path);\n');
            fprintf('New features: Schmitt trigger event detection, modular quality assessment\n');
            
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
        fprintf('    │   ├── schmitt_event_detector.m    ← NEW\n');
        fprintf('    │   └── quality_assessor.m          ← NEW\n');
        fprintf('    ├── visualization/\n');
        fprintf('    │   └── baseline_plotter.m\n');
        fprintf('    └── io/\n');
        fprintf('        └── csv_loader_v2.m\n');
    end
end