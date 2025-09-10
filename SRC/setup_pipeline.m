function setup_pipeline()
    % SETUP_PIPELINE - Initialize the fluorescent imaging analysis environment
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
        'main_pipeline', 'main_pipeline.m'
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
            
            % Test main pipeline function (without running it)
            if exist('main_pipeline', 'file')
                fprintf('  ✓ Main pipeline accessible\n');
            end
            
            fprintf('\n=== Setup Complete ===\n');
            fprintf('Ready to run: results = main_pipeline(folder_path);\n');
            
        catch ME
            fprintf('  ✗ Error testing functionality: %s\n', ME.message);
            allModulesFound = false;
        end
    end
    
    if ~allModulesFound
        fprintf('\n=== Setup Issues Found ===\n');
        fprintf('Please check that you are in the SRC directory and all modules exist.\n');
        fprintf('Current folder structure should be:\n');
        fprintf('  SRC/\n');
        fprintf('    ├── main_pipeline.m\n');
        fprintf('    ├── io/\n');
        fprintf('    │   └── csv_loader_v2.m\n');
        fprintf('    └── (other subfolders)\n');
    end
end