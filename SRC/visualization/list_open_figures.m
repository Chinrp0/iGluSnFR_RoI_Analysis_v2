function list_open_figures()
    % LIST_OPEN_FIGURES - Show all open figure windows with their names
    %
    % Usage:
    %   list_open_figures()
    %
    % This helps diagnose why 26 figures are being created
    
    fprintf('\n=== OPEN FIGURES DIAGNOSTIC ===\n');
    
    % Get all figure handles
    fig_handles = findall(0, 'Type', 'figure');
    num_figs = length(fig_handles);
    
    fprintf('Total open figures: %d\n\n', num_figs);
    
    if num_figs == 0
        fprintf('No figures are currently open.\n');
        return;
    end
    
    % Sort by figure number
    fig_numbers = zeros(num_figs, 1);
    for i = 1:num_figs
        fig_numbers(i) = fig_handles(i).Number;
    end
    [~, sort_idx] = sort(fig_numbers);
    fig_handles = fig_handles(sort_idx);
    
    % List each figure
    fprintf('Figure List:\n');
    fprintf('%-8s %-50s %-15s\n', 'Number', 'Name', 'Size');
    fprintf('%s\n', repmat('-', 1, 80));
    
    for i = 1:num_figs
        fig = fig_handles(i);
        fig_num = fig.Number;
        fig_name = fig.Name;
        
        % Get figure size
        fig_pos = fig.Position;
        fig_size = sprintf('[%d×%d]', fig_pos(3), fig_pos(4));
        
        % Truncate long names
        if length(fig_name) > 50
            fig_name = [fig_name(1:47), '...'];
        end
        
        fprintf('%-8d %-50s %-15s\n', fig_num, fig_name, fig_size);
    end
    
    fprintf('\n');
    
    % Group by name patterns
    fprintf('Grouped by Module:\n');
    fprintf('%s\n', repmat('-', 1, 80));
    
    % Define expected module patterns
    patterns = {
        'Event Frequency', 'Core frequency comparison';
        'Event Amplitude', 'Core amplitude comparison';
        'Active ROI', 'Core activity comparison';
        'File Summary', 'Core file summary';
        'Amplitude-Frequency', 'Amplitude-frequency correlation';
        'Event Timing Raster', 'Event timing raster';
        'ROI Recruitment', 'ROI recruitment curves';
        'Cumulative Events', 'Cumulative events over time';
        'ROI Trace', 'ROI trace comparison';
        'Inter-Event Interval', 'IEI comparison';
        'Biological Variability', 'Biological variability';
        'WT vs R213W', 'Trace comparison';
    };
    
    for p = 1:size(patterns, 1)
        pattern = patterns{p, 1};
        module_name = patterns{p, 2};
        
        % Count figures matching this pattern
        count = 0;
        matching_figs = [];
        for i = 1:num_figs
            if contains(fig_handles(i).Name, pattern, 'IgnoreCase', true)
                count = count + 1;
                matching_figs(end+1) = fig_handles(i).Number;
            end
        end
        
        if count > 0
            fprintf('  %-40s: %d figure(s) [%s]\n', module_name, count, ...
                mat2str(matching_figs));
        end
    end
    
    % Check for unnamed/empty figures
    unnamed_count = 0;
    unnamed_figs = [];
    for i = 1:num_figs
        if isempty(fig_handles(i).Name)
            unnamed_count = unnamed_count + 1;
            unnamed_figs(end+1) = fig_handles(i).Number;
        end
    end
    
    if unnamed_count > 0
        fprintf('  %-40s: %d figure(s) [%s]\n', 'Unnamed/blank figures', ...
            unnamed_count, mat2str(unnamed_figs));
    end
    
    fprintf('\n=== END OF DIAGNOSTIC ===\n\n');
    
    % Suggest cleanup if too many figures
    if num_figs > 15
        fprintf('⚠ You have %d figures open, which may indicate duplicates.\n', num_figs);
        fprintf('To close all figures: close all\n');
        fprintf('To close specific figures: close([fig_numbers])\n\n');
    end
end
