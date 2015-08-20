%%%%%%
% 
% Plots varts data saved by PHASTA. Supports...
%   - Multiple varts directories plotted with phase offsets.
%   - Phase averaging. 
% 
% Ryan Skinner, August 2015
% 
%%%

Set_Default_Plot_Properties();

%%%
% Initialize variables.
%%%

dirs = {};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Directories containing plot data. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirs{end+1} = 'varts_1';
dirs{end+1} = 'varts_2';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
% Sanitize directory inputs.
%%%

dirs = Sanitize_Paths(dirs);

%%%
% Load data from directories.
%%%

% Set up plot size and initialize figure/axes.
fig_width  = 1000;
fig_height = 300;
figure('Position',Centered_Figure_Position(fig_width,fig_height));
hax = axes();
hold on;

return

% Set up color order and line styles to use.
c_order = get(hax,'ColorOrder');
set(hax,'ColorOrder',c_order([5,7,1,2,3,4,6],:));
styles = {'-','--','-.',':'};

% Plot in order of z-coordinate as sorted above.
for dir_i = 1:length(dir_data);
    for file_i = 1:length(dir_data{dir_i})
        
        % If the slice is centered on z=0...
        if abs(z_avg{dir_i}(file_i)) < 1e-15
            % Set the slice center to zero identically.
            z_avg{dir_i}(file_i) = 0;
            % And color its line black.
            color_str = 'k';
        else
            color_str = '';
        end
        
        % Cycle through the line styles, one for each directory.
        style_str = styles{1+mod(dir_i-1,length(styles))};
        
        % Determine name to appear on legend.
        m2in = 39.3701;
        z_avg_inches = m2in * z_avg{dir_i}(file_i);
        dir_str = strsplit(dirs{dir_i},'/');
        if strcmp(dir_str{end-1},'AIP_Plot')
            dir_str = dir_str{end-2};
        else
            dir_str = dir_str{end-1};
        end
        display_name = sprintf('%s  %.1f in',dir_str,z_avg_inches);
        
        % Plot in order of z-coordinate.
        plot(dir_data{dir_i}(file_i).pr, ...
             dir_data{dir_i}(file_i).yoverh, ...
             [style_str,color_str], ...
             'DisplayName', display_name);
    end
    % Reset color plotting order for the next file.
    set(hax,'ColorOrderIndex',1);
    
end

xlabel('Pressure Recovery');
ylabel('y/h');
hleg = legend('show');
set(hleg,'Location','west');
hold off;






