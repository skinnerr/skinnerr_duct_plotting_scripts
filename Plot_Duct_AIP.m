%%%%%%
% 
% Plots the pressure recovery (PR) in a duct 
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
dirs{end+1} = 'TestDir1/AIP_Plot';
dirs{end+1} = 'TestDir2/AIP_Plot';
dirs{end+1} = '/home/ryan/Documents/CU_Boulder/Research/NGC/NGC_Ryan_Plotting_Scripts/TestDir3';
dirs{end+1} = '/home/ryan/Documents/CU_Boulder/Research/NGC/NGC_Ryan_Plotting_Scripts/TestDir4';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
% Sanitize directory inputs.
%%%

dirs = Sanitize_Paths(dirs);

%%%
% Load data from directories.
%%%

%   Note: Data from a specific line is indexed by dir_data{dir_index}(file_index).field
%   Note: Matlab allows you to say str='my_field' and access it by my_struct.(my_field)

dir_data = cell(length(dirs),1);
for i = 1:length(dirs)
    dir_data{i} = Load_CFD_Directory(dirs{i},'AIP');
end

% Sort data by y-coordinate, as it is saved in node-number order or somesuch.
sort_field = 'y';
for dir_i = 1:length(dir_data)
    for file_i = 1:length(dir_data{dir_i})
        % Sort sort_field data, and save the sort order.
        [dir_data{dir_i}(file_i).(sort_field), ...
                                 sort_order] = sort(dir_data{dir_i}(file_i).y);
        % Sort all other field data, using saved sort order.
        fields = fieldnames(dir_data{dir_i}(file_i));
        for i = 1:length(fields)
            field = fields{i};
            if any(strcmp(field,{sort_field,'ref','filename'}))
                continue
            end
            dir_data{dir_i}(file_i).(field) = ...
                dir_data{dir_i}(file_i).(field)(sort_order);
        end
    end
end

% Create new field for non-dimensionalized ertical position y/h.
duct_max_y =  0.1121537;
duct_min_y = -0.0333375;
assert(duct_max_y > duct_min_y, 'Duct min/max are out of order.');
for dir_i = 1:length(dir_data)
    for file_i = 1:length(dir_data{dir_i})
        dir_data{dir_i}(file_i).yoverh = ...
            (dir_data{dir_i}(file_i).y - duct_min_y) / (duct_max_y - duct_min_y);
    end
end

%%% TODO: The following is a HACK
% Normalize PR so it has a maximum value of 1.
for dir_i = 1:length(dir_data)
    for file_i = 1:length(dir_data{dir_i})
        dir_data{dir_i}(file_i).pr = ...
            dir_data{dir_i}(file_i).pr / max(dir_data{dir_i}(file_i).pr);
    end
end

% Determine plot order based on average z-coordinates of data files.
z_avg = cell(1,length(dir_data));
for dir_i = 1:length(dir_data)
    for file_i = 1:length(dir_data{dir_i})
        z_avg{dir_i}(file_i) = mean(dir_data{dir_i}(file_i).z);
    end
    [~,sort_order] = sort(z_avg{dir_i});
    dir_data{dir_i} = dir_data{dir_i}(sort_order);
end

% Set up plot size and initialize figure/axes.
fig_width  = 450;
fig_height = 700;
figure('Position',Centered_Figure_Position(fig_width,fig_height));
hax = axes();
hold on;

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
        
        % Cycle through the line styles, one for each file.
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
    % Reset color plotting order for the next directory.
    set(hax,'ColorOrderIndex',1);
    
end

xlabel('Pressure Recovery');
ylabel('y/h');
hleg = legend('show');
set(hleg,'Location','west');
hold off;






