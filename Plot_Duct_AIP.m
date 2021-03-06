%%%%%%
% 
% Plots the pressure recovery (PR) in a duct 
% 
%%%

Set_Default_Plot_Properties();

%%%
% Initialize variables.
%%%

m2in = 39.3701;
h_duct = 0.1524;
duct_max_y =  0.1121537;
duct_min_y = -0.0333375;

x_limits = [0.78,1.01];
y_limits = [0.00,1.00];

dirs = {};
dir_names = {};
styles = {};
markers = {};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Directories containing plot data. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dirs{end+1} = 'series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-DDES-Baseline/ts19.7k-54.1k_CpPressureRecoveryPlots'; dir_names{end+1} = 'Baseline'; styles{end+1} = '-'; markers{end+1} = '';

% dirs{end+1} = 'series14.4.2/meshing-Coarse2-Par/Run-DDES-100Hz-ExpMatch150803-LBOff-P101.3k/ts426.1k-460.1k_CpPressureRecoveryPlots'; dir_names{end+1} = 'mdot = 0.90%';
% dirs{end+1} = 'series14.4.2/meshing-Coarse2-Par/Run-DDES-100Hz-ExpMatch150803-LBOff/ts363.6k-375.6k_CpPressureRecoveryPlots'; dir_names{end+1} = 'mdot = 1.18%';

% dirs{end+1} = 'series14.4.2/meshing-Coarse2-Par/Run-RANS-Steady-mdot1.25/ts189k-193k_data/'; dir_names{end+1} = 'RANS 1.15% (74a)'; styles{end+1} = ':'; markers{end+1} = 'o';
% dirs{end+1} = 'series14.4.2/meshing-Coarse2-Par/Run-RANS-Steady-mdot1.5/ts191k-195k_data/'; dir_names{end+1} = 'RANS 1.55% (75a)'; styles{end+1} = ':'; markers{end+1} = 'x';
dirs{end+1} = 'series14.4.2/meshing-Coarse2-Par/Run-DDES-Steady-mdot1.25/77a_Data'; dir_names{end+1} = 'DDES 1.12% (77a)'; styles{end+1} = '--'; markers{end+1} = 'o';
dirs{end+1} = 'series14.4.2/meshing-Coarse2-Par/Run-DDES-Steady-mdot1.5/78a_Data'; dir_names{end+1} = 'DDES 1.57% (78a)'; styles{end+1} = '--'; markers{end+1} = 'x';

% dirs{end+1} = 'series14.5.1/A0/A0-Run-RANS-Steady-mdot1.25/ts105k-115k_data'; dir_names{end+1} = 'RANS 1.24% (73a)'; styles{end+1} = ':'; markers{end+1} = 'o';
% dirs{end+1} = 'series14.5.1/A0/A0-Run-RANS-Steady-mdot1.5/ts105k-119k_data'; dir_names{end+1} = 'RANS 1.83% (72b)'; styles{end+1} = ':'; markers{end+1} = 'x';
% dirs{end+1} = 'series14.5.1/A0/A0-Run-DDES-Steady-mdot1.25/79a_Data'; dir_names{end+1} = 'DDES 1.22% (79a)'; styles{end+1} = '--'; markers{end+1} = 'o';
% dirs{end+1} = 'series14.5.1/A0/A0-Run-DDES-Steady-mdot1.5/80a_Data'; dir_names{end+1} = 'DDES 1.80% (80a)'; styles{end+1} = '--'; markers{end+1} = 'x';
% dirs{end+1} = 'series14.5.1/A0/A0-Run-RANS-Steady-mdot1.25-LB1.72/81a_Data'; dir_names{end+1} = 'RANS 1.26% / 1.95% (81a)'; styles{end+1} = ':'; markers{end+1} = 'o';
% dirs{end+1} = 'series14.5.1/A0/A0-Run-RANS-Steady-mdot1.5-LB1.72/82a_Data'; dir_names{end+1} = 'RANS 2.32% / 1.89% (82a)'; styles{end+1} = ':'; markers{end+1} = 'x';
% dirs{end+1} = 'series14.5.1/A0/A0-Run-DDES-Steady-mdot1.25-LB1.72/83a_Data'; dir_names{end+1} = 'DDES 1.27% / 1.90% (83a)'; styles{end+1} = '--'; markers{end+1} = 'o';
% dirs{end+1} = 'series14.5.1/A0/A0-Run-DDES-Steady-mdot1.5-LB1.72/84a_Data'; dir_names{end+1} = 'DDES 2.33% / 1.91% (84a)'; styles{end+1} = '--'; markers{end+1} = 'x';

directory = 'AIP_Plot';
for i = 1:length(dirs)
    dirs{i} = ['/projects/ngc/probePoints/', dirs{i}, '/', directory, '/'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

styles = {'--',':'};
markers = {'x','o','o','x'};

%%%
% Sanitize directory inputs.
%%%

dirs = Slash_Terminate_Paths(dirs);

% Directory names.
if length(dirs) ~= length(dir_names)
    error('Number of CFD paths is not equal to number of named directories.');
end

% Line styles.
if length(dirs) ~= length(styles)
    fprintf('Warning: Number of CFD paths is not equal to number of line styles.');
end

%%%
% Load data from directories.
%%%

% Note: Data from a specific line is indexed by dir_data{dir_index}(file_index).field
% Note: Matlab allows you to say str='my_field' and access it by my_struct.(my_field)

dir_data = cell(length(dirs),1);
for i = 1:length(dirs)
    dir_data{i} = Load_CFD_Directory(dirs{i},'AIP');
end

% Sort data by y-coordinate, as it is saved in node-number order or somesuch.
sort_field = 'y';
for dir_i = 1:length(dir_data)
    for file_i = 1:length(dir_data{dir_i})
        
        % Sort sort_field data, and save the sort order.
        [dir_data{dir_i}(file_i).(sort_field), sort_order] = ...
            sort(dir_data{dir_i}(file_i).(sort_field));
        
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

% Create new field for non-dimensionalized vertical position y/h.
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

% Plot in order of z-coordinate as sorted above.
for dir_i = 1:length(dir_data);
    for file_i = 1:length(dir_data{dir_i})
        
        % Color centerline black.
        color_str = '';
        if abs(z_avg{dir_i}(file_i)) < 10*eps
            color_str = 'k';
            % Set the slice center to zero identically.
            z_avg{dir_i}(file_i) = 0;
        end
        
        % Cycle through the line styles, one for each file.
         style_str =  styles{1+mod(dir_i-1,length(styles))};
        marker_str = markers{1+mod(dir_i-1,length(markers))};
        
        % Localize some variables for convenience.
        yoverh = dir_data{dir_i}(file_i).yoverh;
        pr = dir_data{dir_i}(file_i).pr;
        
        % Plot data lines.
        coi = get(hax,'ColorOrderIndex');
        plot(pr, yoverh, [style_str,color_str]);
        
        % Plot data markers separately so we can control their spacing.
        if ~isempty(marker_str)
            set(hax,'ColorOrderIndex',coi);
            n_markers = 30;
            marker_yoverh = linspace(min(yoverh),max(yoverh), n_markers);
            marker_pr = linterp(yoverh, pr, marker_yoverh);
            plot(marker_pr, marker_yoverh, [marker_str,color_str]);
        end
        
        % Set limits.
        xlim(x_limits);
        ylim(y_limits);
    end
    % Reset color plotting order for the next directory.
    set(hax,'ColorOrderIndex',1);
    
end

% Finally, label the lines.
% We cannot do this above, since for each data set, we plot lines and
%   data markers separately.
legend_entries = [];
for dir_i = 1:length(dir_data);
    for file_i = 1:length(dir_data{dir_i})
        
        % Color centerline black.
        color_str = '';
        if abs(z_avg{dir_i}(file_i)) < 10*eps
            color_str = 'k';
            % Set the slice center to zero identically.
            z_avg{dir_i}(file_i) = 0;
        end
        
        % Determine name to appear on legend.
        dir_str = dir_names{dir_i};
        %z_avg_inches = m2in * z_avg{dir_i}(file_i);
        %display_name = sprintf('%s (%.1f in)', dir_str, z_avg_inches);
        z_avg_over_h = z_avg{dir_i}(file_i) / h_duct;
        display_name = sprintf('%s (%.2f z/h)', dir_str, z_avg_over_h);
        
        % Cycle through the line styles, one for each file.
         style_str =  styles{1+mod(dir_i-1,length(styles))};
        marker_str = markers{1+mod(dir_i-1,length(markers))};
        
        % Plot a dummy line to appear in the legend.
        legend_entries(end+1) = plot(x_limits(1)-1, y_limits(1)-1, ...
                                     [style_str, marker_str, color_str], ...
                                     'DisplayName', display_name);
        
    end
    % Reset color plotting order for the next directory.
    set(hax,'ColorOrderIndex',1);
end

xlabel('Pressure Recovery');
ylabel('y/h');
hleg = legend(legend_entries);
set(hleg,'Location','southwest');
hold off;






