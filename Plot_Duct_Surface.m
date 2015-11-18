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
L_ramp = 0.217342;

dirs = {};
dir_names = {};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Directories containing plot data. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dirs{end+1} = 'series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-DDES-Baseline/ts19.7k-54.1k_CpPressureRecoveryPlots'; dir_names{end+1} = 'Baseline';

% dirs{end+1} = 'series14.4.2/meshing-Coarse2-Par/Run-DDES-100Hz-ExpMatch150803-LBOff-P101.3k/ts426.1k-460.1k_CpPressureRecoveryPlots'; dir_names{end+1} = 'mdot = 0.90%';
% dirs{end+1} = 'series14.4.2/meshing-Coarse2-Par/Run-DDES-100Hz-ExpMatch150803-LBOff/ts363.6k-375.6k_CpPressureRecoveryPlots'; dir_names{end+1} = 'mdot = 1.18%';

% dirs{end+1} = 'series14.4.2/meshing-Coarse2-Par/Run-RANS-Steady-mdot1.25/ts189k-193k_data/'; dir_names{end+1} = 'RANS 1.15% (74a)';
% dirs{end+1} = 'series14.4.2/meshing-Coarse2-Par/Run-RANS-Steady-mdot1.5/ts191k-195k_data/'; dir_names{end+1} = 'RANS 1.55% (75a)';
% dirs{end+1} = 'series14.4.2/meshing-Coarse2-Par/Run-DDES-Steady-mdot1.25/77a_Data'; dir_names{end+1} = 'DDES 1.12% (77a)';
% dirs{end+1} = 'series14.4.2/meshing-Coarse2-Par/Run-DDES-Steady-mdot1.5/78a_Data'; dir_names{end+1} = 'DDES 1.57% (78a)';

% dirs{end+1} = 'series14.5.1/A0/A0-Run-RANS-Steady-mdot1.25/ts105k-115k_data'; dir_names{end+1} = 'RANS 1.24% (73a)';
% dirs{end+1} = 'series14.5.1/A0/A0-Run-RANS-Steady-mdot1.5/ts105k-119k_data'; dir_names{end+1} = 'RANS 1.83% (72b)';
% dirs{end+1} = 'series14.5.1/A0/A0-Run-DDES-Steady-mdot1.25/79a_Data'; dir_names{end+1} = 'DDES 1.22% (79a)';
% dirs{end+1} = 'series14.5.1/A0/A0-Run-DDES-Steady-mdot1.5/80a_Data'; dir_names{end+1} = 'DDES 1.80% (80a)';
% dirs{end+1} = 'series14.5.1/A0/A0-Run-RANS-Steady-mdot1.25-LB1.72/81a_Data'; dir_names{end+1} = 'RANS 1.26% / 1.95% (81a)';
% dirs{end+1} = 'series14.5.1/A0/A0-Run-RANS-Steady-mdot1.5-LB1.72/82a_Data'; dir_names{end+1} = 'RANS 2.32% / 1.89% (82a)';
dirs{end+1} = 'series14.5.1/A0/A0-Run-DDES-Steady-mdot1.25-LB1.72/83a_Data'; dir_names{end+1} = 'DDES 1.27% / 1.90% (83a)';
dirs{end+1} = 'series14.5.1/A0/A0-Run-DDES-Steady-mdot1.5-LB1.72/84a_Data'; dir_names{end+1} = 'DDES 2.33% / 1.91% (84a)';

% directory = 'dpdxPlot_Upper'; custom_ylim = [-0.5, 0.8];
directory = 'dpdxPlot_Lower'; custom_ylim = [-0.1, 0.7];

for i = 1:length(dirs)
    dirs{i} = ['/projects/ngc/probePoints/', dirs{i}, '/', directory, '/'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
% Sanitize directory inputs.
%%%

dirs = Slash_Terminate_Paths(dirs);

% Directory names.
if length(dirs) ~= length(dir_names)
    error('Number of CFD paths is not equal to number of named directories.');
end

%%%
% Load data from directories.
%%%

% Note: Data from a specific line is indexed by dir_data{dir_index}(file_index).field
% Note: Matlab allows you to say str='my_field' and access it by my_struct.(my_field)

dir_data = cell(length(dirs),1);
for i = 1:length(dirs)
    dir_data{i} = Load_CFD_Directory(dirs{i});
end

% Sort data by x-coordinate, as it is saved in node-number order or somesuch.
sort_field = 'x';
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

% Create new field for non-dimensionalized stream-wise position x/L.
for dir_i = 1:length(dir_data)
    for file_i = 1:length(dir_data{dir_i})
        dir_data{dir_i}(file_i).xoverL = dir_data{dir_i}(file_i).x / L_ramp;
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
fig_width  = 900;
fig_height = 350;
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
        if abs(z_avg{dir_i}(file_i)) < 1e-5
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
        dir_str = dir_names{dir_i};
%         z_avg_inches = m2in * z_avg{dir_i}(file_i);
%         display_name = sprintf('%s (%.1f in)', dir_str, z_avg_inches);
        z_avg_over_h = z_avg{dir_i}(file_i) / h_duct;
        display_name = sprintf('%s (%.2f z/h)', dir_str, z_avg_over_h);
        
        % Plot in order of z-coordinate.
        plot(dir_data{dir_i}(file_i).xoverL, ...
             dir_data{dir_i}(file_i).cp, ...
             [style_str,color_str], ...
             'DisplayName', display_name);
    end
    % Reset color plotting order for the next directory.
    set(hax,'ColorOrderIndex',1);
    
end

xlabel('x/L');
ylabel('Cp');
xlim([0,1.1]);
ylim(custom_ylim);
hleg = legend('show');
set(hleg,'Location','southeast');
hold off;






