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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Directories containing plot data. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirs = {};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dirs{end+1} = 'varts_1';
dirs{end+1} = 'varts_2';
% dirs{end+1} = 'varts_3';
% dirs{end+1} = 'varts_4';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase average settings. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
avg_settings = {};
do_phase_avg = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
s.start_ts = 61000;
s.end_ts   = 65000;
s.freq_Hz  = 300;
avg_settings{end+1} = s;
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Probe point selection. %
%%%%%%%%%%%%%%%%%%%%%%%%%%
probeIDs = [1, 261, 264, 267, 270, 273];
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
% Sanitize directory inputs.
%%%

dirs = Sanitize_Paths(dirs);

%%%
% Load data from directories.
%%%

varts = cell(1,length(dirs));

for dir_i = 1:length(dirs)
    
    [dt, t, xyz, p, u, T, nu] = Load_Varts_Directory(dirs{dir_i});
    
    % Retain only requested probe points.
    p  =  p(probeIDs,:);
    u  =  u(probeIDs,:,:);
    T  =  T(probeIDs,:);
    nu = nu(probeIDs,:);
    
    % Save data to structure
    varts{dir_i}.dt  = dt;
    varts{dir_i}.t   = t;
    varts{dir_i}.xyz = xyz;
    varts{dir_i}.p   = p;
    varts{dir_i}.u   = u;
    varts{dir_i}.T   = T;
    varts{dir_i}.nu  = nu;
    
end

%%%
% Phase-average data if requested.
%%%

if do_phase_avg
    if length(avg_settings) ~= length(dirs)
        error('Phase average settings are not specified for all directories requested.');
    end
    for dir_i = 1:length(dirs)
        varts{dir_i} = Phase_Average_Varts(varts{dir_i}, avg_settings{dir_i});
    end
end

%%%
% Plot data.
%%%

% Set up plot size and initialize figure/axes.
fig_width  = 1000;
fig_height = 300;
figure('Position',Centered_Figure_Position(fig_width,fig_height));
hax = axes();
hold on;

for dir_i = 1:length(dirs);
    for probe_i = 1:length(probeIDs)
        plot(varts{dir_i}.t, varts{dir_i}.u(probe_i,:,2));
    end
end

% Set up color order and line styles to use.
c_order = get(hax,'ColorOrder');
set(hax,'ColorOrder',c_order([5,7,1,2,3,4,6],:));
styles = {'-','--','-.',':'};