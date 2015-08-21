%%%%%%
% 
% Plots varts data saved by PHASTA. Supports phase averaging.
% 
% Ryan Skinner, August 2015
% 
%%%

Set_Default_Plot_Properties();

%%%
% Initialize variables.
%%%

% Physical constants.
gamma = 1.4;
R     = 288.294;

% Line styles to cycle through for each directory.
styles = {'-','--','-.',':'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Directories containing plot data. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirs = {};
dir_names = {};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirs{end+1} = 'varts_1'; dir_names{end+1} = 'varts_1';
% dirs{end+1} = 'varts_2'; dir_names{end+1} = 'varts_2';
% dirs{end+1} = 'varts_3'; dir_names{end+1} = 'varts_3';
% dirs{end+1} = 'varts_4'; dir_names{end+1} = 'varts_4';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase-average settings. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
avg_settings = {};
do_phase_avg = false;
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
% Sanitize inputs.
%%%

% Directories.
dirs = Sanitize_Paths(dirs);

% Directory names.
if length(dirs) ~= length(dir_names)
    error('Number of paths is not equal to number of named directories.');
end

% Phase-averaging.
if ~do_phase_avg
    avg_settings = {};
end
if do_phase_avg && (length(avg_settings) ~= length(dirs))
    error('Phase-average settings are not specified for all directories requested.');
end

% Probe points.
if length(probeIDs) < 1
    error('No probe IDs were specified.');
end

%%%
% Load data from directories.
%%%

varts = cell(1,length(dirs));

for dir_i = 1:length(dirs)
    
    [dt, ts, xyz, p, u, T, nu] = Load_Varts_Directory(dirs{dir_i});
    
    % Retain only requested probe points.
    p  =  p(probeIDs,:);
    u  =  u(probeIDs,:,:);
    T  =  T(probeIDs,:);
    nu = nu(probeIDs,:);
    
    % Save data to structure
    varts{dir_i}.dt  = dt;
    varts{dir_i}.ts  = ts;
    varts{dir_i}.xyz = xyz;
    varts{dir_i}.p   = p;
    varts{dir_i}.u   = u;
    varts{dir_i}.T   = T;
    varts{dir_i}.nu  = nu;
    
end

%%%
% Sanitize dt inputs to plotting functions.
%%%

% Set dt to NaN if any directory has non-uniform dt or if two directory dt's don't match.
dt = nan(1,length(dirs));
for dir_i = 1:length(dirs)
    if abs(mean(varts{dir_i}.dt) - varts{dir_i}.dt(1)) < 1e-15
        dt(dir_i) = varts{dir_i}.dt(1);
    end
end
if mean(dt) ~= dt(1)
    dt = nan(1:length(dirs));
end

%%%
% Loop through the directories and plot their contents.
%%%

for dir_i = 1:length(dirs)
    
    % Pass an empty structure of settings if phase-averaging is disabled.
    settings = {};
    if ~isempty(avg_settings)
        settings = avg_settings{dir_i};
    end
    
    % Set properties for line names and styles.
    name = dir_names{dir_i};
    style_i = 1 + mod(dir_i-1, length(styles));
    style   = styles{style_i};
    
    v = varts{dir_i};
    % Calculate extra fields.
    tmp_ts = varts{dir_i}.ts;
    speed  = sqrt(sum(v.u(:, :, :) .* v.u(:, :, :), 3));
    mach   =  sqrt(sum(v.u(:, :, :) .* v.u(:, :, :), 3) ./ (gamma*R*v.T(:, :)));
    rho    = v.p(:,:,:) ./ (R*v.T(:, :));
    dynp   = 0.5 * rho .* speed.^2;
    totp   = v.p(:,:) + dynp;
    
    %%%
    % Plot time histories.
    %%%
    
    Plot_Varts_Data(dt, tmp_ts, v.u(:,:,1), 'u1', settings, name, probeIDs, style);
    Plot_Varts_Data(dt, tmp_ts, v.u(:,:,2), 'u2', settings, name, probeIDs, style);
    Plot_Varts_Data(dt, tmp_ts, v.u(:,:,3), 'u3', settings, name, probeIDs, style);
    Plot_Varts_Data(dt, tmp_ts, v.T,         'T', settings, name, probeIDs, style);
    Plot_Varts_Data(dt, tmp_ts, speed,    'magu', settings, name, probeIDs, style);
    Plot_Varts_Data(dt, tmp_ts, mach,        'M', settings, name, probeIDs, style);
    Plot_Varts_Data(dt, tmp_ts, rho,       'rho', settings, name, probeIDs, style);
    Plot_Varts_Data(dt, tmp_ts, v.p,         'p', settings, name, probeIDs, style);
    Plot_Varts_Data(dt, tmp_ts, dynp,     'dynp', settings, name, probeIDs, style);
    Plot_Varts_Data(dt, tmp_ts, totp,     'totp', settings, name, probeIDs, style);
    
    %%%
    % Plot Fourier transforms.
    %%%
    
    for probe_i = 1:length(probeIDs)
        fig_number = 100 + probeIDs(probe_i);
        figure(fig_number);
        hold on;
        time_step = varts{dir_i}.dt(1);
        [freq, spectral_power] = Fourier_Transform(v.u(:, probe_i, 2), time_step);
        loglog(freq, spectral_power, 'DisplayName', name);
        title(['Power Spectrum for Probe ', num2str(probeIDs(probe_i))]);
        xlabel('Frequency (Hz)');
        ylabel('Spectral Power');
    end
    
end


















