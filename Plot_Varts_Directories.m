%%%%%%
% 
% Plots varts data saved by PHASTA. Supports phase averaging.
% 
% Ryan Skinner, August 2015
% 
%%%

Set_Default_Plot_Properties();

% close all;

%%%
% Initialize variables.
%%%

% Physical constants.
gamma = 1.4;
R     = 288.294;

% Line styles to cycle through for each directory.
styles = {'-','--','-.',':'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Directories containing CFD data. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirs = {};
dir_names = {};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dirs{end+1} = '../series14.4.2/meshing-fine_blowerOff_unsteady-100Hz';
% dirs{end+1} = '../series14.4.2/meshing-fine_blowerOn';
% dirs{end+1} = '../series14.4.2/meshing-fine_blowerOff';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-TrueDDES-300Hz-LowBlow1.72';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-TrueDDES-300Hz-LowBlowOff';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-TrueDDES-NoBlowing';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-DDES-Steady-LBOff';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-DDES-300Hz-LowBlowOn/From128600';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-DDES-300Hz';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-DDES-600Hz';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-DDES-1000Hz';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-DDES-NoBlowing';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-DDES-NoBlowing/ForkedFrom';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-RANS-300Hz-prof2-LowBlowOff';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-RANS-300Hz-prof2-LowBlow1.72';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-RANS-300Hz-LowBlow1.72';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-RANS-300Hz-LowBlowOff';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-RANS-NoBlowing';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-RANS-MatchExpBlow300Hz';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-RANS-MatchExpBlow300HzTry2';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-DDES-300HzLB-ExpBlow';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-DDES-300HzLBOff-ExpBlow';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-DDES-300Hz-LowBlowOff-ExpBlow-UR';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-DDES-150Hz-ExpBlowHalfNumSlits-LowBlowOff';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-DDES-150Hz-ExpBlowHalfNumSlits-LowBlowOff-UR';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-DDES-300Hz-ExpMatch150803-LBOff';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-DDES-Steady-LBOff';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-DDES-Baseline';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-DDES-Baseline25k-recent';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-DDES-300Hz-prof2-LowBlow1.72';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-DDES-300Hz-prof2-LowBlow0.86';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-DDES-300Hz-prof2-LowBlowOff';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-DDES-600Hz-prof2-LowBlowOff';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-DDES-1000Hz-prof2-LowBlowOff';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-DDES-600Hz-LowBlowOff';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-DDES-300Hz-LowBlow1.72';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-DDES-300Hz-LowBlow1.72_Dt2pt5em6';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-DDES-300Hz-LowBlowOff';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-DDES-300Hz-ExpMatch150803-LBOff'; dir_names{end+1} = 'S14.4.2 A1 Trap 300Hz';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-DES97-300Hz-LowBlow1.72';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-DES97-300Hz-LowBlowOff';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-RANS-Baseline';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-RANS-300Hz-LowBlowOff';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-RANS-300Hz-LowBlow1.72';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-RANS-300Hz-prof2-LowBlowOff';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-DDES-100Hz-LowBlowOff';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-DDES-100Hz-LowBlow0.86';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-DDES-100Hz-LowBlow1.72';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-DDES-134Hz-LowBlowOff';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-DDES-134Hz-LowBlow0.86';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-DDES-134Hz-LowBlow1.72';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-DDES-200Hz-LowBlowOff';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-DDES-200Hz-LowBlow1.72';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par-waggingjet/Run-RANS-Start';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par-waggingjet/Run-DDES-Wag-300Hz-LowBlowOff';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par-waggingjet/Run-DDES-Wag-300Hz-LowBlowOff-HalfSlitWidth';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par-waggingjet-fullUBmesh-Coarse2/Run-RANS-StartingUp';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par-waggingjet-fullUBmesh-Coarse2/Run-DDES-Wag-300Hz-LowBlowOff';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par-waggingjet-fullUBmesh-Coarse2/Run-DDES-Wag-300Hz-ExpMatch150803-LowBlowOff'; dir_names{end+1} = 'S14.4.2 A0-Wag-Fine 300Hz';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par-waggingjet-fullUBmesh-Coarse4/Run-DDES-Wag-300Hz-LowBlowOff';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par-waggingjet-fullUBmesh-Coarse4/Run-DDES-Wag-300Hz-ExpMatch150803-LowBlowOff'; dir_names{end+1} = 'S14.4.2 A0-Wag-Coarse 300Hz';
% 
% dirs{end+1} = '../series14.5.1/A0/A0-Run-RANS-Start2';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-DDES-Base-KJ';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-DDES-Steady-LBOff';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-DDES-Steady-LBOff-UR';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-DDES-100Hz-LB0.86-KJ';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-DDES-100Hz-LB1.72-KJ';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-DDES-100Hz-LB1.72-KJ-UR';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-DDES-100Hz-LB1.72-KJ-NoNormalUR';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-DDES-100Hz-KJ';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-RANS-100Hz-LB1.72';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-RANS-100Hz-LB1.72-UR';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dirs{end+1} = 'varts_1'; dir_names{end+1} = 'varts_1';
% dirs{end+1} = 'varts_2'; dir_names{end+1} = 'varts_2';
% dirs{end+1} = 'varts_3'; dir_names{end+1} = 'varts_3';
% dirs{end+1} = 'varts_4'; dir_names{end+1} = 'varts_4';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase-average and offset settings. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
avg_settings_cfd = {};
DO_PHASE_AVG_CFD = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % A1-Iso trapezoidal expmatch
% s.start_ts = 62000;
% s.end_ts   = 68800;
% s.freq_Hz  = 300;
% avg_settings{end+1} = s;
% % Wag mesh coarse2 expmatch
% s.start_ts = 18135;
% s.end_ts   = 42500;
% s.freq_Hz  = 300;
% avg_settings{end+1} = s;
% % Wag mesh coarse4 expmatch
% s.start_ts = 13086;
% s.end_ts   = 35100;
% s.freq_Hz  = 300;
% avg_settings{end+1} = s;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test directory
s.start_ts = 61000;
s.end_ts   = 65000;
s.freq_Hz  = 300;
avg_settings_cfd{end+1} = s;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourier transform flag. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
DO_PLOT_FOURIER_CFD = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Probe point selection. %
%%%%%%%%%%%%%%%%%%%%%%%%%%
% probeIDs = [1, 261, 264, 267, 270, 273];
% probeIDs = [1, 261, 264, 267];
probeIDs = [267];
%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%
% Sanitize CFD inputs.
%%%

if length(dirs) < 1
    warning('No CFD directories specified.');
else

    % Directories.
    dirs = Slash_Terminate_Paths(dirs);

    % Directory names.
    if length(dirs) ~= length(dir_names)
        error('Number of CFD paths is not equal to number of named directories.');
    end

    % Phase-averaging.
    if ~DO_PHASE_AVG_CFD
        avg_settings_cfd = {};
    end
    if DO_PHASE_AVG_CFD && (length(avg_settings_cfd) ~= length(dirs))
        error('Phase-average settings are not specified for all CFD directories requested.');
    end

    % Probe points.
    if length(probeIDs) < 1
        error('No CFD probe IDs were specified.');
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
        if ~isempty(avg_settings_cfd)
            settings = avg_settings_cfd{dir_i};
        end

        % Set properties for line names and styles.
        name = dir_names{dir_i};
        style_i = 1 + mod(dir_i-1, length(styles));
        style   = styles{style_i};

        % Calculate extra fields.
        v = varts{dir_i};
        tmp_dt = dt(dir_i);
        tmp_ts = v.ts;
        speed  = sqrt(sum(v.u(:, :, :) .* v.u(:, :, :), 3));
        mach   = sqrt(sum(v.u(:, :, :) .* v.u(:, :, :), 3) ./ (gamma*R*v.T(:, :)));
        rho    = v.p(:,:,:) ./ (R*v.T(:, :));
        dynp   = 0.5 * rho .* speed.^2;
        totp   = v.p(:,:) + dynp;

        %%%
        % Plot time histories.
        %%%

        Plot_Varts_Data(tmp_dt, tmp_ts, v.u(:,:,1), 'u1', settings, name, probeIDs, style);
        Plot_Varts_Data(tmp_dt, tmp_ts, v.u(:,:,2), 'u2', settings, name, probeIDs, style);
        Plot_Varts_Data(tmp_dt, tmp_ts, v.u(:,:,3), 'u3', settings, name, probeIDs, style);
        Plot_Varts_Data(tmp_dt, tmp_ts, v.T,         'T', settings, name, probeIDs, style);
        Plot_Varts_Data(tmp_dt, tmp_ts, speed,    'umag', settings, name, probeIDs, style);
        Plot_Varts_Data(tmp_dt, tmp_ts, mach,        'M', settings, name, probeIDs, style);
        Plot_Varts_Data(tmp_dt, tmp_ts, rho,       'rho', settings, name, probeIDs, style);
        Plot_Varts_Data(tmp_dt, tmp_ts, v.p,         'p', settings, name, probeIDs, style);
        Plot_Varts_Data(tmp_dt, tmp_ts, dynp,     'dynp', settings, name, probeIDs, style);
        Plot_Varts_Data(tmp_dt, tmp_ts, totp,     'totp', settings, name, probeIDs, style);

        %%%
        % Plot Fourier transforms.
        %%%

        if DO_PLOT_FOURIER_CFD
            for probe_i = 1:length(probeIDs) %#ok<UNRCH>
                fig_number = 100 + probeIDs(probe_i);
                figure(fig_number);
                time_step = varts{dir_i}.dt(1);
                [freq, spectral_power] = Fourier_Transform(v.u(probe_i, :, 2), time_step);
                loglog(freq, spectral_power, 'DisplayName', name);
                title(['Power Spectrum for Probe ', num2str(probeIDs(probe_i))]);
                xlabel('Frequency (Hz)');
                ylabel('Spectral Power');
                hold on;
            end
        end

    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPERIMENTAL Plotting Settings. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experimental data paths. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exp_files = {};
exp_names = {};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exp_files{end+1} = '/users/skinnerr/NGC/Blowing_Profiles_S14.4.2_FromJeremy/P3_mdot1.25_300Hz.txt';
% exp_names{end+1} = 'Exp 300Hz mdot:1.25';
exp_files{end+1} = 'Experimental_Waveforms_2015-08-03/P3_mdot1.25_100Hz_csv.csv';
exp_names{end+1} = 'Exp 100Hz mdot:1.25';
% exp_files{end+1} = 'Experimental_Waveforms_2015-08-03/P3_mdot1.25_300Hz_csv.csv';
% exp_names{end+1} = 'Exp 301.21Hz mdot:1.25';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase-average and offset settings. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
avg_settings_exp = {};
DO_PHASE_AVG_EXP = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experimental_Waveforms_2015-08-03/P3_mdot1.25_100Hz_csv.csv
s.start_ts = 1;
s.end_ts   = 100;
s.freq_Hz  = 100;
avg_settings_exp{end+1} = s;
% % Experimental_Waveforms_2015-08-03/P3_mdot1.25_300Hz_csv.csv
% s.start_ts = 1;
% s.end_ts   = 100000;
% s.freq_Hz  = 301.21;
% avg_settings_exp{end+1} = s;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourier transform flag. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
DO_PLOT_FOURIER_EXP = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(exp_files) < 1
    warning('No experimental data files specified.');
else
    
    %%%
    % Sanitize experimental inputs.
    %%%

    % Directory names.
    if length(exp_files) ~= length(exp_names)
        error('Number of experimental paths is not equal to number of named directories.');
    end

    % Phase-averaging.
    if ~DO_PHASE_AVG_EXP
        avg_settings_exp = {};
    end
    if DO_PHASE_AVG_EXP && (length(avg_settings_exp) ~= length(exp_files))
        error('Phase-average settings are not specified for all experimental directories requested.');
    end
    
    %%%
    % Load and plot experimental data.
    %%%
    
    % Iterate over experimental data files.
    for exp_i = 1:length(exp_files)
        
        name = exp_names{exp_i};
        
        exp_raw = csvread(exp_files{exp_i},1,0);
        exp_n_steps = size(exp_raw,1);
        disp(['There are ', num2str(exp_n_steps), ' time steps ', ...
              'in experimental data ''', name, '''.']);
        
        exp_t    = exp_raw(:,1)';
        exp_umag = exp_raw(:,2)';
        
        if exp_t(1) ~= 0
            warning(['Experimental data does not start at time zero for ', name]);
        end

        % Convert experimental data to a format expected by the plotting code.
        % Assumption: uniform experimental sampling frequency.
        exp_dt = exp_t(2) - exp_t(1);
        exp_ts = 1 + (exp_t ./ exp_dt);

        if DO_PHASE_AVG_EXP
            settings = avg_settings_exp{exp_i};
        else
            settings = {};
        end
        probeIDs = nan;
        style = '-';
        
        Plot_Varts_Data(exp_dt, exp_ts, exp_umag, 'umag', settings, name, probeIDs, style);
    
        %%%
        % Plot Fourier transforms.
        %%%

        if DO_PLOT_FOURIER_EXP
            fig_number = 10000 + exp_i;
            figure(fig_number);
            time_step = exp_dt;
            [freq, spectral_power] = Fourier_Transform(exp_umag, time_step);
            loglog(freq, spectral_power, 'DisplayName', name);
            title(['Power Spectrum for ', name]);
            xlabel('Frequency (Hz)');
            ylabel('Spectral Power');
            hold on;
        end
        
    end

end


















