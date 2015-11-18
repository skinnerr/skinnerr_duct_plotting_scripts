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

min_time = nan;
max_time = nan;

% Physical constants.
gamma = 1.4;
R     = 288.294;

% Line styles to cycle through for each directory.
styles = {'-','-',':','-.','--'};

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
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-DDES-100Hz-ExpMatch150803-LBOff/'; dir_names{end+1} = 'S14.4.2 A0-Plug 100Hz';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-DDES-100Hz-ExpMatch150803-LBOff-P110k/'; dir_names{end+1} = 'S14.4.2 A0-Plug 100Hz P110k';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-DDES-100Hz-ExpMatch150803-LBOff-P100k/'; dir_names{end+1} = 'S14.4.2 A0-Plug 100Hz P100k';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-DDES-100Hz-ExpMatch150803-LBOff-P101.3k/'; dir_names{end+1} = 'S14.4.2 A0-Plug 100Hz P101.3k';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-DDES-Steady-LBOff-MassFlowMatch/'; dir_names{end+1} = 'S14.4.2 A0-Plug Steady MassFlowMatch';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-DDES-Steady-mdot1.25/'; dir_names{end+1} = 'x/L = 0, mdot 1.25%';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-DDES-Steady-mdot1.5/'; dir_names{end+1} = 'x/L = 0, mdot 1.5%';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-RANS-Steady-mdot1.25/'; dir_names{end+1} = 'x/L = 0, mdot 1.25%';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-RANS-Steady-mdot1.5/'; dir_names{end+1} = 'x/L = 0, mdot 1.5%';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-RANS-Steady-mdot1.5-try2/'; dir_names{end+1} = 'x/L = 0, mdot 1.5% (try2)';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-RANS-Steady-mdot2.0/'; dir_names{end+1} = 'x/L = 0, mdot 2.0%';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-RANS-Steady-mdot2.5/'; dir_names{end+1} = 'x/L = 0, mdot 2.5%';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/Run-RANS-Steady-mdot3.0/'; dir_names{end+1} = 'x/L = 0, mdot 3.0%';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-DDES-Steady-LBOff';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-DDES-Baseline';
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-DDES-Baseline25k-recent'; dir_names{end+1} = 'S14.4.2 A1-Iso Baseline';
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
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par/A1-Iso/16k-procs/Run-DDES-100Hz-LowBlowOff'; dir_names{end+1} = 'S14.4.2 A1-Iso DDES 100Hz LBOff';
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
% dirs{end+1} = '../series14.4.2/meshing-Coarse2-Par-waggingjet-fullUBmesh-Coarse2-fixedBL/Run-DDES-Wag-100Hz-ExpMatch150803-LowBlowOff'; dir_names{end+1} = 'S14.4.2 A0-Wag-Fine 100Hz';
% 
% dirs{end+1} = '../series14.5.1/A0/A0-Run-RANS-Start2';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-DDES-Base-KJ';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-DDES-Steady-LBOff'; dir_names{end+1} = 'x/L = 0.25, mdot 1.25% (DDES)';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-DDES-Steady-LBOff-UR';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-DDES-Steady-mdot1.25'; dir_names{end+1} = 'x/L = 0.25, mdot 1.25%';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-DDES-Steady-mdot1.5'; dir_names{end+1} = 'x/L = 0.25, mdot 1.8%';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-DDES-Steady-mdot1.25-LB1.72'; dir_names{end+1} = 'x/L = 0.25, mdot 1.25%';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-DDES-Steady-mdot1.5-LB1.72'; dir_names{end+1} = 'x/L = 0.25, mdot 1.8%';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-DDES-100Hz-LB0.86-KJ';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-DDES-100Hz-LB1.72-KJ';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-DDES-100Hz-LB1.72-KJ-UR';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-DDES-100Hz-LB1.72-KJ-NoNormalUR';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-DDES-100Hz-KJ';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-RANS-100Hz-LB1.72';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-RANS-100Hz-LB1.72-UR';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-RANS-Steady-mdot1.25'; dir_names{end+1} = 'x/L = 0.25, mdot 1.25%';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-RANS-Steady-mdot1.5'; dir_names{end+1} = 'x/L = 0.25, mdot 1.8%';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-RANS-Steady-mdot2.0'; dir_names{end+1} = 'x/L = 0.25, mdot 2.0%';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-RANS-Steady-mdot2.5'; dir_names{end+1} = 'x/L = 0.25, mdot 2.5%';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-RANS-Steady-mdot3.0'; dir_names{end+1} = 'x/L = 0.25, mdot 3.0%';
% dirs{end+1} = '../series14.5.1/A0/A0-Run-RANS-Steady-mdot1.25-LB1.72'; dir_names{end+1} = 'x/L = 0.25, mdot 1.25%';
dirs{end+1} = '../series14.5.1/A0/A0-Run-RANS-Steady-mdot1.5-LB1.72'; dir_names{end+1} = 'x/L = 0.25, mdot 1.8%';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dirs{end+1} = 'varts_1'; dir_names{end+1} = 'varts_1';
% dirs{end+1} = 'varts_2'; dir_names{end+1} = 'varts_2';
% dirs{end+1} = 'varts_3'; dir_names{end+1} = 'varts_3';
% dirs{end+1} = 'varts_4'; dir_names{end+1} = 'varts_4';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Probe point selection. %
%%%%%%%%%%%%%%%%%%%%%%%%%%
probeIDs = [1, 261, 264, 267, 270, 273];
% probeIDs = [11, 10, 9, 8];
% probeIDs = [1, 261, 264, 267]; % Upper blower
% probeIDs = [1, 291, 294, 297, 300]; % Lower blower
% probeIDs = [1, 261, 264, 267, 270, 291, 294, 297, 300]; % Upper AND Lower blowers
% probeIDs = [38:42, 53:57]; % AIP z = -2 in
% probeIDs = [43:47, 58:62]; % AIP z =  0 in
% probeIDs = [48:52, 63:67]; % AIP z = +2 in
% probeIDs = [267];
% probeIDs = [302, 303, 267, 304, 305]; % Five span-wise UB throat probes (horizontal).
% probeIDs = [306, 307, 267, 308, 309]; % Five span-wise UB throat probes (vertical).
%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase-average settings. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
avg_settings_cfd = {};
DO_PHASE_AVG_CFD = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % A1-Iso trapezoidal expmatch
% s.start_ts = 62000;
% s.end_ts   = 68800;
% s.freq_Hz  = 300;
% avg_settings_cfd{end+1} = s;

% % Wag mesh coarse2 expmatch
% s.start_ts = 18135;
% s.end_ts   = 42500;
% s.freq_Hz  = 300;
% avg_settings_cfd{end+1} = s;

% % Wag mesh coarse4 expmatch
% s.start_ts = 13086;
% s.end_ts   = 35100;
% s.freq_Hz  = 300;
% avg_settings_cfd{end+1} = s;

% % S14.4.2 A0 ExpMatch150803 100Hz
% s.start_ts = 363675;
% s.end_ts   = 375600;
% s.freq_Hz  = 100;
% avg_settings_cfd{end+1} = s;

% S14.4.2 A0 ExpMatch150803 100Hz - mdot 1.18%, ts 426.1k-460.1k
s.start_ts = 1573;
s.end_ts   = 32000;
s.freq_Hz  = 100;
avg_settings_cfd{end+1} = s;

% % S14.4.2 A0 meshing-Coarse2-Par-waggingjet-fullUBmesh-Coarse2-fixedBL ExpMatch150803 100Hz-Wag
% s.start_ts = 38261;
% s.end_ts   = 55000;
% s.freq_Hz  = 100;
% avg_settings_cfd{end+1} = s;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Test directory
% s.start_ts = 61000;
% s.end_ts   = 65000;
% s.freq_Hz  = 300;
% avg_settings_cfd{end+1} = s;
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Offset settings. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
offsets_cfd = {};
DO_OFFSET_CFD = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
offsets_cfd{end+1} = 390000;
offsets_cfd{end+1} = 426101;
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Whether to convert CFD to time (sec). %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DO_PLOT_TIME = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourier transform flag. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
DO_PLOT_FOURIER_CFD = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%
% Sanitize CFD inputs.
%%%

if length(dirs) < 1
    disp('No CFD directories specified.');
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
    
    % Offsets.
    if DO_OFFSET_CFD
        if (length(offsets_cfd) ~= length(dirs))
            error('Offsets are not specified for all CFD directories requested.');
        end
    else
        offsets_cfd = cell(1,length(dirs));
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

        [dir_starts, dt, ts, xyz, p, u, T, nu] = ...
            Load_Varts_Directory(dirs{dir_i}, DO_OFFSET_CFD, offsets_cfd{dir_i});

        % Retain only requested probe points.
        p  =  p(probeIDs,:);
        u  =  u(probeIDs,:,:);
        T  =  T(probeIDs,:);
        nu = nu(probeIDs,:);
        
        % Remove all data before the offset.
        if DO_OFFSET_CFD
            ts_offset = offsets_cfd{dir_i};
            if ts(end) < ts_offset
                error('Offset is greater than maximum time step in %s', dirs{dir_i});
            else
                % Deal with cropping; need to modify a few data structures.
                noncropped_dirs_i = (ts(dir_starts) - ts_offset) > 0;
                first_cropped_i = find(noncropped_dirs_i);
                first_cropped_i = first_cropped_i(1) - 1;
                if first_cropped_i > 0
                    % Accout for the case where the start time step of a
                    % file dictates we chop it off, but it still contains
                    % relevant data.
                    noncropped_dirs_i(first_cropped_i) = 1;
                    ts_chopped_in_first_cropped_i = ...
                        ts_offset - ts(dir_starts(first_cropped_i));
                end
                dir_starts = dir_starts(noncropped_dirs_i);
                dir_starts = dir_starts + 1 - dir_starts(1);
                dir_starts(2:end) = dir_starts(2:end) - ts_chopped_in_first_cropped_i;
                dt = dt(noncropped_dirs_i);
                ts = ts - ts_offset;
                noncropped_ts_i = ts > 0;
                p  = p(:,noncropped_ts_i);
                u  = u(:,noncropped_ts_i);
                T  = T(:,noncropped_ts_i);
                nu = nu(:,noncropped_ts_i);
                ts_crop_i = find(ts > 0, 1);
                ts = ts(ts_crop_i:end);
            end
        end

        % Save data to structure
        varts{dir_i}.dir_starts = dir_starts;
        varts{dir_i}.dt         = dt;
        varts{dir_i}.ts         = ts;
        varts{dir_i}.xyz        = xyz;
        varts{dir_i}.p          = p;
        varts{dir_i}.u          = u;
        varts{dir_i}.T          = T;
        varts{dir_i}.nu         = nu;

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
        v      = varts{dir_i};
        ds     = v.dir_starts;
        dt     = v.dt;
        ts     = v.ts;
        speed  = sqrt(sum(v.u(:, :, :) .* v.u(:, :, :), 3));
        mach   = sqrt(sum(v.u(:, :, :) .* v.u(:, :, :), 3) ./ (gamma*R*v.T(:, :)));
        rho    = v.p(:,:,:) ./ (R*v.T(:, :));
        dynp   = 0.5 * rho .* speed.^2;
        totp   = v.p(:,:) + dynp;
        
        if DO_PLOT_TIME
            % Convert time steps to physical times, skipping NaN values,
            % which indicate that some time steps were skipped during load.
            time = nan(length(ts),1);
            n_files = length(ds);
            for file_i = 1:n_files
                start_i = ds(file_i);
                if file_i < n_files
                    end_i = ds(file_i + 1) - 1;
                else
                    end_i = length(ts);
                end
                if isnan(ts(start_i))
                    % Skip this index so it remains a NaN.
                    start_i = start_i + 1;
                end
                current_time = 0;
                if start_i > 1
                    current_time = time(start_i - 1);
                end
                time(start_i:end_i) = current_time + dt(file_i)*(1:(1+end_i-start_i));
            end
        else
            time = ts;
        end
        
        if min(time) < min_time || isnan(min_time)
            min_time = min(time);
        end
        if max(time) > max_time || isnan(max_time)
            max_time = max(time);
        end

        %%%
        % Plot time histories.
        %%%

%         Plot_Varts_Data(dt, time, min_time, max_time, v.u(:,:,1), 'u1', settings, name, probeIDs, style, DO_PLOT_TIME);
%         Plot_Varts_Data(dt, time, min_time, max_time, v.u(:,:,2), 'u2', settings, name, probeIDs, style, DO_PLOT_TIME);
%         Plot_Varts_Data(dt, time, min_time, max_time, v.u(:,:,3), 'u3', settings, name, probeIDs, style, DO_PLOT_TIME);
        Plot_Varts_Data(dt, time, min_time, max_time, v.T,         'T', settings, name, probeIDs, style, DO_PLOT_TIME);
        Plot_Varts_Data(dt, time, min_time, max_time, speed,    'umag', settings, name, probeIDs, style, DO_PLOT_TIME);
        Plot_Varts_Data(dt, time, min_time, max_time, mach,        'M', settings, name, probeIDs, style, DO_PLOT_TIME);
        Plot_Varts_Data(dt, time, min_time, max_time, rho,       'rho', settings, name, probeIDs, style, DO_PLOT_TIME);
        Plot_Varts_Data(dt, time, min_time, max_time, v.p,         'p', settings, name, probeIDs, style, DO_PLOT_TIME);
        Plot_Varts_Data(dt, time, min_time, max_time, dynp,     'dynp', settings, name, probeIDs, style, DO_PLOT_TIME);
        Plot_Varts_Data(dt, time, min_time, max_time, totp,     'totp', settings, name, probeIDs, style, DO_PLOT_TIME);
%         mean(totp')'
        
        % Plot min and max of experiment that we're trying to achieve.
        target_min = 119;
        target_max = 247;
        domain = [min_time, max_time];
%         Plot_Varts_Data(1, domain, min_time, max_time, target_min*[1,1], 'umag', {}, 'Target Min', 267, '--', DO_PLOT_TIME);
%         Plot_Varts_Data(1, domain, min_time, max_time, target_max*[1,1], 'umag', {}, 'Target Max', 267, '--', DO_PLOT_TIME);

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
% exp_files{end+1} = 'Experimental_Waveforms_2015-08-03/P3_mdot1.25_100Hz_csv.csv';
% exp_names{end+1} = 'Exp 100.76Hz mdot:1.25';
% exp_files{end+1} = 'Experimental_Waveforms_2015-08-03/P3_mdot1.25_300Hz_csv.csv';
% exp_names{end+1} = 'Exp 301.21Hz mdot:1.25';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase-average and repeat settings. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
avg_settings_exp = {};
DO_PHASE_AVG_EXP = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experimental_Waveforms_2015-08-03/P3_mdot1.25_100Hz_csv.csv
s.start_ts = 89;
s.end_ts   = 100000;
s.freq_Hz  = 101.76;
avg_settings_exp{end+1} = s;
% % Experimental_Waveforms_2015-08-03/P3_mdot1.25_300Hz_csv.csv
% s.start_ts = 25;
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
    disp('No experimental data files specified.');
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
        exp_ts = 1 + round(exp_t ./ exp_dt);
        % IMPORTANT: round() prevents ts from containing elements like 2.0000000e+5.

        if DO_PHASE_AVG_EXP
            settings = avg_settings_exp{exp_i};
        else
            settings = {};
        end
        probeIDs = nan;
        style = '-';
        
        if min(exp_ts) < min_time
            min_time = min(exp_ts);
        end
        if max(exp_ts) > max_time
            max_time = max(exp_ts);
        end
        
        Plot_Varts_Data(exp_dt, exp_ts, min_time, max_time, exp_umag, 'umag', settings, name, probeIDs, style, DO_PLOT_TIME);
    
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
