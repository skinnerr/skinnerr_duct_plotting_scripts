function [ dir_starts, dt, ts, xyz, p, u, T, nu ] = ...
    Load_Varts_Directory( directory_path, do_offset_data, offset )
%%%
%
% Loads all varts data files from a directory.
% 
% Ryan Skinner, August 2015
%
%%%

    % Flag to insert NaNs into data to bridge non-contiguous time steps.
    insert_nans = true;

    %%%
    % Sanitize inputs.
    %%%
    
    validateattributes(directory_path,{'char'},{'nonempty'});
    
    %%%
    % Obtain a list of all 'varts.*.dat' files in the directory.
    %%%
    
    % List all files.
    files = dir(directory_path);
    % Retain only files that are not directories.
    files = files(find( ~[files.isdir] )); %#ok<*FNDSB>
    % Retain only files that match our pattern.
    files = files(find( ~cellfun(@isempty, regexp({files.name},'^varts.*.dat$')) ));
    
    if isempty(files)
        error('No varts files found in direcotry %s', directory_path);
    end
    
    if do_offset_data
        file_i_keep = [];
        for i = 1:(length(files)-1)
            tokens = strsplit(files(i+1).name,'.');
            next_varts_ts = str2double(tokens(2)) + 1; % This +1 is important, because the
                                                       %  varts.x.dat file begins with the
                                                       %  timestep x+1, rather than x.
            if next_varts_ts <= offset
                disp(['Skipping ', files(i).name, ' because data is before the offset.']);
            else
                file_i_keep = [file_i_keep, i];
            end
        end
        files = files(file_i_keep);
    end
    
    %%%
    % Sort files by starting time step indicated in filename.
    %%%
    
    file_ts = zeros(1,length(files));
    for i = 1:length(files)
        tokens = strsplit(files(i).name,'.');
        file_ts(i) = str2double(tokens(2));
    end
    [~,sort_order] = sort(file_ts);
    files = files(sort_order);
    
    %%%
    % Load varts files and build the output data structure.
    %%%
    
    dir_starts = [];
    dt         = [];
    ts         = [];
    xyz        = [];
    p          = [];
    u          = [];
    T          = [];
    nu         = [];
    
    load_str = 'Loading varts files...';
    hwait = waitbar(0, sprintf('%s %i / %i', load_str, 1, length(files)), ...
                    'CreateCancelBtn', 'setappdata(gcbf, ''cancel_loading'', 1)');
    setappdata(hwait, 'cancel_loading',0);
    
    for i = 1:length(files)
        if getappdata(hwait, 'cancel_loading')
            delete(hwait);
            error('Directory loading canceled by user. Program stopped.');
        end
        waitbar(i/length(files), hwait, ...
                sprintf('%s %i / %i', load_str, i, length(files)));
        data_path = [directory_path, files(i).name];
        
        try
            [tmp_xyz, tmp_dt, tmp_ts, tmp_p, tmp_u, tmp_T, tmp_nu] = ...
                Load_Varts_Data(data_path);
            
            inserted_nans = false;

            if i > 1
                % Warn if time step size is not consistent.
                if tmp_dt ~= dt(end)
                    disp(['Warning: Time step size change', ... 
                          ' from ', num2str(dt(end)), ...
                          ' to ',   num2str(tmp_dt), ...
                          ' in file ', data_path]);
                end
                % Warn if probe locations are not consistent, and pad with NaNs.
                n_xyz_t = size(tmp_xyz, 1);
                n_xyz   = size(xyz, 1);
                if n_xyz_t ~= n_xyz
                    disp(['Warning: Probe locations change in file ', data_path]);
                end
                % Insert NaNs if we are missing some time steps.
                if insert_nans && (ts(end) + 1 ~= tmp_ts(1))
                    disp(['Some time steps skipped before file ', data_path]);
                    ts = cat(1, ts, nan);
                    p  = cat(2,  p, NaN(size( p,1), 1));
                    u  = cat(2,  u, NaN(size( u,1), 1, 3));
                    T  = cat(2,  T, NaN(size( T,1), 1));
                    nu = cat(2, nu, NaN(size(nu,1), 1));
                    inserted_nans = true;
                end
            end
            
            % First index of each file, including the preceeding NaN, if one exists.
            dir_starts(end+1) = 1 + length(ts) - inserted_nans; %#ok<AGROW>

            dt  = cat(1,  dt, tmp_dt );
            ts  = cat(1,  ts, tmp_ts );
            p   = cat(2,   p, tmp_p  );
            u   = cat(2,   u, tmp_u  );
            T   = cat(2,   T, tmp_T  );
            nu  = cat(2,  nu, tmp_nu );
            xyz = tmp_xyz; % See above warning about changed probe xyz.
            
        catch err
            disp(['Skipping malformed file ', data_path]);
        end
        
    end
    delete(hwait);

end
