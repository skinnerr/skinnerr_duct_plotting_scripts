function [ dt, t, xyz, p, u, T, nu ] = Load_Varts_Directory( directory_path )
%%%
%
% Loads all varts data files from a directory.
% 
% Ryan Skinner, August 2015
%
%%%

    % Flag to insert NaNs into data to bridge non-contiguous time steps.
    insertNaNs = true;

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
    
    dt  = [];
    t   = [];
    xyz = [];
    p   = [];
    u   = [];
    T   = [];
    nu  = [];
    
    load_str = 'Loading varts files...';
    hwait = waitbar(0,sprintf('%s %i / %i',load_str,1,length(files)), ...
                    'CreateCancelBtn','setappdata(gcbf,''cancel_loading'',1)');
    setappdata(hwait,'cancel_loading',0);
    for i = 1:length(files)
        if getappdata(hwait,'cancel_loading')
            delete(hwait);
            error('Directory loading canceled by user. Program stopped.');
        end
        waitbar(i/length(files),hwait,sprintf('%s %i / %i',load_str,i,length(files)));
        data_path = [directory_path,files(i).name];
        
        [tmp_xyz, tmp_dt, tmp_t, tmp_p, tmp_u, tmp_T, tmp_nu] = ...
            Load_Varts_Data(data_path);
        
        if i > 1
            % Warn if probe locations are not consistent.
            if any(size(tmp_xyz) ~= size(xyz))
                warning('Probe locations change in file %s', data_path);
            end
            % Warn if time step size is not consistent.
            if tmp_dt ~= dt(end)
                warning('Time step size change in file %s', data_path);
            end
            % Insert NaNs if we are missing some time steps.
            if insertNaNs && (t(end) + 1 ~= tmp_t(1))
                sprintf('Some time steps skipped before file %s', data_path);
                t  = cat(1, t, nan);
                p  = cat(2, p,  NaN(size( p,1), 1));
                u  = cat(2, u,  NaN(size( u,1), 1, 3));
                T  = cat(2, T,  NaN(size( T,1), 1));
                nu = cat(2, nu, NaN(size(nu,1), 1));
                
            end
        end
        
        xyz = tmp_xyz;
        dt  = cat(1,  dt, tmp_dt );
        t   = cat(1,   t, tmp_t  );
        p   = cat(2,   p, tmp_p  );
        u   = cat(2,   u, tmp_u  );
        T   = cat(2,   T, tmp_T  );
        nu  = cat(2,  nu, tmp_nu );
        
    end
    delete(hwait);

end