function [ xyz, dt, t, p, u, T, nu ] = Load_Varts_Data( data_path )
%%%
%
% Loads varts data saved by PHASTA.
% 
% Ryan Skinner, August 2015
%
%%%

    % Assert that requested_fields is a nonempty string or a nonempty cell array.
    validateattributes(data_path,{'char'},{'nonempty'});

    %%%"
    % Open file and read contents.
    %%%
    
    file_ID = fopen(data_path);
    if file_ID == -1
        error('File does not exist: %s',data_path);
    end
    
    try
        
        %%%
        % Read first line and parse delta_t.
        % (First line expected to look like...
        %  ' Time Step:  0.200000000000000016E-04')
        %%%
        
        first_line = fgetl(file_ID);
        if ~ischar(first_line)
            error('File ended before expected; no time step duration found: %s',data_path);
        end
        tokens = strsplit(first_line,': ');
        dt = str2double(tokens(2));
        
        %%%
        % Move to beginning of probe location data.
        %%%
        
        line = fgetl(file_ID);
        while ~any(strfind(line,'Probe ID'))
            line = strtrim(fgetl(file_ID));
        end
        
        %%%
        % Read subequent lines containing probe IDs and locations.
        % (These lines expected to look like...
        %  '    1       -0.215900000000  0.000000000000  0.000000000000')
        %%%
        
        xyz = fscanf(file_ID,'%f');
        
        % Reshape so we're indexed by (probe, field).
        xyz = reshape(xyz, length(xyz)/4, 4);
        
        % Remove the first field, which is just the probe ID.
        xyz = xyz(:,2:4);
        
        n_probes = length(xyz);
        
        %%%
        % Move to beginning of probe trace data.
        %%%
        
        line = fgetl(file_ID);
        while ~any(strfind(line,'Probe Data'))
            line = strtrim(fgetl(file_ID));
        end
        
        %%%
        % Read subequent lines containing probe data at each time step.
        % (These lines expected to look like...
        %  '  281100  0.7908003E+05  0.2385260E+03 -0.5588275E+00 ...'
        %  i.e. timestep, fields for point 1, fields for point 2, etc.)
        %%%
        
        % Each field stores pressure (p), velocity (u,v,w), temperature (T), and kinematic
        % viscosity (nu). This totals six fields per probe point.
        n_fields = 6;
        
        n_vars = n_fields * n_probes;
        
        probe_data = fscanf(file_ID,'%f');
        
        % Ensure data is not truncated or corrupt.
        if mod(length(probe_data), 1 + (n_fields * n_probes)) ~= 0
            error('Probe data is corrupt for %s.', data_path);
        end
        
        increment = 1 + n_vars;
        t = probe_data(1:increment:end);
        probe_data(1:increment:end) = [];
        
        probe_data = reshape(probe_data, n_vars, length(probe_data) / n_vars)';
        
        % These are now indexed by field(probeID, timestep).
        p        = probe_data(:,1:n_fields:end)';
        u(:,:,1) = probe_data(:,2:n_fields:end)';
        u(:,:,2) = probe_data(:,3:n_fields:end)';
        u(:,:,3) = probe_data(:,4:n_fields:end)';
        T        = probe_data(:,5:n_fields:end)';
        nu       = probe_data(:,6:n_fields:end)';
        
    catch exception
        fclose(file_ID);
        throw(exception);
    end
    fclose(file_ID);

end