function [ delta_t, t, probe_locations, probe_data ] = Load_Varts_Data( data_path )
%%%
%
% Loads varts data saved by PHASTA.
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
        delta_t = str2double(tokens(2));
        
        %%%
        % Read subequent lines containing probe IDs and locations.
        % (These lines expected to look like...
        %  '    1       -0.215900000000  0.000000000000  0.000000000000')
        %%%
        
        probe_locations = {};
        
        line = strtrim(fgetl(file_ID));
        while ~strcmp(line,'Probe Data:')
            
            tokens = strsplit(line);
            if length(tokens) ~= 4
                % Ignore lines not in the format of "ID x y z"
                line = strtrim(fgetl(file_ID));
                continue
            end
            x = str2double(tokens(2));
            y = str2double(tokens(3));
            z = str2double(tokens(4));
            probe_locations{end+1} = [x, y, z];
            
            % On to the next line!
            line = strtrim(fgetl(file_ID));
        end
        
        n_probes = length(probe_locations);
        
        %%%
        % Read subequent lines containing probe data at each time step.
        % (These lines expected to look like...
        %  '  281100  0.7908003E+05  0.2385260E+03 -0.5588275E+00 ...')
        %%%
        
        % Each field stores pressure (p), velocity (u,v,w), temperature (T), and kinematic
        % viscosity (nu). This totals six fields per probe point.
        n_fields = 6;
        
        probe_data_raw = {};
        t_start = [];
        
        line = strtrim(fgetl(file_ID));
        while ischar(line)
            
            tokens = strsplit(line);
            tokens(1)
            if length(tokens) ~= 1 + (n_fields * n_probes)
                % Line is truncated or corrupt; time to stop.
                break
            end
            % Save the first timestep.
            if length(t_start) ~= 1
                t_start = str2double(tokens(1));
            end
            % Read all raw numeric data into a buffer.
            blank_indices = strfind(line,' ');
            probe_data_raw{end+1} = str2num(line(blank_indices(1)+1:end));
            
            % On to the next line!
            line = strtrim(fgetl(file_ID));
        end
        
    catch exception
        fclose(file_ID);
        throw(exception);
    end
    fclose(file_ID);

end