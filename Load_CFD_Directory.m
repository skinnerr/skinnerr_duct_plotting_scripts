function [ CFD_data ] = Load_CFD_Directory( varargin )
%%%
% Loads all the .csv files from a directory. Requries a 'ref.csv' file to work:
%   either specify the ref_file_path, or it will look for the file within
%   'directory_path/..'.
%
% Input: Load_CFD_Data(directory_path)
%        Load_CFD_Data(directory_path, filter_string)
%         (where this function ignores all filenames not containing filter_string)
%
% Output: CFD_Data is an array containing a struct for each file loaded in.
%%%

    %%%
    % Sanitize inputs.
    %%%

    if nargin < 0 || nargin > 2
        error('Invalid number of arguments passed in.');
    end
    
    directory_path = varargin{1};
    if directory_path(end) ~= '/'
        directory_path(end+1) = '/';
    end
    
    if nargin < 2
        file_filter_string = '';
    else
        file_filter_string = varargin{2};
    end

    % Assert that inputs are strings.
    validateattributes(directory_path,    {'char'},{'nonempty'});
    validateattributes(file_filter_string,{'char'},{});
    
    %%%
    % Obtain a list of all '.csv' files in the directory.
    %%%
    
    % List all files.
    files = dir(directory_path);
    % Retain only files that are not directories.
    files = files(find( ~[files.isdir] )); %#ok<*FNDSB>
    % Retain only files that end in '.csv'.
    files = files(find( ~cellfun(@isempty, regexp({files.name},'.csv$')) ));
    % Retain only files that contain the filter string, if it is non-empty.
    if ~isempty(file_filter_string)
        files = files(find( ...
                    ~cellfun(@isempty, regexp({files.name},file_filter_string)) ));
    end
    
    %%%
    % Find ref file in the directory. If it's already given by
    %   ref_file_path, remove it from our list of files.
    %%%
    
    % Find the reference probe file in this directory, and remove it from our list.
    ref_file_name = 'ref.csv';
    ref_file_index = find(strcmp({files.name},ref_file_name), 1);
    files(ref_file_index) = [];
    
    %%%
    % Load the remaining '.csv' files and build the output data structure.
    %%%
    
    load_str = 'Loading directory files...';
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
        ref_path = [directory_path,'ref.csv'];
        [CFD_data(i).x, ...
         CFD_data(i).y, ...
         CFD_data(i).z, ...
         CFD_data(i).p, ...
         CFD_data(i).cp, ...
         CFD_data(i).pr, ...
         CFD_data(i).ref] = Load_CFD_Data(data_path,ref_path); %#ok<*AGROW>
        CFD_data(i).filename = files(i).name;
    end
    delete(hwait);

end