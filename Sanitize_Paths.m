%%%%%%
%
% Sanitizes a list of paths by:
%  - Checking for existence.
%  - Appending forward slashes as needed.
%
%%%

function [ dir_paths ] = Sanitize_Paths( dir_paths )

if length(dir_paths) < 1
    error('No directories specified.');
    return
end

for i = 1:length(dir_paths)
    % Append forward-slash to all paths names if needed.
    if dir_paths{i}(end) ~= '/'
        dir_paths{i}(end+1) = '/';
    end
    if ~isdir(dir_paths{i})
        error('Directory does not exist: %s',dir_paths{i});
    end
end
