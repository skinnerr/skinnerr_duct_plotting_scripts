function [ x, y, z, p, cp, pr, ref ] = Load_CFD_Data( data_path, ref_path )
%%%
% Loads data from a CFD run saved in ParaView, and calculates Cp, PR, etc.
%
% Usage:
%   Load_CFD_Data('particular_slice.csv','ref.csv')
%%%

    % Assert that inputs are nonempty strings.
    validateattributes(data_path,{'char'},{'nonempty'});
    validateattributes(ref_path,{'char'},{'nonempty'});
    
    %%%
    % Load raw data from ParaView files.
    %%%
    
    [p, T, ...
     u, v, w, ...
     x, y, z] = ...
        Load_ParaView_Data(data_path, {'p', 'T', 'u:0', 'u:1', 'u:2', ...
                                      'Points:0', 'Points:1', 'Points:2'});
    
    [ref.p, ref.T, ...
     ref.u(1), ref.u(2), ref.u(3), ...
     ref.x, ref.y, ref.z] = ...
        Load_ParaView_Data(ref_path, {'p', 'T', ...
                                      'u:0', 'u:1', 'u:2', ...
                                      'Points:0', 'Points:1', 'Points:2'});
    
	%%%
    % Compute extra fields relevant to flow analysis.
    %%%
    
    % Physical constants.
	gamma = 1.4;
    R = 288.294;
    
    % Reference probe Mach number.
    ref.M = norm(ref.u) / sqrt(gamma * R * ref.T);
    
    % Flow pressure coefficient.
    cp = 2*(p/ref.p - 1) / (gamma * ref.M^2);
    
    % Flow pressure recovery.
    uMag = sqrt(u.^2 + v.^2 + w.^2);
    numerator = p .* (1 + 0.5 * uMag.^2 ./ (R * T));
    denomenator = ref.p * (1 + 0.5 * norm(ref.u)^2 / (R * ref.T));
    pr = numerator / denomenator;

end