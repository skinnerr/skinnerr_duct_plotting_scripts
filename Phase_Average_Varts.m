function [ avg_varts_data ] = Phase_Average_Varts( varts_data, settings )
%%%
%
% Performs phase averaging on a given structure of varts data.
%
% Ryan Skinner, August 2015
%
%%%

if settings.start_ts >= settings.end_ts
    error('Out-of-order of start and end time steps for phase averaging.');
end

% Assume uniform time step size.
if varts_data.dt(1) ~= mean(varts_data.dt)
    error('Phase averaging for variable-timestep varts data not implemented.');
end
dt = varts_data.dt(1);

%%%
% Create bins for phase averaging.
%%%

% Calculate number of bins.
period = 1 / settings.freq_Hz;
n_bins = floor(period / dt);

% Create bin times.
t_bins = dt*((1:n_bins)-1);

% Ensure window contains multiple periods.
if settings.end_ts - settings.start_ts <= n_bins
    warning('Phase averaging requested time window is too small.');
    avg_varts_data = varts_data;
    return
end

%%%
% Determine where to start and stop averaging.
%%%

start_i = find(varts_data.t == settings.start_ts, 1);
if length(start_i) < 1
    error('Start time step not found in data for phase averaging.');
end

end_i = find(varts_data.t == settings.end_ts, 1);
if length(end_i) < 1
    error('Start time step not found in data for phase averaging.');
end

%%%
% Bin and calculate phase averages.
%%%

ts = varts_data.t(start_i:end_i);
for ts_i = 1:length(ts);
    bin_i = mod((ts(ts_i) - settings.start_ts), n_bins);
end


end

