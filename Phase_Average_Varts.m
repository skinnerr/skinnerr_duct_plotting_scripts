function [ avg_t, avg_field, max_count ] = Phase_Average_Varts( dt, ts, field, settings )
%%%
%
% Performs phase averaging on a given set of data.
%
% Ryan Skinner, August 2015
%
%%%

if settings.start_ts >= settings.end_ts
    error('Out-of-order of start and end time steps for phase averaging.');
end

%%%
% Create bins for phase averaging.
%%%

% Calculate number of bins.
period = 1 / settings.freq_Hz;
ts_per_period = period / dt;
n_bins = ceil(ts_per_period);

% Create phase-averaged time bins, field bins, and counter bins.
avg_t = dt*((1:n_bins)-1);
avg_field = zeros(1,n_bins);
avg_count = zeros(1,n_bins);

% Ensure window contains multiple periods.
if settings.end_ts - settings.start_ts <= n_bins
    warning('Phase averaging requested time window is too small.');
    avg_t = ts;
    avg_field = field;
    return
end

%%%
% Determine where to start and stop averaging.
%%%

start_i = find(ts == settings.start_ts, 1);
if length(start_i) < 1
    error('Start time step not found in data for phase averaging.');
end

end_i = find(ts == settings.end_ts, 1);
if length(end_i) < 1
    error('End time step not found in data for phase averaging.');
end

%%%
% Bin and calculate phase averages.
%%%

norm_ts = ts(start_i:end_i) - settings.start_ts;
for nts_i = 1:length(norm_ts);
    if isnan(norm_ts(nts_i))
        continue
    end
    bin_i = 1 + floor(mod(norm_ts(nts_i), ts_per_period));
    avg_field(bin_i) = avg_field(bin_i) + field(nts_i - 1 + start_i);
    avg_count(bin_i) = avg_count(bin_i) + 1;
end
avg_field = avg_field ./ avg_count;

max_count = max(avg_count);


end

