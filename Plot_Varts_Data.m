function [ output_args ] = Plot_Varts_Data( dt, ts, field, field_name, avg_settings, ...
                                            name, probeIDs, style )
%%%
%
% Plots varts data.
%
% Ryan Skinner, August 2015
%
%%%

    %%%
    % Validate inputs.
    %%%
    
    if length(probeIDs) ~= size(field,1)
        error('First dimension of field data must match number of probes.')
    end

    %%%
    % Set up plot size and initialize figure/axes.
    %%%
    
    fig_width  = 1000;
    fig_height = 300;
    
    valid_plot_fields = {'p','u1','u2','u3','T','magu','rho','totp','dynp','M','nu'};
    
    field_index = find(strcmp(valid_plot_fields, field_name));
    
    if length(field_index) < 1
        error(['''%s'' is not a valid field name. ', ...
               'Please add it here and create axis labels.'], field_name);
    end
    if length(field_index) > 1
        error('''%s'' is duplicated in the list of valid fields. Please fix.', ...
              field_name);
    end
    
    fig_num_search = find(strcmp(field_name, valid_plot_fields));
    
    if ~any(fig_num_search)
        %%%
        % If this is our first time accessing the plot...
        %%%
        
        h = figure('Position', Centered_Figure_Position(fig_width,fig_height));
        hax1 = axes();
        c_order = get(hax1,'ColorOrder');
        set(hax1,'ColorOrder',c_order([5,7,1,2,3,4,6],:));
    
    else
        %%%
        % If we have already plotted here...
        %%%
        
        h = figure(fig_num_search);
        set(h, 'Position', Centered_Figure_Position(fig_width,fig_height));
        hax1 = gca();
    end
    
    %%%
    % Perform phase-averaging if requested.
    %%%
    
    if ~isempty(avg_settings)
        if isnan(dt)
            error('Phase averaging requires constant time step size.');
        end
        avg_field = [];
        for probe_i = 1:size(field,1)
            [avg_ts, tmp_field] = Phase_Average_Varts(dt, ts, field(probe_i,:), ...
                                                      avg_settings);
            avg_field = cat(1, avg_field, tmp_field);
        end
        ts    = avg_ts;
        field = avg_field;
    end
    
    %%%
    % Plot the data!
    %%%

    hold on;
    for probe_i = 1:size(field,1)
        plot(ts, field(probe_i,:), ...
             'LineStyle', style, ...
             'DisplayName', [name, ' ID:', num2str(probeIDs(probe_i))]);
    end
    
    %%%
    % Add annotations: legend and labels.
    %%%
    
    legend();
    
    y_labels = {'Pressure (pa)', ...
                'x-Velocity (m/s)', ...
                'y-Velocity (m/s)', ...
                'z-Velocity (m/s)', ...
                'Temperature (K)', ...
                'Speed (m/s)', ...
                'Density (kg/m3)', ...
                'Total Pressure (Pa)', ...
                'Dynamic Pressure (Pa)', ...
                'Mach', ...
                'Kinematic Viscosity'};
	ylabel(y_labels{field_index});
    xlabel('Time Step');
    
    % Reset color plotting order for the next time the plot is drawn to.
    set(hax1,'ColorOrderIndex',1);
    
end

