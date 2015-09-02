function [ output_args ] = Plot_Varts_Data( dt, ts, field, field_name, avg_settings, ...
                                            name, probeIDs, style, label_seconds )
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
    
    plotting_exp = isnan(probeIDs(1));
    do_phase_avg = ~isempty(avg_settings);
    
    if ~plotting_exp && (length(probeIDs) ~= size(field,1))
        error('First dimension of field data must match number of probes.')
    end

    %%%
    % Set up plot size and initialize figure/axes.
    %%%
    
    fig_width  = 1000;
    fig_height = 300;
    
    valid_plot_fields = {'u1','u2','u3','T','umag','M','rho','p','totp','dynp','nu'};
    
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
    
    if do_phase_avg
        if mean(dt) ~= dt(1)
            error('Phase averaging requires constant time step size.');
        end
        avg_field = [];
        for probe_i = 1:size(field,1)
            [avg_t, tmp_field, n_periods] = ...
                Phase_Average_Varts(dt(1), ts, field(probe_i,:), ...
                                    avg_settings);
            avg_field = cat(1, avg_field, tmp_field);
        end
        time  = avg_t;
        field = avg_field;
        x_label = 'Time (sec)';
    else
        time  = ts;
        x_label = 'Time Step';
        if label_seconds
            x_label = 'Time (sec)';
        end
    end
    
    %%%
    % Plot the data!
    %%%
    
    display_append = '';
    if ~isempty(avg_settings)
        display_append = ['Periods:', sprintf('%.1f',n_periods)];
    end
    
    hold on;
    if plotting_exp
        display_name = [name, ' ', display_append];
        plot(time, field, 'LineStyle', style, 'DisplayName', display_name);
    else
        for probe_i = 1:size(field,1)
            display_name = [name, ' ID:', num2str(probeIDs(probe_i)), display_append];
            plot(time, field(probe_i,:), 'LineStyle', style, 'DisplayName', display_name);
        end
    end
    xlim([min(time),max(time)]);
    
    %%%
    % Add annotations: legend and labels.
    %%%
    
    legend();
    
    y_labels = {'x-Velocity (m/s)', ...
                'y-Velocity (m/s)', ...
                'z-Velocity (m/s)', ...
                'Temperature (K)', ...
                'Speed (m/s)', ...
                'Mach', ...
                'Density (kg/m3)', ...
                'Pressure (pa)', ...
                'Total Pressure (Pa)', ...
                'Dynamic Pressure (Pa)', ...
                'Kinematic Viscosity'};
	ylabel(y_labels{field_index});
    xlabel(x_label);
    
    hleg = legend('show');
    set(hleg, 'Location', 'southwest');
    
    % Reset color plotting order for the next time the plot is drawn to.
    set(hax1,'ColorOrderIndex',1);
    
end

