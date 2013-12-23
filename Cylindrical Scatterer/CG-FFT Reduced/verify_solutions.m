% Run all solutions and plot outputs from all on the same graphs
clear all
close all

% Load variables
variables

% Create shape for problem to be solved over
create_shape

% Create VEFIE elements
create_vefie_elements

% Run VEFIE solution
vefie_solution

% Plot results
if 1 == 2
    
    % Run Mie Series solution
    mie_series_solution
    
    for counter = 1:M % Loop through y axis
        h = figure;
        hold on
        title('Comparison of Real Part of Total Electric Field');
        plot(E_total_vefie_real(counter, 1:end), 'b');
        plot(E_total_mie_real(counter, 1:end), 'g');

        xlabel('x-axis');
        ylabel('Re(Ez)');
        hold off
    end
end

% Surface plot of total electric field
if 1 == 1
    h = figure;
    
    temp = max([x_side y_side]);
    axis_lim = [-0.5*temp +0.5*temp];
    
    surf(real(pos), imag(pos), scatterer);
    shading interp;
    view(2);
    xlabel('Distance (m)');
    ylabel('Distance (m)');
    title('Plot of Cylindrical Scatterer');
    xlim(axis_lim);
    ylim(axis_lim);
    
    %filename = sprintf('p_scatterer_r%0.f', radius);
    %hgsave(h, filename);
    
    h = figure;
    
    surf(real(pos), imag(pos), abs(E_total_vefie_real + 1i*E_total_vefie_imag));
    shading interp;
    view(2);
    xlabel('Distance (m)');
    ylabel('Distance (m)');
    title('Plot of Total Electric Field Interacting on Cylindrical Scatterer');
    xlim(axis_lim);
    ylim(axis_lim);
    
    %filename = sprintf('p_electric_field_r%0.f', radius);
    %hgsave(h, filename);
end

% Movie
if 1 == 1
    h = figure;
    
    start_time = 0;
    end_time = 10;
    step_time = (end_time - start_time) / 200;
    
    t = start_time:step_time:end_time;
    
    X = real(pos);
    Y = imag(pos);
    
    for i = 1:length(t)
        Et = E_total_vefie_real + 1i*E_total_vefie_imag;
        wave = exp(1i*t(i));
        
        Z = Et*wave;

        surf(X, Y, real(Z));
        shading interp;
        view(2);
        xlabel('Distance (m)');
        ylabel('Distance (m)');
        title('Plot of Total Electric Field Interacting on Cylindrical Scatterer');
        
        movie_frames(i) = getframe(gcf);
    end
    
    movie2avi(movie_frames, 'scatterer.avi');
end