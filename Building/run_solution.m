% Run VEFIE solution for building
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
if 1 == 1 % Total electric field
    h = figure;
    
    surf(real(pos), imag(pos), abs(E_total));
    shading interp;
    view(2);
    xlabel('Distance (m)');
    ylabel('Distance (m)');
    title('Plot of Total Electric Field Interacting in Building');
end

if 1 == 1 % Incident field
    h = figure;
    
    surf(real(pos), imag(pos), abs(E_inc));
    shading interp;
    view(2);
    xlabel('Distance (m)');
    ylabel('Distance (m)');
    title('Plot of Incident Electric Field in Building');
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
        Et = E_total + 1i*E_total;
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
    
    movie2avi(movie_frames, 'building.avi');
end