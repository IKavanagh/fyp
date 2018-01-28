% Run VEFIE solution for scatterer
clear all
close all
clc

% Load variables
variables
% variables_lossy

% Create shape for problem to be solved over
create_shape

% Create VEFIE elements
create_vefie_elements

% Run VEFIE solution
vefie_solution
% bicgstab_solution

analytic_solution

%% Plot results
max_val = max(max(log10(abs(E_total))));
min_val = min(min(log10(abs(E_total))));

for y = 1:N
    for x = 1:N

        if (scatterer(y,x) == 1) 
            scatterer(y,x) = max_val; 
        end

        if (scatterer(y,x) == 0)
            scatterer(y,x) = min_val; 
        end

    end
end

if 1 == 2 % Total electric field with Scatterer shown
    h = figure;
    
    surf(real(pos), imag(pos), log10(abs(E_total)));
    hold on
    surf(real(pos), imag(pos), scatterer);
    
    shading interp;
    view(2);
    xlabel('Distance (m)');
    ylabel('Distance (m)');
    title('Electric Field Throughout Scatterer');
    hold off
end

if 1 == 2 % Total electric field without Scatterer shown
    h = figure;
    
    surf(real(pos), imag(pos), log10(abs(E_total)));
    
    shading interp;
    view(2);
    xlabel('Distance (m)');
    ylabel('Distance (m)');
    title('Electric Field Throughout Scatterer');
end

if 1 == 2 % Incident field with Scatterer shown
    h = figure;
    
    surf(real(pos), imag(pos), log10(abs(E_inc)));
    hold on
    surf(real(pos), imag(pos), scatterer);
    
    shading interp;
    view(2);
    xlabel('Distance (m)');
    ylabel('Distance (m)');
    title('Plot of Incident Electric Field');
    hold off
end

if 1 == 2 % Incident field without Scatterer shown
    h = figure;
    
    surf(real(pos), imag(pos), log10(abs(E_inc)));
    
    shading interp;
    view(2);
    xlabel('Distance (m)');
    ylabel('Distance (m)');
    title('Plot of Incident Electric Field');
end

if 1 == 2 % Movie
    h = figure;
    
    start_time = 0;
    end_time = 10;
    step_time = (end_time - start_time) / 200;
    
    t = start_time:step_time:end_time;
    
    X = real(pos);
    Y = imag(pos);
    
    for i = 1:length(t)
        wave = exp(1i*t(i));
        
        Z = E_total*wave;

        surf(X, Y, real(Z)); % Doesn't work with scatterer shown
        
        shading interp;
        view(2);
        xlabel('Distance (m)');
        ylabel('Distance (m)');
        title('Plot of Total Electric Field Interacting on Scatterer');
        
        movie_frames(i) = getframe(gcf);
    end
    
    %movie2avi(movie_frames, 'scatterer.avi');
end