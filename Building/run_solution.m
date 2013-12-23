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