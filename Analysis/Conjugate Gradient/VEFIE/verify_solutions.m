% Run all solutions and plot outputs from all on the same graphs
clear all
close all

% Load variables
variables

% Create shape for problem to be solved over
create_shape

% Run Mie Series solution
mie_series_solution

% Create VEFIE elements
create_vefie_elements

% Run VEFIE solution
vefie_solution

% Run Peterson's VEFIE Solution
petersons_solution

% Plot results
for counter = 1:M % Loop through y axis
    h = figure;
    hold on
    title('Comparison of Real Part of Total Electric Field');
    plot(E_total_vefie_real(counter, 1:end), 'b');
    plot(E_total_mie_real(counter, 1:end), 'g');
    plot(E_total_vefie_peter_real(counter, 1:end), 'r');
    
    xlabel('x-axis');
    ylabel('Re(Ez)');
    hold off
end