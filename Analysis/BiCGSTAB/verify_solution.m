% Run VEFIE solution for scatterer
clear all
close all

% Load variables
variables

% Create shape for problem to be solved over
create_shape

% Create VEFIE elements
create_mie_elements

% Run VEFIE solution
bicgstab_solution

% Run Mie Series solution
mie_series_solution

num_figs = 5; % Plot around 5 figures

step = round((M - 1) / num_figs);

for counter = 1:step:M % Loop through y axis
    h = figure;
    hold on
    title('Comparison of electric fields produced by Mie series and VEFIE');
    plot(E_total_vefie_real(counter, 1:end), 'b');
    plot(E_total_mie_real(counter, 1:end), 'g');

    legend('VEFIE', 'Mie Series');

    xlabel('x-axis');
    ylabel('Re(Ez)');
    hold off
end