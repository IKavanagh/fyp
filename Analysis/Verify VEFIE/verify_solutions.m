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

x_axis(1, 1:N) = real(position(1:N)); % Needed here because Peterson overwrites this
x_limit = roundto(max(abs(x_axis)), 2);

% Run Peterson's VEFIE Solution
petersons_solution

% Plot results
for counter = 1:M % Loop through y axis
    h = figure;
    hold on
    title('Comparison of Real Part of Total Electric Field');
    
    plot(x_axis, E_total_vefie_real(counter, :), 'b');
    plot(x_axis, E_total_mie_real(counter, :), 'g');
    plot(x_axis, E_total_vefie_peter_real(counter, :), 'r');

    xlim([-x_limit +x_limit]);
    legend('VEFIE', 'Mie Series', 'Petersons VEFIE');
    xlabel('x-axis');
    ylabel('Re(Ez)');
    hold off
end