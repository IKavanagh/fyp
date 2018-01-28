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
% vefie_solution
bicgstab_solution

% Run Mie Series solution
mie_series_solution

num_figs = 10; % Plot around 10 figures

step = round((M - 1) / num_figs);

for counter = 1:step:M % Loop through y axis
    h = figure;
    hold on
    title('Comparison of electric fields produced by Mie series and VEFIE');
    plot(real(position(1:N)), E_total_vefie_real(counter, 1:end), 'b-');
    plot(real(position(1:N)), E_total_mie_real(counter, 1:end), 'r-.');

    legend('VEFIE', 'Mie Series');

    xlabel('x-axis');
    ylabel('Real Component');
    hold off
    
    str = sprintf('p_mie_real_%.3fr_%dx_%dy_%ddisc_%.2fer_%dyaxis.fig', radius, x_side, y_side, disc_per_lambda, epsilonr, counter);
    
    hgsave(h, str);
    
    h = figure;
    hold on
    title('Comparison of electric fields produced by Mie series and VEFIE');
    plot(real(position(1:N)), E_total_vefie_imag(counter, 1:end), 'b-');
    plot(real(position(1:N)), E_total_mie_imag(counter, 1:end), 'r-.');

    legend('VEFIE', 'Mie Series');

    xlabel('x-axis');
    ylabel('Imaginary Component');
    hold off
    
    str = sprintf('p_mie_imag_%.3fr_%dx_%dy_%ddisc_%.2fer_%dyaxis.fig', radius, x_side, y_side, disc_per_lambda, epsilonr, counter);
    
    hgsave(h, str);
end