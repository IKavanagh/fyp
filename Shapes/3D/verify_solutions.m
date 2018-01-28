% Run VEFIE solution for 3D sphere embedded in a cube and compare with Mie
% Series
clear all
close all

% Load variables
variables

% Create shape
create_shape

% Create VEFIE elements
create_vefie_elements

E = zeros(problem_size, 3);

% Run VEFIE solution
vefie_solution

% Verify results with Mie Series

% Set up variables for Mie Series and plotting results
% Mie Series only valid outside of scatterer

% Keep x-axis constant
num_xaxis = floor((abs(start_point(1)) - radius) / delta_x); % Number of x-axis lines outside of sphere on each side

xaxis = randi([(-num_xaxis+1) num_xaxis]); % Randomise selection for fun

% Ensure valid value either side of sphere
if (xaxis < 1)
    xaxis = xaxis + N;
end

num_figs = 10; % Plot ~5 graphs

step = round((M - 1) / num_figs);

fprintf(1, '\n');

for yaxis = 1:step:M
    selection_start = (yaxis - 1)*N + xaxis;
    selection = selection_start:N*M:problem_size; % Indices of positions to plot
    selection_size = length(selection);

    selection_positions(1:O, 1:3) = position(selection, :);

    if (sum(X(selection)) ~= 0)
        fprintf(1, 'Selection intersects sphere! Mie Series is invalid');
    end

    % Spherical co-ordinates for selection
    rho_selection(1:selection_size, 1) = rho(selection);
    phi_selection(1:selection_size, 1) = phi(selection); 
    theta_selection(1:selection_size, 1) = theta(selection); 

    % Create VEFIE fields for selection
    Einc_vefie = Einc(selection, :);
    Et_vefie = E(selection, :);
    Es_vefie = Et_vefie - Einc_vefie;

    % Run Mie Series solution
    mie_series_solution

    % Plot results
    h = figure;

    fprintf(1, 'Plot %d along x-axis line %.0f, y-axis line %.0f\n', h, xaxis, yaxis);

    hold on
    title('Comparison of electric fields produced by Mie series and VEFIE');

    plot(selection_positions(:,3), real(Et_vefie(:,1)), 'b-', selection_positions(:,3), real(Et_mie(:,1)), 'r-.');
    
    legend('VEFIE', 'Mie Series');

    xlabel('z-axis');
    ylabel('Real Component');
    hold off
    
    str = sprintf('p_3d_mie_real_verif_%.3fr_%dx_%dy_%dz_%ddisc_%.2fer_%dxaxis_%dyaxis.fig', radius, x_side, y_side, z_side, disc_per_lambda, epsilonr, xaxis, yaxis);
    
    hgsave(h, str);

    % Plot results
    h = figure;

    fprintf(1, 'Plot %d along x-axis line %.0f, y-axis line %.0f\n', h, xaxis, yaxis);

    hold on
    title('Comparison of electric fields produced by Mie series and VEFIE');

    plot(selection_positions(:,3), imag(Et_vefie(:,1)), 'b-', selection_positions(:,3), imag(Et_mie(:,1)), 'r-.');
    
    legend('VEFIE', 'Mie Series');

    xlabel('z-axis');
    ylabel('Imaginary Component');
    hold off
    
    str = sprintf('p_3d_mie_imag_verif_%.3fr_%dx_%dy_%dz_%ddisc_%.2fer_%dxaxis_%dyaxis.fig', radius, x_side, y_side, z_side, disc_per_lambda, epsilonr, xaxis, yaxis);
    
    hgsave(h, str);
end