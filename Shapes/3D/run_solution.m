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

% Plot results

% Keep x-axis constant
num_xaxis = floor((abs(start_point(1)) - radius) / delta_x); % Number of x-axis lines outside of sphere on each side

xaxis = randi([(-num_xaxis+1) num_xaxis]); % Randomise selection for fun

% Ensure valid value either side of sphere
if (xaxis < 1)
    xaxis = xaxis + N;
end

num_figs = 5; % Plot ~5 graphs

step = round((M - 1) / num_figs);

fprintf(1, '\n');

for yaxis = 1:step:M
    selection_start = (yaxis - 1)*N + xaxis;
    selection = selection_start:N*M:problem_size; % Indices of positions to plot

    selection_positions(1:O, 1:3) = position(selection, :);

    % Create VEFIE fields for selection
    Einc_vefie = Einc(selection, :);
    Et_vefie = E(selection, :);

    % Plot results
    h = figure;

    fprintf(1, 'Plot %d along x-axis line %.0f, y-axis line %.0f\n', h, xaxis, yaxis);

    plot(selection_positions(:,3), Einc_vefie(:,1), selection_positions(:,3), Et_vefie(:,1));
    
    title('Plot of Fields');
    xlabel('z-axis');
    ylabel('Re(E)');
end