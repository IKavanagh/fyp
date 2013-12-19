% Run all solutions and plot outputs from all on the same graphs
clear all
close all

% Load variables
variables

tol = 1e-6; % Tolerance to solve for
points = 1; % Distance between consecutive points when plotting

% Create shape for problem to be solved over
create_shape

% Run VEFIE solution
create_vefie_elements

% Run CG
[e1, r1] = cg(Z, V, tol, points);

% Run BiCG
[e2, r2] = bicg(Z, V, tol, points);

% Run CGNE
[e3, r3] = cgne(Z, V, tol, points);

it1 = 0:points:(length(r1)-1)*points;
it2 = 0:points:(length(r2)-1)*points;
it3 = 0:points:(length(r3)-1)*points;

h = figure;
hold on
plot(it1, log10(r1), 'r-');
plot(it2, log10(r2), 'b-.');
plot(it3, log10(r3), 'g--');
title('Plot of Convergence Rate of CG Methods for VEFIE Problem');
xlabel('Number of Iterations');
ylabel('Log10 Norm of Residual');
legend('CG', 'BiCG', 'CGNE');
hold off