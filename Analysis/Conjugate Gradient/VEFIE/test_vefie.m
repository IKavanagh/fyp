% Run all solutions and plot outputs from all on the same graphs
clear all
close all

% Load variables
variables

tol = 1e-6; % Tolerance to solve for
points = 5; % Distance between consecutive points when plotting

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
plot(it1, log10(r1), 'r--', it2, log10(r2), 'b:', it3, log10(r3), 'g-');

xlim([0 2*N]);
ylim([log10(1e-3) log10(1e4)]);

legend('CG', 'BiCG', 'CGNE');
titleStr = sprintf('Plot of Convergence Rate of CG Methods\nVEFIE Problem');
title(titleStr);
xlabel('Number of Iterations');
ylabel('Order of Error');