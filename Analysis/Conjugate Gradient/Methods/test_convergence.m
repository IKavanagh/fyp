% Test the convergence of CG, BICG and CGNE
clear all
close all

N = 1e2; % Size of matrices
tol = 1e-6; % Tolerance to solve for
points = N / 10; % Distance between consecutive points when plotting

complex = 1; % Whether matrix should be complex or real
posdef = 0; % Whether matrix should be positive definite or not

[A, b] = create_matrices(N, complex, posdef);

% Run CG
[x1, r1] = cg(A, b, tol, points);

% Run BiCG
[x2, r2] = bicg(A, b, tol, points);

% Run CGNE
[x3, r3] = cgne(A, b, tol, points);

it1 = 0:points:(length(r1)-1)*points;
it2 = 0:points:(length(r2)-1)*points;
it3 = 0:points:(length(r3)-1)*points;

h = figure;
hold on
plot(it1, log10(r1), 'r-');
plot(it2, log10(r2), 'b-.');
plot(it3, log10(r3), 'g--');
title('Plot of Convergence Rate of CG Methods');
xlabel('Number of Iterations');
ylabel('Log10 Norm of Residual');
axis([0 2*N log10(tol*0.1) log10(1e5)]);
legend('CG', 'BiCG', 'CGNE');
hold off