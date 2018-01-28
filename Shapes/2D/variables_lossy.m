% Load variables
variables

% Dielectric Scatterer %
sigma_d = 0;

propagation_constant_d = sqrt(1i*omega*mu_d*(sigma_d + 1i*omega*epsilon_d));

k_d = -1i*propagation_constant_d;
eta_d = 1i*omega*mu_d/propagation_constant_d;
lambda_d = 2*pi/(real(k_d));