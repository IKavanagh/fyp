% General Variables %
MHz = 1000000.0; 
f = 1500.0*MHz; 
omega = 2.0*pi*f;

% Define Mie Series limit
% Depends on frequency in question
mie_limit = 100;

% Materials %
% Free Space (Vacuum) %
epsilon0 = 8.854e-12; % Permittivity
mu0 = 4.0*pi*1.0e-7; % Permeability
c0 = 1.0/sqrt(epsilon0*mu0) ; % Speed of EM wave (light)
k0 = omega*sqrt(mu0*epsilon0) ; % Wave number
eta0 = sqrt(mu0/epsilon0) ;
lambda0 = c0 / f ; % Wavelength

% Dielectric Scatterer %
epsilonr = 2.62 - 1i;
epsilon_d = epsilonr*epsilon0;
mu_d = mu0;
k_d = omega*sqrt(epsilon_d*mu_d);
eta_d = sqrt(mu_d/epsilon_d);
c_d = 1.0/sqrt(epsilon_d*mu_d); 
lambda_d = c_d / f;

% Shape %
disc_per_lambda = 10; % Recommended peterson (Pg. 62)

% Containing Square/Rectangle
x_side = 1;
y_side = 1;

centre = 0.0 + 0.0 * 1i;

start_point = centre - (x_side * 0.5 + y_side * 0.5 * 1i);

% Determine N, multiple of 2
N = floor(x_side / (abs(lambda_d) / disc_per_lambda));
while (mod(N, 2) ~= 0)
    N = N + 1;
end
delta_x = x_side / N;

% Determine M, multiple of 2
M = floor(y_side / (abs(lambda_d) / disc_per_lambda));
while (mod(M, 2) ~= 0)
    M = M + 1;
end
delta_y = y_side / M;

problem_size = N*M;

% Determine radius (a) of equivalent circle for discretised area
equiv_a = sqrt(delta_x * delta_y / pi);