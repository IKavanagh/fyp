% General Variables %
MHz = 1000000.0; 
f = 700.0*MHz; 
omega = 2.0*pi*f;

tol = 1e-3;

% Define Mie Series limits
% Depend heavily on frequency in question
high_limit = 100;
low_limit = -high_limit;

% Materials %
% Free Space (Vacuum) %
epsilon0 = 8.854e-12; % Permittivity
mu0 = 4.0*pi*1.0e-7; % Permeability
c0 = 1.0/sqrt(epsilon0*mu0); % Speed of EM wave (light)
k0 =  omega*sqrt(mu0*epsilon0); % Wave number
eta0 =  sqrt(mu0/epsilon0);
lambda0 = c0/f; % Wavelength

% Homogenous Cylnder's Material %
epsilonrd = 10 + (0.0*1i);
epsilond = epsilonrd*epsilon0; 
mud = mu0; % No magnetic permeability in concrete
c_d = 1.0/sqrt(epsilond*mud); 
kd = omega*sqrt(mud*epsilond);
etad =  sqrt(mud/epsilond);
lambdad = c_d/f;

% Shape %
% Containing Rectangle
centre = 0.0 + 0.0 * 1i;
% x_side = 1; % Defined outside this script
% y_side = 1; % Defined outside this script

% Cylinder %
radius = min([x_side y_side]) * 0.5 * 0.5;

start_point = centre - (x_side * 0.5 + y_side * 0.5 * 1i);

disc_per_lambda = 10; % Recommended peterson (Pg. 62)

% Determine N, multiple of 2
N = floor(x_side / (abs(lambdad) / disc_per_lambda));
while (mod(N, 2) ~= 0)
    N = N + 1;
end
delta_x = x_side / N;

% Determine M, multiple of 4
M = floor(y_side / (abs(lambdad) / disc_per_lambda));
while (mod(M, 2) ~= 0)
    M = M + 1;
end
delta_y = y_side / M;

problem_size = N*M;

% Determine radius (a) of equivalent circle for discretised area
equiv_a = sqrt(delta_x * delta_y / pi);

% Determine area of discretised area
area = delta_x * delta_y;