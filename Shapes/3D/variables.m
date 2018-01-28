% General Variables %
MHz = 1000000.0; 
f = 1500.0*MHz; 
omega = 2.0*pi*f;

% Define Mie Series limit
% Depend heavily on frequency in question
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
epsilonr = 2.62;
epsilon_d = epsilonr*epsilon0;
mu_d = mu0;
k_d = omega*sqrt(epsilon_d*mu_d);
eta_d = sqrt(mu_d/epsilon_d);
c_d = 1.0/sqrt(epsilon_d*mu_d); 
lambda_d = c_d / f;

% Shape %
disc_per_lambda = 10; % Recommended peterson (Pg. 62)

% Containing Cube
x_side = 1;
y_side = 1;
z_side = 1;

centre = [0.0, 0.0, 0.0]; % x, y, z

% Sphere %
radius_wrt_lambda = 2;
radius = radius_wrt_lambda * lambda_d;

if (x_side < 2*radius)
    x_side = 3*radius;
end

if (y_side < 2*radius)
    y_side = 3*radius;
end

if (z_side < 2*radius)
    z_side = 3*radius;
end

start_point = centre - [x_side * 0.5, y_side * 0.5, z_side * 0.5];

% N (M or O) is number of discretisations along x (y or z) side
N = floor((x_side / lambda_d) * disc_per_lambda);
while (mod(N, 2) ~= 0)
    N = N + 1;
end
delta_x = x_side / N;

M = floor((y_side / lambda_d) * disc_per_lambda);
while (mod(M, 2) ~= 0)
    M = M + 1;
end
delta_y = y_side / M;

O = floor((z_side / lambda_d) * disc_per_lambda);
while (mod(O, 2) ~= 0)
    O = O + 1;
end
delta_z = z_side / O;

problem_size = N*M*O;

% Determine radius (a) of equivalent circle for discretised area
equiv_a = nthroot(0.75*(delta_x*delta_y*delta_z) / pi, 3);

% Determine volume of discretised cube
delta_v = delta_x * delta_y * delta_z;