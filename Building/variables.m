% Define general, material and shape variables

% General Variables %
MHz = 1000000.0; 
f = 700.0*MHz; 
omega = 2.0*pi*f;

% Materials %
% Free Space (Vacuum) %
epsilon0 = 8.854e-12; % Permittivity
mu0 = 4.0*pi*1.0e-7; % Permeability
c0 = 1.0/sqrt(epsilon0*mu0); % Speed of EM wave (light)
k0 =  omega*sqrt(mu0*epsilon0); % Wave number
eta0 =  sqrt(mu0/epsilon0);
lambda0 = c0/f; % Wavelength

% Concrete %
epsilonrc = 6.0 - 1i;
epsilonc = epsilonrc*epsilon0;
muc = mu0; % No magnetic permeability in concrete
cc = 1.0/sqrt(epsilonc*muc);
kc = omega*sqrt(epsilonc*muc);
etac = sqrt(muc/epsilonc);
lambdac = cc/f;

% sigmac = 100000.001 ;
% muc = mu0;
% epsilonc = 3.5*epsilon0;
% 
% propagation_constant_c = sqrt(1i*omega*muc*(sigmac +1i*omega*epsilonc));
% 
% kc = -1i*propagation_constant_c; 
% etac = 1i*omega*muc/propagation_constant_c; 
% lambdac = 2*pi/(real(kc));

% Glass %
epsilonrg = 4.0;
epsilong = epsilonrg*epsilon0;
mug = mu0; % No magnetic permeability in glass
cg = 1.0/sqrt(epsilong*mug);
kg = omega*sqrt(epsilong*mug);
etag = sqrt(mug/epsilong);
lambdag = cg/f;

% Wood %
epsilonrw = 4.0;
epsilonw = epsilonrw*epsilon0;
muw = mu0; % No magnetic permeability in wood
cw = 1.0/sqrt(epsilonw*muw);
kw = omega*sqrt(epsilonw*muw);
etaw = sqrt(muw/epsilonw);
lambdaw = cw/f;

% Define a general lambda for dielectrics
lambdad = lambdac; % Use free space because it has a smaller lambda

disc_per_lambda = 10; % Recommended peterson (Pg. 62)

% Shape %
% Outline %
x_side = 2; % metres
y_side = 2; % metres

wall = min([x_side y_side]) / 40; % Roughly 0.25 metres
door = min([x_side y_side]) / 20; % Roughly 0.5 metres

% Define antennae location
antennae = -0.8 -0.8*1i;