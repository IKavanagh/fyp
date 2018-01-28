% Create VEFIE elements (V = ZE, Z = I + GD) for defined area and materials
tic;

% Incident field
Einc = zeros(problem_size, 3);
E0inc = 1; % Amplitude of current/incident wave

l = lambda_d / 4;
Ie = E0inc*0.01;

% Contrast function, contains material information
X = zeros(problem_size, 1);

% Greens vector
G = zeros(problem_size, 1);
G(1) = (exp(-1i*k0*equiv_a)*(1+1i*k0*equiv_a) - 1) / (delta_v*k0*k0);

for counter = 1:problem_size
    k = wave_number(counter);
    
    X(counter) = (k*k) / (k0*k0) - 1;
    
    % Plane wave
    Einc(counter, 1:3) = [E0inc*exp(-1i*k0*position(counter, 3)), 0 + 1i*0, 0 + 1i*0];
    
    % Electric Dipole
%     r = sqrt(sum(position(counter, 1:3).^2));
%     
%     Einc(counter, 1) = eta_d * ((Ie*l*cos(theta(counter)) / (2*pi*(r^2)))) * exp(-1i*k*r); % rho
%     Einc(counter, 2) = 1i*eta_d * ((Ie*l*cos(theta(counter)) / (2*pi*(r^2)))) * exp(-1i*k*r); % theta
%     Einc(counter, 3) = 0; % phi
    
    if counter ~= 1
        Rtemp = sqrt(sum((position(1, 1:3) - position(counter, 1:3)).^2));
        G(counter) = (exp(-1i*k0*Rtemp)) / (4.0*pi*Rtemp);
    end
end

% Convert from spherical to cartesian co-ordinates
% Einc_x = Einc(:,1) .* sin(theta) .* cos(phi) + ...
%          Einc(:,2) .* cos(theta) .* cos(phi) - ...
%          Einc(:,3) .* sin(phi);
%          
% Einc_y = Einc(:,1) .* sin(theta) .* sin(phi) + ...
%          Einc(:,2) .* cos(theta) .* sin(phi) + ...
%          Einc(:,3) .* cos(phi);
%          
% Einc_z = Einc(:,1) .* cos(theta) - ...
%          Einc(:,2) .* sin(theta);
%          
% Einc = [Einc_x, Einc_y, Einc_z];

t = toc;
fprintf(1, 'Created Einc, X and G for a problem size of %.0f in %.4fs\n\n', problem_size, t);