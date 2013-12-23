% Create VEFIE elements (V = ZE, Z = I + GD) for defined area and materials
tic;

% Incident field
V = zeros(problem_size, 1);

% Diagonal contrast matrix, contains material information
D = zeros(problem_size, 1);

% Dense matrix, contains scatterer interaction information
G = zeros(1, problem_size);

leading_G_factor = ((2*pi*equiv_a) / k0);
for counter = 1:problem_size
    
    %if(mod(counter,100) == 0) 
    %     fprintf(1,'Counter equals %d of %d \n', counter, problem_size); 
    %end
    
    % besselh(0, 2, 0) = NaN -> Need to put antennae in a position other
    % than the centre of a cell
    R = abs(position(counter, 1) - antennae);
    V(counter, 1) = besselh(0, 2, k0 * R); %H0(2)(k0*R)
    
    k = wave_number(counter, 1);
    D(counter, 1) = (k*k - k0*k0);
    
    if (counter == 1) % Self term
        G(1, counter) = (1i / 4) * (leading_G_factor * besselh(1, 2, k0*equiv_a) - 4.0*1i / (k0*k0));
    else
        Rmn = abs(position(1) - position(counter));
        
        G(1, counter) = (1i / 4) * (leading_G_factor * besselj(1, k0 * equiv_a) * besselh(0, 2, k0*Rmn));
    end
end

t = toc;
fprintf(1, 'Created V, G and D for a problem size of %.0f in %.4f seconds\n', problem_size, t);