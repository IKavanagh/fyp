% Create VEFIE elements (V = ZE, Z = I + GD) for defined area and materials
tic;

% Incident field
V = zeros(problem_size, 1);

I = eye(problem_size, problem_size);

% Diagonal contrast matrix, contains material information
D = zeros(problem_size, problem_size);

% Dense matrix, contains scatterer interaction information
G = zeros(1, problem_size);

leading_G_factor = ((2*pi*equiv_a) / k0);
for counter = 1:problem_size
    % besselh(0, 2, 0) = NaN -> Need to put antennae in a position other
    % than the centre of a cell
    R = abs(position(counter, 1) - antenna_location);
    V(counter, 1) = besselh(0, 2, k0 * R); %H0(2)(k0*R)
    
    D(counter, counter) = (wave_number(counter, 1)*wave_number(counter, 1) - k0*k0);
    
    % problem_size refers to the cell that is interacting on other cells
    for row = 1:problem_size
        % G has different values for M = N and M ~= N
        if (row == counter) % M = N => Self term
            G(row, counter) = (1i / 4) * (leading_G_factor * besselh(1, 2, k0*equiv_a) - 4.0*1i / (k0*k0));
        else % M ~= N => Interaction of cell (row) on cell (counter)
            Rmn = abs(position(row) - position(counter));
            
            G(row, counter) = (1i / 4) * (leading_G_factor * besselj(1, k0 * equiv_a) * besselh(0, 2, k0*Rmn));
        end
        
    end
end

Z = I + G*D;

t = toc;
fprintf(1, 'Created V, G and D for a problem size of %.0f in %.4f seconds\n', problem_size, t);