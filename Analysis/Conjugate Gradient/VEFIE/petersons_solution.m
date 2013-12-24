% Solve Peterson's VEFIE (V = ZE) for a defined area and materials
tic;

% Only use discretised areas fully within the scatterer
new_problem_size = 0;

for counter = 1:problem_size
   k = wave_number(counter, 1);
   
   % Only want cells fully within the scatterer
   if (k == kd)
       new_problem_size = new_problem_size + 1;
       new_position(new_problem_size, 1) = position(counter);
   end
end

position = new_position;
problem_size = new_problem_size;

% Create matrices for VEFIE (V = ZE, Z = (I + GD))
% Incident field
V = zeros(problem_size, 1);

Z = zeros(problem_size, problem_size);

leading_Z_factor = (eta0 * pi * equiv_a) / 2;

for counter = 1:problem_size
    if( mod(counter, 100) == 0 ) 
         fprintf(1,'Counter equals %d of %d \n', counter, problem_size); 
    end
    
    V(counter, 1) = exp(-1i*k0*real(position(counter)));
    
    % problem_size refers to the cell that is interacting on other cells
    for row = 1:problem_size
        % Z has different values for M = N and M ~= N
        if (row == counter) % M = N => Self term
            
            Z(row, counter) = leading_Z_factor * (besselj(1, k0 * equiv_a) - 1i * bessely(1, k0 * equiv_a)) - (1i*eta0*epsilonrd) / (k0*(epsilonrd - 1));
        else % M ~= N => Interaction of cell (row) on cell (counter)
            Rmn = abs(position(row) - position(counter));
            
            Z(row, counter) = leading_Z_factor * besselj(1, k0 * equiv_a) * (besselj(0, k0*Rmn) - 1i * bessely(0, k0*Rmn));
        end
        
    end
end

t = toc;
fprintf(1, 'Created V and Z = (I + GD) for a problem size of %.0f in %.4f seconds\n', problem_size, t);

% Solves V = ZJ, using J = Z^-1*V
tic;

J = inv(Z)*V;

E = J ./ (1i*omega*epsilon0*(epsilonrd - 1));

t = toc;
fprintf(1, 'Solved Petersons VEFIE for a problem size of %.0f in %.4f seconds\n', problem_size, t);

E_total_vefie_peter_real = zeros(M, N);
E_total_vefie_peter_imag = zeros(M, N);

E_inc_vefie_peter_real = zeros(M, N);
E_inc_vefie_peter_imag = zeros(M, N);

% Put E back into original problem and convert to x and y
new_E = zeros(M, N);

counter = 0;

for y = 1:M
    for x = 1:N
        % Determine if we are inside or outside of scatterer
        if (scatterer(y, x) == 1)
            % Inside scatterer we have a result
            counter = counter + 1;
            
            E_total_vefie_peter_real(y, x) = real(E(counter));
            E_total_vefie_peter_imag(y, x) = real(E(counter));
   
            E_inc_vefie_peter_real(y, x) = real(V(counter));
            E_inc_vefie_peter_imag(y, x) = imag(V(counter));
        end
    end
end