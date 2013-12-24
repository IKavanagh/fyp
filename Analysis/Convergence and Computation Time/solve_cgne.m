% Solve V = ZE, where Z = I + GD, using CGNE
n = 0;

% Create an initial guess
E = zeros(problem_size, 1);

Zt = Z';

rn_1 = Z*E - V; % r0 = Z*E - V

p = -Zt*rn_1; % p0 = -Z'*r0

while (norm(rn_1) > tol && n < problem_size)
    n = n + 1;

    if (mod(n, 10) == 0)
        fprintf(1, 'Number of iterations %.f \n', n);
    end
    
    alpha = (norm(Zt*rn_1) / norm(Z*p))^2; % norm(Z'*rn_1)^2 / norm(Z*p)^2;
    
    E = E + alpha * p;
    
    rn = rn_1 + alpha*Z*p; % rn = rn_1 + alpha*Z*p;
    
    beta = (norm(Zt*rn) / norm(Zt*rn_1))^2; % beta = norm(Z'*rn)^2 / norm(Z'*rn_1)^2;
    
    p = -Zt*rn + beta * p; % p = -Z'*rn + beta*p;

    % Move rn back to rn_1
    rn_1 = rn;
end