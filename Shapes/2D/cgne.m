function [E, SE] = cgne(V, Z, tol)
    % Solve V = ZE, where Z = I + GD, using CGNE-FFT
    tic;
    
    n = 0;
    
    problem_size = length(V);
    
    SE = zeros(problem_size, 1);

    % Create an initial guess
    E = zeros(problem_size, 1);

    Zt = Z';

    rn_1 = Z*E - V; % r0 = Z*E - V

    p = -Zt*rn_1; % p0 = -Z'*r0
    
    error = sqrt(rn_1'*rn_1)/sqrt(V'*V); % e = ||rn_1||2 / ||V||2
    SE(1) = error;

    while (norm(rn_1) > tol && n < problem_size)
        n = n + 1;

        alpha = (norm(Zt*rn_1) / norm(Z*p))^2; % norm(Z'*rn_1)^2 / norm(Z*p)^2;

        E = E + alpha * p;

        rn = rn_1 + alpha*Z*p; % rn = rn_1 + alpha*Z*p;

        beta = (norm(Zt*rn) / norm(Zt*rn_1))^2; % beta = norm(Z'*rn)^2 / norm(Z'*rn_1)^2;

        p = -Zt*rn + beta * p; % p = -Z'*rn + beta*p;

        % Move rn back to rn_1
        rn_1 = rn;
        
        error = sqrt(rn_1'*rn_1)/sqrt(V'*V); % e = ||rn_1||2 / ||V||2
        SE(n+1) = error;
    end
    
    SE = SE(SE~=0); % Remove additional zeros

    t = toc;
    fprintf(1, 'Solved VEFIE with CGNE for a problem size of %.0f in %.4f seconds and %.0f iterations\n', problem_size, t, n);
end