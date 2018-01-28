function [E, SE] = cgne(V, Z, tol)
    % Solves V = ZE using the Conjugate Gradient applied to the Normal
    % Equations method
    tic;

    N = length(V);

    n = 0;

    Zt = Z';

    % Create an initial guess
    E = zeros(N, 1);
    
    rn_1 = Z*E - V; % r0
    p = -Zt * rn_1; % p0

    SE = zeros(N, 1);
    SE(1) = norm(rn_1);
    while (norm(rn_1) > tol && n < N)
        n = n + 1;

        alpha = norm(Zt * rn_1)^2 / norm(Z * p)^2;
        E = E + alpha * p;

        rn = rn_1 + alpha * Z * p;

        beta = norm(Zt * rn)^2 / norm(Zt * rn_1)^2;

        p = -Zt * rn + beta * p;

        % Move rn back to rn_1
        rn_1 = rn;
        
        if (mod(n, 10) == 0)
            fprintf(1, 'Number of iterations %.f \n', n);
        end
        
        SE(n+1) = norm(rn_1);
    end

    t = toc;
    fprintf(1, 'Solved Ax = b using CGNE in %.4f seconds and %d iterations\n', t, n);
end