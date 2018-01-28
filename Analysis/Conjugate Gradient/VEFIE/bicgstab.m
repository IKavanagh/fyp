function [E, SE] = bicgstab(Z, V, tol, points)
    % Solves V = ZE using the BiConjugate Gradient Stabilised
    tic;

    N = length(V);

    n = 0;
    
    % Create an initial guess
    E = zeros(N, 1);
    
    r = V - Z*E;
    
    rho_n = 1;
    alpha = 1;
    w = 1;

    v = zeros(N, 1);
    p = zeros(N, 1);

    r0 = r / sqrt(r'*r); % Arbitrary choice such that <r0, r> ~= 0

    error = sqrt(r'*r)/sqrt(V'*V); % e = ||r||2 / ||V||2
    
    SE = zeros(N, 1);
    SE(1) = error;
    while (error > tol && n < N)
        n = n + 1;
        
        rho_n_1 = rho_n;
    
        rho_n = r0'*r;

        beta = (rho_n/rho_n_1) * (alpha/w);

        p = r + beta*(p-w*v);

        v = Z*p; % v = Z*p

        alpha = rho_n / (r0'*v);

        s = r - alpha*v;

        t = Z*s; % t = Z*s

        w = (t'*s) / (t'*t); % w = <t, s> / <t, t>

        E = E + alpha*p + w*s;

        r = s - w*t;

        error = sqrt(r'*r)/sqrt(V'*V); % e = ||r||2 / ||V||2
        
        if (mod(n, points) == 0)
            if (mod(n, 100) == 0)
                fprintf(1, 'Number of iterations %.f \n', n);
            end
            SE((n / points) + 1) = norm(rn);
        end
        
    end
        
    t = toc;
    fprintf(1, 'Solved Ax = b using BiCGSTAB in %.4f seconds and %d iterations\n', t, n);
end