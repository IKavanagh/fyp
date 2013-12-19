function [x, r] = cg(A, b, tol, points)
    % Solves Ax = b using the preconditioned Conjugate Gradient method
    tic;
    
    N = length(b);
    
    % Calculates preconditioner as I
    M = eye(N, N);

    n = 0;

    % Create initial guess x0
    x0 = zeros(N, 1);
    xn = x0;

    % Compute r0 = b - Ax
    rn = b - A*x0;

    r = zeros(N, 1);
    r(1) = norm(rn);
    while (norm(rn) > tol && n < 2*N)
        n = n + 1;
        
        % Preconditioner calculations %
        % Solve Mz^(i-1) = r^(i-1)
        % z = inv(M) * r;
        zn = M \ rn;

        rhon_1 = rn'*zn;

        if (n == 1)
            p = zn;
        else
            beta = rhon_1 / rhon_2;
            p = zn + beta * p;
        end

        q = A*p;
        alpha = rhon_1 / ((p') * q);
        xn = xn + alpha * p;
        rn = rn - alpha * q;

        % Move rho_n_1 to rho_n_2
        rhon_2 = rhon_1;

        if (mod(n, points) == 0)
            if (mod(n, 100) == 0)
                fprintf(1, 'Number of iterations %.f \n', n);
            end
            r((n / points) + 1) = norm(rn);
        end
    end
    x = xn;

    t = toc;
    fprintf(1, 'Solved Ax = b using CG in %.4f seconds and %d iterations\n', t, n);
end