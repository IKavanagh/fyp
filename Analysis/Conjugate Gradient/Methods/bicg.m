function [x, r] = bicg(A, b, tol, points)
    % Solves Ax = b using the preconditioned BiConjugate Gradient method
    tic;

    N = length(b);
    
    % Calculates preconditioner as I
    M = eye(N, N);

    n = 0;

    % Create initial guess x0
    x0 = zeros(N, 1);
    xn = x0;

    % Compute r0 = b - Ax0
    rn = b - A*x0;
    rtn = rn;

    At = A';

    r = zeros(N, 1);
    r(1) = norm(rn);
    while (norm(rn) > tol && n < 2*N)
        n = n + 1;
        
        % Preconditioner calculations %
        % Solve Mz^(i-1) = r^(i-1)
        % z = inv(M) * r;
        zn = M \ rn;
        
        ztn = (M') \ rtn;

        rho_i_1 = zn'*rtn;
        if (rho_i_1 == 0)
            % KABOOM
            fprintf(1, 'BiCG Exploded');
            break
        end

        if (n == 1)
            p = zn;
            pt = ztn;
        else
            beta = rho_i_1 / rho_i_2;
            p = zn + beta * p;
            pt = ztn + beta * pt;
        end

        q = A*p;
        qt = At*pt;
        alpha = rho_i_1 / ((pt') * q);
        xn = xn + alpha * p;
        rn = rn - alpha * q;
        rtn = rtn - alpha * qt;

        % Move rho_i_1 to rho_i_2
        rho_i_2 = rho_i_1;

        if (mod(n, points) == 0)
            if (mod(n, 100) == 0)
                fprintf(1, 'Number of iterations %.f \n', n);
            end
            r((n / points) + 1) = norm(rn);
        end
    end
    x = xn;

    t = toc;
    fprintf(1, 'Solved Ax = b using BiCG in %.4f seconds and %d iterations\n', t, n);
end