function [x, r] = cgne(A, b, tol, points)
    % Solves Ax = b using the Conjugate Gradient applied to the Normal
    % Equations method
    tic;

    N = length(b);

    n = 0;

    At = A';

    % Create an initial guess
    x0 = zeros(N, 1);

    xn = x0;
    rn_1 = A * x0 - b; % r0
    p = -At * rn_1; % p0

    r = zeros(N, 1);
    r(1) = norm(rn_1);
    while (norm(rn_1) > tol && n < 2*N)
        n = n + 1;

        alpha = norm(At * rn_1)^2 / norm(A * p)^2;
        xn = xn + alpha * p;

        rn = rn_1 + alpha * A * p;

        beta = norm(At * rn)^2 / norm(At * rn_1)^2;

        p = -At * rn + beta * p;

        % Move rn back to rn_1
        rn_1 = rn;
        
        if (mod(n, points) == 0)
            if (mod(n, 100) == 0)
                fprintf(1, 'Number of iterations %.f \n', n);
            end
            r((n / points) + 1) = norm(rn);
        end
    end
    x = xn;

    t = toc;
    fprintf(1, 'Solved Ax = b using CGNE in %.4f seconds and %d iterations\n', t, n);
end