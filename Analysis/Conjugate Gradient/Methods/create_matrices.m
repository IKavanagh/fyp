function [A, b] = create_matrices(N, complex, posdef)
    % Creates random matrix A and b for use in solving Ax = b using the
    % Conjugate Gradient method

    % Make A and b matrices
    % A must be symmetric positve definite
    rng('shuffle');

    if (complex == 0)
        A = randn(N, N);
        b = randn(N, 1);
    else
        A = randn(N, N) + randn(N, N)*1i;
        b = randn(N, 1) + randn(N, 1)*1i;
    end

    if (posdef == 1)
        % Symmetric positive definite
        A = A+A'; 
    end

end