function [E, SE] = bicgstab_fft(Einc, X, G, position, N, M, O, delta_x, delta_y, delta_z, k0, tol)
    % Sinc = L[E] where Einc = E - iwA - Eirr
    tic;
    
    n = 0;
    
    problem_size = length(Einc);
    
    SE = zeros(problem_size, 1);

    % Required constants
    delta_v = delta_x * delta_y * delta_z;
    
    Co = delta_v;
    kb_2 = k0*k0;

    % Create an initial guess
    E = zeros(problem_size, 3);

    A = Co*conv_fft_3D(G, [X, X, X].*E, N, M, O);
    Eirr = central_finite_differencing(A, position, N, M, O, delta_x, delta_y, delta_z);

    L = E - kb_2*A - Eirr; % Linear operator, L = L[E]

    % BiCGSTAB Initiation
    r = Einc - L;

    rho_n = 1;
    alpha = 1;
    w = 1;

    v = zeros(problem_size, 3);
    p = zeros(problem_size, 3);
    s = zeros(problem_size, 3);

    r0 = r / sqrt(r(:)'*r(:)); % Arbitrary choice such that <r0, r> ~= 0

    error = sqrt(r(:)'*r(:))/sqrt(Einc(:)'*Einc(:)); % e = ||r||2 / ||Einc||2
    
    SE(1) = error;

    while (error > tol && n < problem_size)
        n = n + 1;

%         if (mod(n, 10) == 0)
%             fprintf(1, 'Number of iterations %.f \n', n);
%             fprintf(1, 'Error is %2.3f - Aiming for %2.3f \n', norm(rn_1), tol);
%         end

        rho_n_1 = rho_n;

        rho_n = r0(:)'*r(:); % rho = <r0, r>
        beta = (rho_n/rho_n_1) * (alpha/w);

        p = r + beta*(p-w*v);

            A = Co*conv_fft_3D(G, [X, X, X].*p, N, M, O);
            Eirr = central_finite_differencing(A, position, N, M, O, delta_x, delta_y, delta_z);
            L = p - kb_2*A - Eirr;

        v = L; % v = L[p]

        alpha = rho_n / (r0(:)'*v(:)); % alpha = rho / <r0, v>

        s = r - alpha*v;

            A = Co*conv_fft_3D(G, [X, X, X].*s, N, M, O);
            Eirr = central_finite_differencing(A, position, N, M, O, delta_x, delta_y, delta_z);
            L = s - kb_2*A - Eirr;

        t = L; % t = L[s]

        w = (t(:)'*s(:)) / (t(:)'*t(:)); % w = <t, s> / <t, t>

        E = E + alpha*p + w*s;

        r = s - w*t;

        error = sqrt(r(:)'*r(:))/sqrt(Einc(:)'*Einc(:)); % e = ||r||2 / ||Einc||2
        
        SE(n+1) = error;
    end
    
    SE = SE(SE~=0); % Remove additional zeros

    t = toc;
    fprintf(1, 'Solved VEFIE with BiCGSTAB-FFT for a problem size of %.0f in %.4f seconds and %.0f iterations\n', problem_size, t, n);
end