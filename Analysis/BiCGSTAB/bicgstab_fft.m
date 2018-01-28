function [E, SE] = bicgstab_fft(V, D, G, N, tol)
    % Solve V = ZE, where Z = I + GD, using CGNE-FFT
    tic;
    
    n = 0;
    
    problem_size = length(V);
    
    SE = zeros(problem_size, 1);

    % Create an initial guess
    E = zeros(problem_size, 1);

    r = V - E + conv_fft(G, D.*E, N); % r = V - Z*E

    rho_n = 1;
    alpha = 1;
    w = 1;

    v = zeros(problem_size, 1);
    p = zeros(problem_size, 1);

    r0 = r / sqrt(r'*r); % Arbitrary choice such that <r0, r> ~= 0

    error = sqrt(r'*r)/sqrt(V'*V); % e = ||r||2 / ||V||2
    
    SE(1) = error;

    while (error > tol && n < problem_size)
        n = n + 1;

%         if (mod(n, 10) == 0)
%             fprintf(1, 'Number of iterations %.f \n', n);
%             fprintf(1, 'Error is %2.3f - Aiming for %2.3f \n', norm(rn_1), tol);
%         end

        rho_n_1 = rho_n;
    
        rho_n = r0'*r;

        beta = (rho_n/rho_n_1) * (alpha/w);

        p = r + beta*(p-w*v);

        v = p + conv_fft(G, D.*p, N); % v = Z*p

        alpha = rho_n / (r0'*v);

        s = r - alpha*v;

        t = s + conv_fft(G, D.*s, N); % t = Z*s

        w = (t'*s) / (t'*t); % w = <t, s> / <t, t>

        E = E + alpha*p + w*s;

        r = s - w*t;

        error = sqrt(r'*r)/sqrt(V'*V); % e = ||r||2 / ||V||2
        
        SE(n+1) = error;
        
        if (log10(error) <= -6)
            SE(n+1) = 0;
            fprintf(1, 'Number of iterations %.0f\n', n');
        end
    end
    
    SE = SE(SE~=0); % Remove additional zeros

    t = toc;
    fprintf(1, 'Solved VEFIE with BiCGSTAB-FFT for a problem size of %.0f in %.4f seconds and %.0f iterations\n', problem_size, t, n);
end