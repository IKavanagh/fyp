function [E, SE] = cgne_fft(V, D, G, N, tol)
    % Solve V = ZE, where Z = I + GD, using CGNE-FFT
    tic;
    
    n = 0;
    
    problem_size = length(V);
    
    SE = zeros(problem_size, 1);

    % Create an initial guess
    E = zeros(problem_size, 1);

    Dconj = conj(D);
    Gconj = conj(G.');

    % Begin CG
    rn_1 = (E + conv_fft(G, D.*E, N)) - V; % r0 = Z*E - V

    p = -(rn_1 + Dconj.*conv_fft(Gconj, rn_1, N)); % p0 = -Z'*r0
    
    error = sqrt(rn_1'*rn_1)/sqrt(V'*V); % e = ||rn_1||2 / ||V||2
    
    SE(1) = error;

    while (error > tol && n < problem_size)
        n = n + 1;

        if (mod(n, 10) == 0)
            fprintf(1, 'Number of iterations %.f \n', n);
            fprintf(1, 'Error is %2.8f - Aiming for %2.8f \n', error, tol);
        end

        alpha = (norm(rn_1 + Dconj.*conv_fft(Gconj, rn_1, N)) / norm(p + conv_fft(G, D.*p, N)))^2; % norm(Z'*rn_1)^2 / norm(Z*p)^2;

        E = E + alpha * p;

        rn = rn_1 + alpha*(p + conv_fft(G, D.*p, N)); % rn = rn_1 + alpha*Z*p;

        beta = (norm(rn + Dconj.*conv_fft(Gconj, rn, N)) / norm(rn_1 + Dconj.*conv_fft(Gconj, rn_1, N)))^2; % beta = norm(Z'*rn)^2 / norm(Z'*rn_1)^2;

        p = -(rn + Dconj.*conv_fft(Gconj, rn, N)) + beta * p; % p = -Z'*rn + beta*p;

        % Move rn back to rn_1
        rn_1 = rn;
        
        error = sqrt(rn_1'*rn_1)/sqrt(V'*V); % e = ||rn_1||2 / ||V||2
        
        SE(n+1) = error;
        
        if (log10(error) <= -6)
            fprintf(1, 'Number of iterations %.0f\n', n');
        end
    end
    
    SE = SE(SE~=0); % Remove additional zeros

    t = toc;
    fprintf(1, 'Average time per iteration is %.4fs\n', (t / n));
    fprintf(1, 'Solved VEFIE with CGNE-FFT for a problem size of %.0f in %.4f seconds and %.0f iterations\n', problem_size, t, n);
end