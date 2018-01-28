% Solve V = ZE, where Z = I + GD, using CGNE-FFT
n = 0;

% Create an initial guess
E = zeros(problem_size, 1);

rn_1 = (E + conv_fft(G, D.*E, N)) - V; % r0 = Z*E - V

p = -(rn_1 + conj(D).*conv_fft(conj(G.'), rn_1, N)); % p0 = -Z'*r0

while (norm(rn_1) > tol && n < problem_size)
    n = n + 1;

    if (mod(n, 10) == 0)
        fprintf(1, 'Number of iterations %.f \n', n);
    end
    
    alpha = (norm(rn_1 + conj(D).*conv_fft(conj(G.'), rn_1, N)) / norm(p + conv_fft(G, D.*p, N)))^2; % norm(Z'*rn_1)^2 / norm(Z*p)^2;
    
    E = E + alpha * p;
    
    rn = rn_1 + alpha*(p + conv_fft(G, D.*p, N)); % rn = rn_1 + alpha*Z*p;
    
    beta = (norm(rn + conj(D).*conv_fft(conj(G.'), rn, N)) / norm(rn_1 + conj(D).*conv_fft(conj(G.'), rn_1, N)))^2; % beta = norm(Z'*rn)^2 / norm(Z'*rn_1)^2;
    
    p = -(rn + conj(D).*conv_fft(conj(G.'), rn, N)) + beta * p; % p = -Z'*rn + beta*p;

    % Move rn back to rn_1
    rn_1 = rn;
end