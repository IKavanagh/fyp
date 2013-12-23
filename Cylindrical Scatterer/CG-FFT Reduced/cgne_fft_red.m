% Solve V = ZE, where Z = I + GD, using CGNE-FFT
n = 0;
tol = 1e-6;

rfo = logical(D); % Reduced forward operator
Vred = rfo.*V;

% Create an initial guess
E = zeros(problem_size, 1);

rn_1 = (rfo.*E + rfo.*conv_fft(G, D.*E, N)) - Vred; % r0 = Z*E - V

p = -(rfo.*rn_1 + conj(D).*conv_fft(conj(G.'), rfo.*rn_1, N)); % p0 = -Z'*r0

while (norm(rn_1) > tol && n < problem_size)
    n = n + 1;

    if (mod(n, 100) == 0)
        fprintf(1, 'Number of iterations %.f \n', n);
    end
    
    alpha = (norm(rfo.*rn_1 + conj(D).*conv_fft(conj(G.'), rfo.*rn_1, N)) / norm(rfo.*p + rfo.*conv_fft(G, D.*p, N)))^2; % norm(Z'*rn_1)^2 / norm(Z*p)^2;
    
    E = E + alpha * p;
    
    rn = rn_1 + alpha*(rfo.*p + rfo.*conv_fft(G, D.*p, N)); % rn = rn_1 + alpha*Z*p;
    
    beta = (norm(rfo.*rn + conj(D).*conv_fft(conj(G.'), rfo.*rn, N)) / norm(rfo.*rn_1 + conj(D).*conv_fft(conj(G.'), rfo.*rn_1, N)))^2; % beta = norm(Z'*rn)^2 / norm(Z'*rn_1)^2;
    
    p = -(rfo.*rn + conj(D).*conv_fft(conj(G.'), rfo.*rn, N)) + beta * p; % p = -Z'*rn + beta*p;

    % Move rn back to rn_1
    rn_1 = rn;
end

% Solve V = Z*E for points ignored by use of reduced forward operator
E = E + ~rfo.*(V - conv_fft(G, D.*E, N));