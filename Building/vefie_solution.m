% Solve V = ZE, where Z = I + GD, using CGNE-FFT
tic;

n = 0;
tol = 1e-3;

rfo = logical(real(D)); % Reduced forward operator
Vred = rfo.*V;

% Create an initial guess
E = zeros(problem_size, 1);

Dconj = conj(D);
Gconj = conj(G.');

rn_1 = (rfo.*E + rfo.*conv_fft(G, D.*E, N)) - Vred; % r0 = Z*E - V

p = -(rfo.*rn_1 + Dconj.*conv_fft(Gconj, rfo.*rn_1, N)); % p0 = -Z'*r0

while (norm(rn_1) > tol && n < problem_size)
    n = n + 1;

    if (mod(n, 100) == 0)
        fprintf(1, 'Number of iterations %.f \n', n);
    end
    
    alpha = (norm(rfo.*rn_1 + Dconj.*conv_fft(Gconj, rfo.*rn_1, N)) / norm(rfo.*p + rfo.*conv_fft(G, D.*p, N)))^2; % norm(Z'*rn_1)^2 / norm(Z*p)^2;
    
    E = E + alpha * p;
    
    rn = rn_1 + alpha*(rfo.*p + rfo.*conv_fft(G, D.*p, N)); % rn = rn_1 + alpha*Z*p;
    
    beta = (norm(rfo.*rn + Dconj.*conv_fft(Gconj, rfo.*rn, N)) / norm(rfo.*rn_1 + Dconj.*conv_fft(Gconj, rfo.*rn_1, N)))^2; % beta = norm(Z'*rn)^2 / norm(Z'*rn_1)^2;
    
    p = -(rfo.*rn + Dconj.*conv_fft(Gconj, rfo.*rn, N)) + beta * p; % p = -Z'*rn + beta*p;

    % Move rn back to rn_1
    rn_1 = rn;
end

% Solve V = Z*E for points ignored by use of reduced forward operator
E = E + ~rfo.*(V - conv_fft(G, D.*E, N));

t = toc;
fprintf(1, 'Solved VEFIE for a problem size of %.0f in %.4f seconds and %.0f iterations\n', problem_size, t, n);

% Flatten out solution for plotting
E_total = zeros(M, N);
E_inc = zeros(M, N);

E_total_vefie_real = zeros(M, N);
E_total_vefie_imag = zeros(M, N);

E_inc_vefie_real = zeros(M, N);
E_inc_vefie_imag = zeros(M, N);

for y = 1:M
    for x = 1:N
        position_counter = (y - 1)*N + x;
        
        E_total(y, x) = E(position_counter);
        E_inc(y, x) = V(position_counter);
   
        E_total_vefie_real(y, x) = real(E(position_counter));
        E_total_vefie_imag(y, x) = imag(E(position_counter));

        E_inc_vefie_real(y, x) = real(V(position_counter));
        E_inc_vefie_imag(y, x) = imag(V(position_counter));
    end
end