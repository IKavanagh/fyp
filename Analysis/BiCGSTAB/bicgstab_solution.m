% Solve V = ZE, where Z = I + GD, using BICGSTAB-FFT
tic;

n = 0;
tol = 1e-3;

% Create an initial guess
E = zeros(problem_size, 1);

r = V - E + conv_fft(G, D.*E, N); % r = V - Z*E

rho_n = 1;
alpha = 1;
w = 1;

v = zeros(problem_size, 1);
p = zeros(problem_size, 1);
s = zeros(problem_size, 1);

r0 = r / sqrt(r'*r); % Arbitrary choice such that <r0, r> ~= 0

error = sqrt(r'*r)/sqrt(V'*V); % e = ||r||2 / ||V||2

while (error > tol && n < problem_size)
    n = n + 1;
    
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

    fprintf(1, 'Iteration %.f \n', n);
    fprintf(1, 'Error is %2.3f - Aiming for %2.3f \n', error, tol);
end

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