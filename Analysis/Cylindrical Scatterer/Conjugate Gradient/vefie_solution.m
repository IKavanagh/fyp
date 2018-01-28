% Solves V = ZE, where Z = I + GD, using the CGNE
tic;

E = cgne(Z, V, 1e-6);

t = toc;
fprintf(1, 'Solved VEFIE for a problem size of %.0f in %.4f seconds\n', problem_size, t);

E_total_vefie_real = zeros(M, N);
E_total_vefie_imag = zeros(M, N);

E_inc_vefie_real = zeros(M, N);
E_inc_vefie_imag = zeros(M, N);

for y = 1:M
    for x = 1:N
        position_counter = (y - 1)*N + x;
   
        E_total_vefie_real(y, x) = real(E(position_counter));
        E_total_vefie_imag(y, x) = imag(E(position_counter));

        E_inc_vefie_real(y, x) = real(V(position_counter));
        E_inc_vefie_imag(y, x) = imag(V(position_counter));
    end
end