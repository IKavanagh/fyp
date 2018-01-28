% Solve the Mie Series for a defined area and materials
tic;

E_inc_mie = zeros(problem_size, 1);
E_total_mie = zeros(problem_size, 1);

for counter = 1:problem_size
    if(mod(counter,100) == 0) 
        fprintf(1,'Counter equals %d of %d \n', counter, problem_size); 
    end
    
   k = wave_number(counter, 1);
   
   if (k == k0)
      % We are outside of scatterer and total E field is made up of
      % incident field and Mie Series scattering field
      E_scat_mie_sum = 0;
      E_inc_mie_sum = 0;
      
      % NB: k below should be k_d... 
      for n = -100:100
         An_top = (eta0 / eta_d) * besselj(n, k0*radius) * 0.5 * (besselj(n-1, k_d*radius) - besselj(n+1, k_d*radius)) - 0.5 * (besselj(n-1, k0*radius) - besselj(n+1, k0*radius)) * besselj(n, k_d*radius);
         
         An_bottom = besselj(n, k_d*radius) * 0.5 * (besselh(n-1, 2, k0*radius) - besselh(n+1, 2, k0*radius)) - (eta0 / eta_d) * 0.5 * (besselj(n-1, k_d*radius) - besselj(n+1, k_d*radius)) * besselh(n, 2, k0*radius);
         
         An = An_top / An_bottom;
          
         E_scat_mie_sum = E_scat_mie_sum + power(1i,-n) * An * besselh(n, 2, k0*rho(counter)) * exp(1i*n*phi(counter));
         
         E_inc_mie_sum = E_inc_mie_sum + power(1i,-n)*besselj(n, k0*rho(counter))*exp(1i*n*phi(counter));
      end
      
      E_inc_mie(counter, 1) = E_inc_mie_sum;
      E_total_mie(counter, 1) = E_inc_mie_sum + E_scat_mie_sum;
   else
      % We are inside scatterer and total E field is given by Mie Series
      E_total_mie_sum = 0;
      E_inc_mie_sum = 0;
      
      for n = -100:100
          Bn_top = -2*1i;
          
          Bn_bottom = (pi*k0*radius) * (besselj(n, k_d*radius) * 0.5 * (besselh(n-1, 2, k0*radius) - besselh(n+1, 2, k0*radius)) - (eta0 / eta_d) * 0.5 * (besselj(n-1, k_d*radius) - besselj(n+1, k_d*radius)) * besselh(n, 2, k0*radius));
          
          Bn = Bn_top / Bn_bottom;
          
          E_total_mie_sum = E_total_mie_sum + power(1i,-n) * Bn * besselj(n, k_d*rho(counter)) * exp(1i*n*phi(counter));
         
          E_inc_mie_sum = E_inc_mie_sum + (power(1i, n)*besselj(n, k0*rho(counter))*exp(1i*n*phi(counter)));
      end
      
      E_inc_mie(counter, 1) = E_inc_mie_sum;
      E_total_mie(counter, 1) = E_total_mie_sum;
   end
end

t = toc;
fprintf(1, 'Solved Mie Series for a problem size of %.0f in %.4f seconds\n', problem_size, t);

E_total_mie_real = zeros(M, N);
E_total_mie_imag = zeros(M, N);

E_inc_mie_real = zeros(M, N);
E_inc_mie_imag = zeros(M, N);

for y = 1:M
    for x = 1:N
        position_counter = (y - 1)*N + x;
        
        E_total_mie_real(y, x) = real(E_total_mie(position_counter));
        E_total_mie_imag(y, x) = imag(E_total_mie(position_counter));
        
        E_inc_mie_real(y, x) = real(E_inc_mie(position_counter));
        E_inc_mie_imag(y, x) = imag(E_inc_mie(position_counter));
    end
end