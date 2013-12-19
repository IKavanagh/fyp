% Investigate the use of the Mie Series for an incident wave by comparing
% the standard exp(-jk0x) representation to the Mie Series representation
% of an infinite sum of bessel functions
clear all
close all

% Load variables
variables

% Incident wave travelling from (-10, 0) to (10, 0)
E_inc = zeros(1, 1001);
E_inc_mie = zeros(1, 1001);

for p = 1:1001
   x(p) = -5 + 0.01*(p - 1);
   
   E_inc(p) = exp(-1i*k0*x(p));
   
   E_inc_mie(p) = 0;
   for n = low_limit:high_limit
      [phi, rho] = cart2pol(x(p), 0); 
      
      E_inc_mie(p) = E_inc_mie(p) + (power(1i, n)*besselj(n, k0*rho)*exp(1i*n*phi)); 
   end
end

figure
plot(x, real(E_inc), x, real(E_inc_mie));
xlabel('x-axis');
ylabel('Re(Incident Field)');
title('Comparison of Incident Fields');
legend('Exact', 'Mie Series');