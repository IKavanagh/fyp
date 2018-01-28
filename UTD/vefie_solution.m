% Solve V = ZE, where Z = I + GD, using CGNE-FFT
tic;


do_comparison = 0 ; 
tol = 1e-5;
max_it_counter = 1000; 

n = 0;

rfo = logical(real(D)); % Reduced forward operator
Vred = rfo.*V;

% Create an initial guess
E = zeros(problem_size, 1);
the_residual  = zeros(problem_size,1) ;

E1 = zeros(problem_size,1) ; 
E2 = zeros(problem_size,1) ; 
E3 = zeros(problem_size,1) ; 
E4 = zeros(problem_size,1) ; 


Dconj = conj(D);
Gconj = conj(G.');

rn_1 = (rfo.*E + rfo.*conv_fft(G, D.*E, N)) - Vred; % r0 = Z*E - V

p = -(rfo.*rn_1 + Dconj.*conv_fft(Gconj, rfo.*rn_1, N)); % p0 = -Z'*r0


while (norm(rn_1) > tol && n <= max_it_counter)
    n = n + 1;

     the_residual(n,1) = log10(norm(rn_1)) ; 
    
    if (mod(n, 10) == 0)
        fprintf(1, 'Using reduced operator: Number of iterations %d \n', n);
        fprintf(1, 'Error is %2.5f - Aiming for %2.5f \n',norm(rn_1),tol );
    end

    
    
    alpha = (norm(rfo.*rn_1 + Dconj.*conv_fft(Gconj, rfo.*rn_1, N)) / norm(rfo.*p + rfo.*conv_fft(G, D.*p, N)))^2; % norm(Z'*rn_1)^2 / norm(Z*p)^2;
    
    E = E + alpha * p;
    
    rn = rn_1 + alpha*(rfo.*p + rfo.*conv_fft(G, D.*p, N)); % rn = rn_1 + alpha*Z*p;
    
    beta = (norm(rfo.*rn + Dconj.*conv_fft(Gconj, rfo.*rn, N)) / norm(rfo.*rn_1 + Dconj.*conv_fft(Gconj, rfo.*rn_1, N)))^2; % beta = norm(Z'*rn)^2 / norm(Z'*rn_1)^2;
    
    p = -(rfo.*rn + Dconj.*conv_fft(Gconj, rfo.*rn, N)) + beta * p; % p = -Z'*rn + beta*p;

    % Move rn back to rn_1
    rn_1 = rn;
    
    % Make an intermediate result - to check out convergence 
    if( n == 0.2*max_it_counter ) 
        E1 = E + ~rfo.*(V - conv_fft(G, D.*E, N));
    end
   if( n == 0.4*max_it_counter ) 
        E2 = E + ~rfo.*(V - conv_fft(G, D.*E, N));
    end
    if( n == 0.6*max_it_counter ) 
        E3 = E + ~rfo.*(V - conv_fft(G, D.*E, N));
    end
   if( n == 0.8*max_it_counter ) 
        E4 = E + ~rfo.*(V - conv_fft(G, D.*E, N));
    end
    
        
end

% Solve V = Z*E for points ignored by use of reduced forward operator
E = E + ~rfo.*(V - conv_fft(G, D.*E, N));
n_reduced = n; 

t = toc;
fprintf(1, 'Solved VEFIE for a problem size of %.0f in %.4f seconds and %.0f iterations\n', problem_size, t, n);

E_reduced = E ; 




if do_comparison == 1 

    max_it_counter = n_reduced; 
n = 0;

rfo = ones(problem_size,1) ; % Reduced forward operator = identity, so no reduction. 

Vred = rfo.*V;

% Create an initial guess
E = zeros(problem_size, 1);
 
non_reduced_residual = zeros(problem_size,1) ;

Dconj = conj(D);
Gconj = conj(G.');

rn_1 = (rfo.*E + rfo.*conv_fft(G, D.*E, N)) - Vred; % r0 = Z*E - V

p = -(rfo.*rn_1 + Dconj.*conv_fft(Gconj, rfo.*rn_1, N)); % p0 = -Z'*r0


while (norm(rn_1) > tol && n <= max_it_counter)
    n = n + 1;

    non_reduced_residual(n,1) = log10(norm(rn_1)) ; 
    
    if (mod(n, 10) == 0)
        fprintf(1, 'Using non-reduced operator: Number of iterations %d \n', n);
        fprintf(1, 'Error is %2.5f - Aiming for %2.5f \n',norm(rn_1),tol );
    end
    
    alpha = (norm(rfo.*rn_1 + Dconj.*conv_fft(Gconj, rfo.*rn_1, N)) / norm(rfo.*p + rfo.*conv_fft(G, D.*p, N)))^2; % norm(Z'*rn_1)^2 / norm(Z*p)^2;
    
    E = E + alpha * p;
    
    rn = rn_1 + alpha*(rfo.*p + rfo.*conv_fft(G, D.*p, N)); % rn = rn_1 + alpha*Z*p;
    
    beta = (norm(rfo.*rn + Dconj.*conv_fft(Gconj, rfo.*rn, N)) / norm(rfo.*rn_1 + Dconj.*conv_fft(Gconj, rfo.*rn_1, N)))^2; % beta = norm(Z'*rn)^2 / norm(Z'*rn_1)^2;
    
    p = -(rfo.*rn + Dconj.*conv_fft(Gconj, rfo.*rn, N)) + beta * p; % p = -Z'*rn + beta*p;

    % Move rn back to rn_1
    rn_1 = rn;
    
    
        
end
end



if do_comparison == 1
figure
plot(the_residual(1:n_reduced,1)) ;
hold on
plot(non_reduced_residual(1:n_reduced,1),'r') ;
xlabel('Iterations') 
ylabel(' Error ') 
legend('Reduced','Non-reduced') 
end


% Flatten out solution for plotting
E_total = zeros(M, N);
E1_total = zeros(M, N);
E2_total = zeros(M, N);
E3_total = zeros(M, N);
E4_total = zeros(M, N);

E_inc = zeros(M, N);

E_total_vefie_real = zeros(M, N);
E_total_vefie_imag = zeros(M, N);

E_inc_vefie_real = zeros(M, N);
E_inc_vefie_imag = zeros(M, N);

for y = 1:M
    for x = 1:N
        position_counter = (y - 1)*N + x;
        
        E_total(y, x) = E_reduced(position_counter);
        
        E_inc(y, x) = V(position_counter);
   end
end