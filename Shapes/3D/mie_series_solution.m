% Solve the Mie Series for a defined area and materials
a = radius;

Einc_sum = zeros(selection_size, 3);
Es_sum = zeros(selection_size, 3);

i_term = 1.0;
minus_term = -1.0 ;

for n = 1:mie_limit
    
    i_term = i_term*1j;
    n_term = (2.0 * n + 1.0) / (n*(n+1));
    minus_term = minus_term * (-1) ; 
    
    % Legendre values
    legendre_values = legendre(n,cos(theta_selection)).'; % 2nd column => Order 1
    n_minus_1_legendre_values = legendre(n-1,cos(theta_selection)).';
    
    if (n == 1) 
        legendre_deriv  = n*cos(theta_selection) .* legendre_values(:,2);
    else
        legendre_deriv = n*cos(theta_selection) .* legendre_values(:,2) - (n+1)*n_minus_1_legendre_values(:,2);
    end
    
    legendre_deriv  = legendre_deriv ./ sqrt(1.0 - cos(theta_selection).*cos(theta_selection)) ;          
    
    % Bessel values (spherical)
    bessel_k0R = sqrt(pi ./ (2*k0*rho_selection)) .* besselj(n+0.5, k0*rho_selection); % Vector
    bessel_k0R_deriv = 0.5*bessel_k0R + 0.5*sqrt(pi*k0*rho_selection / 2.0) .* ...
                       (besselj(n-0.5, k0*rho_selection) - besselj(n+1.5, k0*rho_selection));
    
    bessel_kda = sqrt(pi / (2*k_d*a)) * besselj(n+0.5,k_d*a); % Scalar 
    bessel_kda_deriv = 0.5*bessel_kda + 0.5*sqrt(pi*k_d*a / 2.0) * ...
                       (besselj(n-0.5, k_d*a) - besselj(n+1.5, k_d*a) );

    bessel_k0a = sqrt(pi/(2*k0*a)) * besselj(n+0.5, k0*a); % Scalar 
    bessel_k0a_deriv = 0.5*bessel_k0a + 0.5*sqrt(pi*k0*a / 2.0) * ...
                       (besselj(n-0.5, k0*a) - besselj(n+1.5, k0*a));

    % Hankel values (spherical)
    hankel_k0R = sqrt(pi ./ (2*k0*rho_selection)) .* besselh(n+0.5, 2, k0*rho_selection); % Vector
    hankel_k0R_deriv = 0.5 * hankel_k0R + 0.5*sqrt(pi*k0*rho_selection/2.0) .* ...
                       (besselh(n-0.5, 2, k0*rho_selection) - besselh(n+1.5, 2, k0*rho_selection));
    
    hankel_k0a = sqrt(pi / (2*k0*a)) * besselh(n+0.5, 2, k0*a); % Scalar
    hankel_k0a_deriv = 0.5*hankel_k0a + 0.5*sqrt(pi*k0*a / 2.0) * ...
                       (besselh(n-0.5, 2, k0*a) - besselh(n+1.5, 2, k0*a));
    
    % m odd 1 n vector (Einc)
    m_odd_rho_1 = zeros(selection_size, 1);
    m_odd_theta_1 = (1.0 ./ sin(theta_selection)) .* bessel_k0R .* legendre_values(:,2) .* cos(phi_selection);
    m_odd_phi_1 = -1*bessel_k0R .* legendre_deriv .* sin(phi_selection);
    m_odd_1 = [m_odd_rho_1, m_odd_theta_1, m_odd_phi_1];
    
    % n even 1 n vector (Einc)
    n_even_rho_1 = (n*(n+1) ./ (k0*rho_selection)) .* bessel_k0R .* legendre_values(:,2) .* cos(phi_selection);
    n_even_theta_1 = (1.0 ./ (k0*rho_selection)) .* bessel_k0R_deriv .* legendre_deriv .* sin(phi_selection);
    n_even_phi_1 = (-1.0 ./ (k0*rho_selection .* sin(theta_selection))) .* bessel_k0R_deriv .* legendre_values(:,2) .* sin(phi_selection);
    n_even_1 = [n_even_rho_1,n_even_theta_1,n_even_phi_1];
    
    % m odd 3 n vector (Es)
    m_odd_rho_3 = zeros(selection_size, 1);
    m_odd_theta_3 = (1.0 ./ sin(theta_selection)) .* hankel_k0R .* legendre_values(:,2) .* cos(phi_selection);
    m_odd_phi_3 = -1*hankel_k0R .* legendre_deriv .* sin(phi_selection);
    m_odd_3 = [m_odd_rho_3, m_odd_theta_3, m_odd_phi_3];
    
    % n even 3 n vector (Es)
    n_even_rho_3 = (n*(n+1) ./ (k0*rho_selection)) .* hankel_k0R .* legendre_values(:,2) .* cos(phi_selection);
    n_even_theta_3 = (1.0 ./ (k0*rho_selection)) .* hankel_k0R_deriv .* legendre_deriv .* sin(phi_selection);
    n_even_phi_3 = (-1.0./(k0*rho_selection .* sin(theta_selection))) .* hankel_k0R_deriv .* legendre_values(:,2) .* sin(phi_selection);
    n_even_3 = [n_even_rho_3, n_even_theta_3, n_even_phi_3];
    
    % a_r_n and b_r_n
    a_r_n = -(mu_d * bessel_kda * bessel_k0a_deriv - mu0 * bessel_k0a * bessel_kda_deriv) / ...
            (mu_d * bessel_kda * hankel_k0a_deriv - mu0 * hankel_k0a * bessel_kda_deriv);
        
    b_r_n = -(mu_d * bessel_k0a * bessel_kda_deriv - mu0 * ((k_d^2) / (k0^2)) * bessel_kda * bessel_k0a_deriv) / ...
            (mu_d * hankel_k0a * bessel_kda_deriv - mu0 * ((k_d^2) / (k0^2)) * bessel_kda * hankel_k0a_deriv);
    
    % Create incident and reflected the fields (needed for Er)
    Einc_sum = Einc_sum + i_term*n_term*minus_term*(m_odd_1 + 1j*n_even_1);    
    Es_sum = (Es_sum + i_term*n_term*minus_term*(a_r_n*m_odd_3 + 1j*b_r_n*n_even_3));
end

% Convert from spherical to cartesian co-ordinates
% Einc
Einc_mie = E0inc * Einc_sum;

Einc_mie_x = Einc_mie(:,1) .* sin(theta_selection) .* cos(phi_selection) + ...
             Einc_mie(:,2) .* cos(theta_selection) .* cos(phi_selection) - ...
             Einc_mie(:,3) .* sin(phi_selection);
         
Einc_mie_y = Einc_mie(:,1) .* sin(theta_selection) .* sin(phi_selection) + ...
             Einc_mie(:,2) .* cos(theta_selection) .* sin(phi_selection) + ...
             Einc_mie(:,3) .* cos(phi_selection);
         
Einc_mie_z = Einc_mie(:,1) .* cos(theta_selection) - ...
             Einc_mie(:,2) .* sin(theta_selection);
         
Einc_mie = [Einc_mie_x, Einc_mie_y, Einc_mie_z];

% Es
Es_mie = E0inc * Es_sum;

Es_mie_x = Es_mie(:,1) .* sin(theta_selection) .* cos(phi_selection) + ...
             Es_mie(:,2) .* cos(theta_selection) .* cos(phi_selection) - ...
             Es_mie(:,3) .* sin(phi_selection);
         
Es_mie_y = Es_mie(:,1) .* sin(theta_selection) .* sin(phi_selection) + ...
             Es_mie(:,2) .* cos(theta_selection) .* sin(phi_selection) + ...
             Es_mie(:,3) .* cos(phi_selection);
         
Es_mie_z = Es_mie(:,1) .* cos(theta_selection) - ...
             Es_mie(:,2) .* sin(theta_selection);
         
Es_mie = [Es_mie_x, Es_mie_y, Es_mie_z];

% Et
Et_mie = Einc_vefie + Es_mie;