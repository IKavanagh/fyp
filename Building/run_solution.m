% Run VEFIE solution for building
clear all
close all
clc

% Load variables
variables

% Create shape for problem to be solved over
create_shape

% Create VEFIE elements
create_vefie_elements

% Plot results
if 1 == 1 % Total electric field
    % Run VEFIE solution
    vefie_solution

    % Plot results
    max_val = max(max(log10(abs(E_total))));
    min_val = min(min(log10(abs(E_total))));
    
    range = max_val - min_val;
    
    addition = range / 10;

    for y = 1:N
        for x = 1:N

            if (building(y,x) == 1 || building(y,x) == 2 || building(y,x) == 3) 
                building(y,x) = max_val; 
            end

            if (building(y,x) == 0 || building(y,x) == 4)
                building(y,x) = min_val; 
            end

        end
    end
    
    h = figure;
    
    surf(real(pos), imag(pos), log10(abs(E_total)));
    hold on
    %surf(real(pos), imag(pos), building);
    
    shading interp;
    view(2);
    xlabel('Distance (m)');
    ylabel('Distance (m)');
    title('Electric Field Throughout Building');
end

if 1 == 2 % Incident field
    h = figure;
    
    surf(real(pos), imag(pos), log10(abs(E_inc)));
    hold on
    surf(real(pos), imag(pos), building);
    
    shading interp;
    view(2);
    xlabel('Distance (m)');
    ylabel('Distance (m)');
    title('Plot of Incident Electric Field in Building');
end

% Movie
if 1 == 1
    h = figure;
    
    start_time = 0;
    end_time = 10;
    step_time = (end_time - start_time) / 200;
    
    t = start_time:step_time:end_time;
    
    X = real(pos);
    Y = imag(pos);
    
    for i = 1:length(t)
        Et = E_total + 1i*E_total;
        wave = exp(1i*t(i));
        
        Z = Et*wave;

        surf(X, Y, real(Z));
        
        shading interp;
        view(2);
        xlabel('Distance (m)');
        ylabel('Distance (m)');
        title('Plot of Total Electric Field Interacting on Cylindrical Scatterer');
        
        movie_frames(i) = getframe(gcf);
    end
    
    movie2avi(movie_frames, 'building.avi');
end

% Plot convergence rate
if 1 == 2
    tol = 1e-7;

    [~, BICGSTAB_FFT] = bicgstab_fft(V, D, G, N, tol);
    [~, CGNE_FFT_REDUCED] = cgne_fft_reduced(V, D, G, N, tol);
    
%     CGNE_FFT_REDUCED1(1:109) = CGNE_FFT_REDUCED(1:109);
%     CGNE_FFT_REDUCED = CGNE_FFT_REDUCED1;
    
%     BICGSTAB_FFT1(1:115) = BICGSTAB_FFT(1:115);
%     BICGSTAB_FFT = BICGSTAB_FFT1;
    
    max_length = max([length(CGNE_FFT_REDUCED) length(BICGSTAB_FFT)]);
    
    n = 0:max_length-1;

    CGNE_FFT_REDUCED(length(CGNE_FFT_REDUCED):max_length) = 0;
    BICGSTAB_FFT(length(BICGSTAB_FFT):max_length) = 0;
    
    h = figure;
    
    plot(n, log10(CGNE_FFT_REDUCED), 'b-', n, log10(BICGSTAB_FFT), 'r--');
    
    xlabel('Number of Iterations');
    ylabel('Error');
    title('Error Compared to Number of Iterations');
    
    legend('CG-NE with Reduced Operator', 'BiCGSTAB');
    ylim([log10(tol*10) log10(0.5e1)]);
    hgsave(h, 'p_convergence_10m_10m_10disc_6er_1sig');
end