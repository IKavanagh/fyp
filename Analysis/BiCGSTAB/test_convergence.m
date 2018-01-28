% Run VEFIE solution for scatterer
clear all
close all

% Load variables
variables

% Create shape for problem to be solved over
create_shape

% Create VEFIE elements
create_vefie_elements

tol = 1e-7;

y_lim_max = log10(1e0);
y_lim_min = log10(tol*10);

if 1 == 1 % Plot convergence rate with BiCGSTAB
    [~, BICGSTAB_FFT] = bicgstab_fft(V, D, G, N, tol);
    [~, BICGSTAB_FFT_REDUCED] = bicgstab_fft_reduced(V, D, G, N, tol);
    
    % Create VEFIE elements for scatterer only
    create_vefie_elements_scatterer
    
    [~, BICGSTAB] = bicgstab(V, Z, tol);
    
    max_length = max([length(BICGSTAB) length(BICGSTAB_FFT) length(BICGSTAB_FFT_REDUCED)]);
    
    n = 0:max_length-1;

    BICGSTAB(length(BICGSTAB):max_length) = 0;
    BICGSTAB_FFT(length(BICGSTAB_FFT):max_length) = 0;
    BICGSTAB_FFT_REDUCED(length(BICGSTAB_FFT_REDUCED):max_length) = 0;
    
    n = 0:29;
    
    BICGSTAB1(1:30) = BICGSTAB(1:30);
    BICGSTAB_FFT1(1:30) = BICGSTAB_FFT(1:30);
    BICGSTAB_FFT_REDUCED1(1:30) = BICGSTAB_FFT_REDUCED(1:30);
    
    h = figure;
    
    plot(n, log10(BICGSTAB1), 'b-', n, log10(BICGSTAB_FFT_REDUCED1), 'r--');
    
    xlabel('Number of Iterations');
    ylabel('Error');
    title('Iteration Error');
    
    legend('Reduced Problem', 'Original Problem with Reduced Operator');
    ylim([y_lim_min y_lim_max]);
    xlim([0 30]);
end

if 1 == 2 % Plot convergence rate without BiCGSTAB
    [~, BICGSTAB_FFT] = bicgstab_fft(V, D, G, N, tol);
    [~, BICGSTAB_FFT_REDUCED] = bicgstab_fft_reduced(V, D, G, N, tol);
    
    max_length = max([length(BICGSTAB_FFT) length(BICGSTAB_FFT_REDUCED)]);
    
    n = 0:max_length-1;

    BICGSTAB_FFT(length(BICGSTAB_FFT):max_length) = 0;
    BICGSTAB_FFT_REDUCED(length(BICGSTAB_FFT_REDUCED):max_length) = 0;
    
    h = figure;
    
    plot(n, log10(BICGSTAB_FFT), n, log10(BICGSTAB_FFT_REDUCED));
    
    xlabel('Number of Iterations');
    ylabel('Error');
    title('Iteration Error');
    
    legend('BICGSTAB', 'BICGSTAB-FFT');
    ylim([y_lim_min y_lim_max]);
end