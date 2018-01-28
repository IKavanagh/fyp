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

y_lim_max = log10(0.5e1);
y_lim_min = log10(tol*10);

if 1 == 2 % Plot convergence rate (all methods) - Small problems
    [~, CGNE_FFT] = cgne_fft(V, D, G, N, tol);
    [~, CGNE_FFT_REDUCED] = cgne_fft_reduced(V, D, G, N, tol);
    [~, BICGSTAB_FFT] = bicgstab_fft(V, D, G, N, tol);
    [~, BICGSTAB_FFT_REDUCED] = bicgstab_fft_reduced(V, D, G, N, tol);
    
    max_length = max([length(CGNE_FFT) length(CGNE_FFT_REDUCED) length(BICGSTAB_FFT) length(BICGSTAB_FFT_REDUCED)]);
    
    n = 0:max_length-1;

    CGNE_FFT(length(CGNE_FFT):max_length) = 0;
    CGNE_FFT_REDUCED(length(CGNE_FFT_REDUCED):max_length) = 0;
    BICGSTAB_FFT(length(BICGSTAB_FFT):max_length) = 0;
    BICGSTAB_FFT_REDUCED(length(BICGSTAB_FFT_REDUCED):max_length) = 0;
    
    h = figure;
    
    plot(n, log10(CGNE_FFT), n, log10(CGNE_FFT_REDUCED), n, log10(BICGSTAB_FFT), n, log10(BICGSTAB_FFT_REDUCED));
    
    xlabel('Number of Iterations');
    ylabel('Error');
    title('Error Compared to Number of Iterations');
    
    legend('CGNE-FFT', 'CGNE-FFT Reduced', 'BICGSTAB-FFT', 'BICGSTAB-FFT Reduced');
    ylim([y_lim_min y_lim_max]);
end

if 1 == 2 % Plot convergence rate (no CGNE-FFT) - Large problems
    [~, CGNE_FFT_REDUCED] = cgne_fft_reduced(V, D, G, N, tol);
    [~, BICGSTAB_FFT] = bicgstab_fft(V, D, G, N, tol);
    [~, BICGSTAB_FFT_REDUCED] = bicgstab_fft_reduced(V, D, G, N, tol);
    
    max_length = max([length(CGNE_FFT_REDUCED) length(BICGSTAB_FFT) length(BICGSTAB_FFT_REDUCED)]);
    
    n = 0:max_length-1;

    CGNE_FFT_REDUCED(length(CGNE_FFT_REDUCED):max_length) = 0;
    BICGSTAB_FFT(length(BICGSTAB_FFT):max_length) = 0;
    BICGSTAB_FFT_REDUCED(length(BICGSTAB_FFT_REDUCED):max_length) = 0;
    
    h = figure;
    
    plot(n, log10(CGNE_FFT_REDUCED), n, log10(BICGSTAB_FFT), n, log10(BICGSTAB_FFT_REDUCED));
    
    xlabel('Number of Iterations');
    ylabel('Error');
    title('Error Compared to Number of Iterations');
    
    legend('CGNE-FFT Reduced', 'BICGSTAB-FFT', 'BICGSTAB-FFT Reduced');
    ylim([y_lim_min y_lim_max]);
end

if 1 == 2 % Plot convergence rate (no CGNE-FFT or CGNE-FFT Reduced) - Very large problems
    [~, BICGSTAB_FFT] = bicgstab_fft(V, D, G, N, tol);
    [~, BICGSTAB_FFT_REDUCED] = bicgstab_fft_reduced(V, D, G, N, tol);
    
    BICGSTAB_FFT1(1:802) = BICGSTAB_FFT(1:802);
    BICGSTAB_FFT = BICGSTAB_FFT1;
    
    BICGSTAB_FFT_REDUCED1(1:813) = BICGSTAB_FFT_REDUCED(1:813);
    BICGSTAB_FFT_REDUCED = BICGSTAB_FFT_REDUCED1;
    
    max_length = max([length(BICGSTAB_FFT) length(BICGSTAB_FFT_REDUCED)]);
    
    n = 0:max_length-1;

    BICGSTAB_FFT(length(BICGSTAB_FFT):max_length) = 0;
    BICGSTAB_FFT_REDUCED(length(BICGSTAB_FFT_REDUCED):max_length) = 0;
    
    h = figure;
    
    plot(n, log10(BICGSTAB_FFT), 'b-', n, log10(BICGSTAB_FFT_REDUCED), 'r--');
    
    xlabel('Number of Iterations');
    ylabel('Error');
%     str = sprintf('Error Based on Relative Error\nScatter Unknowns %.2f%%, Free Space Unknowns %.2f%%, Total Unknowns %d', ...
%                    active / problem_size, (problem_size - active) / problem_size, problem_size);
%     % title('Error Based on Relative Error');
%     title(str);
    
    legend('Without Reduced Operator', 'With Reduced Operator');
    ylim([y_lim_min y_lim_max]);
    xlim([0 850]);
end

if 1 == 2 % Plot convergence rate (Non-reduced methods)
    [~, CGNE_FFT] = cgne_fft(V, D, G, N, tol);
    [~, BICGSTAB_FFT] = bicgstab_fft(V, D, G, N, tol);
    
    max_length = max([length(CGNE_FFT) length(BICGSTAB_FFT)]);
    
    n = 0:max_length-1;

    CGNE_FFT(length(CGNE_FFT):max_length) = 0;
    BICGSTAB_FFT(length(BICGSTAB_FFT):max_length) = 0;
    
    h = figure;
    
    plot(n, log10(CGNE_FFT), n, log10(BICGSTAB_FFT));
    
    xlabel('Number of Iterations');
    ylabel('Error');
    title('Error Compared to Number of Iterations');
    
    legend('CGNE-FFT', 'BICGSTAB-FFT');
    ylim([y_lim_min y_lim_max]);
end

if 1 == 2 % Plot convergence rate (Reduced methods)
    [~, CGNE_FFT_REDUCED] = cgne_fft_reduced(V, D, G, N, tol);
    [~, BICGSTAB_FFT_REDUCED] = bicgstab_fft_reduced(V, D, G, N, tol);
    
    max_length = max([length(CGNE_FFT_REDUCED) length(BICGSTAB_FFT_REDUCED)]);
    
    n = 0:max_length-1;

    CGNE_FFT_REDUCED(length(CGNE_FFT_REDUCED):max_length) = 0;
    BICGSTAB_FFT_REDUCED(length(BICGSTAB_FFT_REDUCED):max_length) = 0;
    
    h = figure;
    
    plot(n, log10(CGNE_FFT_REDUCED), n, log10(BICGSTAB_FFT_REDUCED));
    
    xlabel('Number of Iterations');
    ylabel('Error');
    title('Error Compared to Number of Iterations');
    
    legend('CGNE-FFT Reduced', 'BICGSTAB-FFT Reduced');
    ylim([y_lim_min y_lim_max]);
end

if 1 == 1 % Plot convergence rate (CGNE Reduced v BICGSTAB)
    tol = 1e-4;
    
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
    ylim([y_lim_min y_lim_max]);
end

if 1 == 2 % Plot convergence rate (CGNE Non Reduced v CGNE Reduced)
    [~, CGNE_FFT] = cgne_fft(V, D, G, N, tol);
    [~, CGNE_FFT_REDUCED] = cgne_fft_reduced(V, D, G, N, tol);
    
    max_length = max([length(CGNE_FFT) length(CGNE_FFT_REDUCED)]);
    
    n = 0:max_length-1;

    CGNE_FFT(length(CGNE_FFT):max_length) = 0;
    CGNE_FFT_REDUCED(length(CGNE_FFT_REDUCED):max_length) = 0;
    
    h = figure;
    
    plot(n, log10(CGNE_FFT), 'b-', n, log10(CGNE_FFT_REDUCED), 'r--');
    
    xlabel('Number of Iterations');
    ylabel('Error');
    title('Analysis of reduced forward operator applied to CG-NE');
    
    legend('Without Reduced Operator', 'With Reduced Operator');
    ylim([y_lim_min y_lim_max]);
end

if 1 == 2 % Plot convergence rate (no CGNE-FFT or CGNE-FFT Reduced) - Very large problems
    [~, NON_LOSSY] = bicgstab_fft(V, D, G, N, tol);

    % Load variables
    variables_lossy
    
    % Create shape for problem to be solved over
    create_shape

    % Create VEFIE elements
    create_vefie_elements
    
    [~, LOSSY] = bicgstab_fft(V, D, G, N, tol);
    
%     NON_LOSSY1(1:109) = NON_LOSSY(1:109);
%     NON_LOSSY = NON_LOSSY1;
%     
    LOSSY1(1:153) = LOSSY(1:153);
    LOSSY = LOSSY1;
    
    max_length = max([length(NON_LOSSY) length(LOSSY)]);
    
    n = 0:max_length-1;

    NON_LOSSY(length(NON_LOSSY):max_length) = 0;
    LOSSY(length(LOSSY):max_length) = 0;
    
    h = figure;
    
    plot(n, log10(NON_LOSSY), 'b-', n, log10(LOSSY), 'r--');
    
    xlabel('Number of Iterations');
    ylabel('Error');
%     str = sprintf('Error Based on Relative Error\nScatter Unknowns %.2f%%, Free Space Unknowns %.2f%%, Total Unknowns %d', ...
%                    active / problem_size, (problem_size - active) / problem_size, problem_size);
%     % title('Error Based on Relative Error');
%     title(str);
    
    legend('Non lossy', 'Lossy');
    ylim([y_lim_min y_lim_max]);
end