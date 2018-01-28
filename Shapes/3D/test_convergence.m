% Run VEFIE solution for 3D sphere embedded in a cube with reduced and non
% reduced BiCGSTAB
clear all
close all

% Load variables
variables

% Create shape
create_shape

% Create VEFIE elements
create_vefie_elements

tol = 1e-7;

y_lim_max = log10(1e1);
y_lim_min = log10(tol*10);

[~, BICGSTAB_FFT] = bicgstab_fft(Einc, X, G, position, N, M, O, delta_x, delta_y, delta_z, k0, tol);
[~, BICGSTAB_FFT_REDUCED] = bicgstab_fft_reduced(Einc, X, G, position, N, M, O, delta_x, delta_y, delta_z, k0, tol);

max_length = max([length(BICGSTAB_FFT) length(BICGSTAB_FFT_REDUCED)]);

n = 0:max_length-1;

BICGSTAB_FFT(length(BICGSTAB_FFT):max_length) = 0;
BICGSTAB_FFT_REDUCED(length(BICGSTAB_FFT_REDUCED):max_length) = 0;

h = figure;

plot(n, log10(BICGSTAB_FFT), n, log10(BICGSTAB_FFT_REDUCED));

xlabel('Number of Iterations');
ylabel('Error');
title('Error Compared to Number of Iterations');

legend('BICGSTAB-FFT', 'BICGSTAB-FFT Reduced');
ylim([y_lim_min y_lim_max]);