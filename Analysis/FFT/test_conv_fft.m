% Test whether matrix multiplication using convolution and FFT is correct
clear all
close all

% Load variables
variables

% Create shape for problem to be solved over
create_shape

% Create VEFIE elements
create_vefie_elements

% Test 1D problem (i.e. multiplication of a toeplitz matrix and column
% vector)
% Extract 1D problem (N x N from G and 1 x 2N from V)
Gm = G(1:N, 1:N);
Vm = V(1:N, 1);

% Solve by matrix multiplication
Am = Gm*Vm;

% Solve by convolution and fft
Ac = conv_fft_1d(G(1, 1:N), V(1:N, 1));

A = Am - Ac;
if (orderof(max(abs(A))) <= -10)
    fprintf(1, '1D Multiplication is correct\n');
else
    fprintf(1, '1D Multiplication is wrong\n');
end

% Test 2D problem (i.e. multiplication of a block toeplitz mirrored matrix
% and column vector)
% Solve exactly

% Create D and G vectors
D_vec = zeros(problem_size, 1);

for ct = 1:problem_size
    D_vec(ct) = D(ct, ct);
end

G_vec = G(1, :);

Am = G*D_vec;

Ac = conv_fft_2d(G_vec, D_vec, N);

A = Am - Ac;
if (orderof(max(abs(A))) <= -10)
    fprintf(1, '2D Multiplication is correct\n');
else
    fprintf(1, '2D Multiplication is wrong\n');
end