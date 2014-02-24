% Test convergence of CGNE, CGNE-FFT and CGNE-FFT Reduced
clear all
close all

% Define increasing scale for x_side and y_side
start_point = 0;
end_point = 10;
step_point = 1;

x_side_vec = start_point:step_point:end_point;

% Define vectors to hold problem_size, time and number of iterations
l = length(x_side_vec);

size = zeros(1, l);
time = zeros(2, l);
iter = zeros(2, l);

for i = 2:l
    % Define x and y for this iteration
    x_side = x_side_vec(i);
    y_side = x_side;
    
    % Load variables
    variables
    
    % Create shape
    create_shape
    
    size(i) = problem_size;
    
    % CGNE-FFT %
    % Begin timer
    tic;
    
    % Create vefie elements for CGNE-FFT
    create_vefie_elements_fft
    
    % Solve by CGNE-FFT
    solve_cgne_fft
    
    % Stop timer and store number of iterations
    t = toc;
    
    time(1, i) = t;
    iter(1, i) = n;
    
    % CGNE-FFT Reduced %
    % Begin timer
    tic;
    
    % Create vefie elements for CGNE-FFT Reduced
    create_vefie_elements_fft
    
    % Solve by CGNE-FFT Reduced
    solve_cgne_fft_red
    
    % Stop timer and store number of iterations
    t = toc;
    
    time(2, i) = t;
    iter(2, i) = n;
end

% Plot results
if 1 == 1
    h = figure;
    
    plot(size, time(1, :), size, time(2, :));
    
    legend('CGNE-FFT', 'CGNE-FFT Reduced');
    title('Time Take to Solve VEFIE Problem');
    xlabel('Size of Problem');
    ylabel('Time (s)');
    xlim([min(size) max(size)]);
    
    h = figure;
    
    plot(size, iter(1, :), size, iter(2, :));
    
    legend('CGNE-FFT', 'CGNE-FFT Reduced');
    title('Number of Iterations Taken to Solve VEFIE Problem');
    xlabel('Size of Problem');
    ylabel('Number of Iterations');
    xlim([min(size) max(size)]);
end