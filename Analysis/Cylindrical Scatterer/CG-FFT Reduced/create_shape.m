% Create shape (area with circle inside, centered on the origin)
tic;

position = zeros(problem_size, 1);
scatterer = zeros(M, N);
pos = zeros(M, N);

wave_number(1:problem_size, 1) = k0;

position_counter = 0;

for y = 1:M % y-direction
    for x = 1:N % x-direction
        % Position counter traverses along x-axis
        position_counter = (y - 1)*N + x;
        
        % Define position as center of cell i.e. (x - 1/2)
        position(position_counter, 1) = start_point + ((x - 0.5)*delta_x + (y - 0.5)*delta_y*1i);
        pos(y, x) = position(position_counter);
        
        % If cell is within circle its wave number is different to free
        % space
        if (abs(position(position_counter, 1) - centre) <= radius)
           wave_number(position_counter, 1) = kd; 
           scatterer(y, x) = 1;
        end
    end
end

[phi, rho] = cart2pol(real(position), imag(position));

t = toc;
fprintf(1, 'Created shape of width %.0fm and length %.0fm and discretised into %.0f cells in %.4fs\n', x_side, y_side, problem_size, t);