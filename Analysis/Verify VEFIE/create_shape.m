% Create shape (area with circle inside, centered on the origin)
% Define vectors to store position in cartesion and polar coordinates
tic;

position = zeros(problem_size, 1);

wave_number(1:problem_size, 1) = k0;

position_counter = 0;

for y = 1:M % y-direction
    for x = 1:N % x-direction
        % Position counter traverses along x-axis
        position_counter = (y - 1)*N + x;
        
        % Define position as center of cell i.e. (x - 1/2)
        position(position_counter, 1) = start_point + ((x - 0.5)*delta_x + (y - 0.5)*delta_y*1i);
        
        % If cell is within circle its wave number is different to free
        % space
        if (abs(position(position_counter, 1) - centre) <= radius)
           wave_number(position_counter, 1) = kd; 
        end
    end
end

[phi, rho] = cart2pol(real(position), imag(position));

t = toc;
fprintf(1, 'Created shape of width %.4f metres and length %.4f meters and discretised into %.0f cells in %.4f seconds\n', length_x_side, length_y_side, problem_size, t);

% Flatten out scatterer so it can be viewed visually
scatterer = zeros(M, N);
for y = 1:M
    for x = 1:N
        position_counter = (y - 1)*N + x;

        if (wave_number(position_counter) ~= k0)
            scatterer(y, x) = 1;
        end
    end
end