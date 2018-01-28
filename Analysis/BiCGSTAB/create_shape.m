% Create shape (area with circle inside, centered on the origin)
tic;

antenna_location = start_point + (x_side * 0.1 + y_side * 0.9 * 1i);
        
radius = min([x_side y_side]) * 0.5 * 0.5;

% Define Shape
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

        % If cell is within circle its wave number is different to free space
        if (abs(position(position_counter, 1) - centre) <= radius)
            wave_number(position_counter, 1) = k_d; 
            scatterer(y, x) = 1;
        end 
    end
end

[phi, rho] = cart2pol(real(position), imag(position));

t = toc;
fprintf(1, 'Created shape of length %.0fm and width %.0fm with %.0f discretisations in %.4fs\n', ...
        x_side, y_side, problem_size, t);
    
active = sum(sum(scatterer));

fprintf(1, 'Scatterer unknowns %d, free space unknowns %d, total unknowns %d\n', active, (problem_size - active), problem_size);