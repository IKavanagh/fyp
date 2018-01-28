% Create shape (area with circle inside, centered on the origin)
tic;

position = zeros(problem_size, 3); % x, y , z

% Spherical co-ordinates
rho = zeros(problem_size, 1); % sqrt(x^2 + y^2 + z^2)
phi = zeros(problem_size, 1); % atan2(y, x);
r_xy = zeros(problem_size, 1); % sqrt(x^2 + y^2)
theta = zeros(problem_size, 1); % atan2(r_xy, z);

wave_number(1:problem_size, 1) = k0;

position_counter = 0;

for z = 1:O % z-direction
    for y = 1:M % y-direction
        for x = 1:N % x-direction
            % Position counter traverses along x-axis and then y-axis
            position_counter = (z - 1)*N*M + (y - 1)*N + x;
            
            % Define position as centre of each cell
            position(position_counter, 1:3) = start_point + [(x - 0.5)*delta_x, (y - 0.5)*delta_y, (z - 0.5)*delta_z];
            
            rho(position_counter) = sqrt(sum(position(position_counter, 1:3).^2));
            phi(position_counter) = atan2(position(position_counter, 2), position(position_counter, 1));
            r_xy(position_counter) = sqrt(sum(position(position_counter, 1:2).^2));
            theta(position_counter) = atan2(r_xy(position_counter), position(position_counter,3));
            
            % Check cell is within sphere
            if (rho(position_counter) < radius)
                % Scatterer can't be on boundary for Eirr
                if (sum(position(position_counter, 1:3) > (start_point + [delta_x, delta_y, delta_z])) == 3) && ...
                    (sum(position(position_counter, 1:3)  < (start_point + [x_side, y_side, z_side] - [delta_x, delta_y, delta_z])) == 3)
                    wave_number(position_counter) = k_d;
                end
            end
        end
    end
end

t = toc;
fprintf(1, 'Created shape of width %.2fm, length %.2fm and height %.2fm and discretised into %.0f cells in %.4fs\n', x_side, y_side, z_side, problem_size, t);