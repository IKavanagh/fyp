% Create building shape
tic;

% General Information %
% Define information needed for discretising shape
centre = 0.0 + 1i*0.0;
start_point = centre - (x_side * 0.5 + y_side * 0.5 * 1i); % Lower left

% Determine N, multiple of 4
N = floor(x_side / (abs(lambdad) / disc_per_lambda));
while (mod(N, 4) ~= 0)
    N = N + 1;
end
delta_x = x_side / N;

% Determine M, multiple of 4
M = floor(y_side / (abs(lambdad) / disc_per_lambda));
while (mod(M, 4) ~= 0)
    M = M + 1;
end
delta_y = y_side / M;

problem_size = N*M;

% Determine radius (a) of equivalent circle for discretised area
equiv_a = sqrt(delta_x * delta_y / pi);

% Determine dx and dy for walls and doors
w_dx = round(wall / delta_x);
w_dy = round(wall / delta_y);

d_dx = round(door / delta_x);
d_dy = round(door / delta_y);

% Define Building %
create_building

% Define positions %
position = zeros(problem_size, 1);
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
        
        % Define wave numbers
        if (building(y, x) == 1) % Concrete
            wave_number(position_counter, 1) = kc;
        end
        
        if (building(y, x) == 2) % Glass
            wave_number(position_counter, 1) = kg;
        end
        
        if (building(y, x) == 3) % Wood
            wave_number(position_counter, 1) = kw;
        end
        
    end
end

[phi, rho] = cart2pol(real(position), imag(position));

t = toc;
fprintf(1, 'Created building of length %.0fm and width %.0fm with %.0f discretisations in %.4fs\n', x_side, y_side, problem_size, t);

if 1 == 1
    h = figure;
    
    surf(real(pos), imag(pos), building);
    shading interp;
    view(2);
    xlabel('Distance (m)');
    ylabel('Distance (m)');
    title('Plot of Building');
    half_x = x_side * 0.5;
    xlim([-half_x half_x]);
    half_y = y_side * 0.5;
    ylim([-half_y half_y]);
end