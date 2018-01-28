% Create shape (square/rectangle of free space enclosing dielectric
% square/rectange)
tic;

shape = zeros(M, N);

shape_slab % If not using need to comment out

antenna_location = start_point + (x_side * 0.1 + y_side * 0.9 * 1i);

% Define Shape
position = zeros(problem_size, 1);
scatterer = zeros(M, N);
pos = zeros(M, N);

wave_number(1:problem_size, 1) = k0;

position_counter = 0;

wall_start = 100000.0; 

for y = 1:M % y-direction
    for x = 1:N % x-direction
        % Position counter traverses along x-axis
        position_counter = (y - 1)*N + x;
        
        % Define position as center of cell i.e. (x - 1/2)
        position(position_counter, 1) = start_point + ((x - 0.5)*delta_x + (y - 0.5)*delta_y*1i);
        pos(y, x) = position(position_counter);
        
        if (shape(y, x) == 1)
           wave_number(position_counter, 1) = k_d; 
           scatterer(y, x) = 1;
           
           if(real(position(position_counter)) < wall_start)  
                wall_start = real(position(position_counter)); 
           end
        end
        
        % Choose shape
%         shape_block
%         shape_cylinder
    end
end

[phi, rho] = cart2pol(real(position), imag(position));

t = toc;
fprintf(1, 'Created shape of length %.0fm and width %.0fm with %.0f discretisations in %.4fs\n', ...
        x_side, y_side, problem_size, t);
    
active = sum(sum(scatterer));

fprintf(1, 'Scatterer unknowns %d, free space unknowns %d, total unknowns %d\n', active, (problem_size - active), problem_size);


if 1 == 1
    h = figure;
    
    surf(real(pos), imag(pos), shape);
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