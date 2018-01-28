% Create shape (square/rectangle of free space enclosing dielectric
% square/rectange)

x_side_block = 0.5*x_side;
y_side_block = 0.5*y_side;

begin_location = start_point + (x_side_block * 0.5 + y_side_block * 0.5 * 1i);
end_location = begin_location + (x_side_block + y_side_block * 1i);

% If cell is within dielectric its wave number is different to free
% space
if real(begin_location) <= real(position(position_counter)) && ...
     real(end_location) >= real(position(position_counter)) && ...
   imag(begin_location) <= imag(position(position_counter)) && ...
     imag(end_location) >= imag(position(position_counter))
   wave_number(position_counter, 1) = k_d; 
   scatterer(y, x) = 1;
end