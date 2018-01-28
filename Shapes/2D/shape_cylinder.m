% Create shape (area with circle inside, centered on the origin)

radius = min([x_side y_side]) * 0.5 * 0.5;

% If cell is within circle its wave number is different to free space
if (abs(position(position_counter, 1) - centre) <= radius)
    wave_number(position_counter, 1) = k_d; 
    scatterer(y, x) = 1;
end 