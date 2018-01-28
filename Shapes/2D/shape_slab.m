% Create shape (square/rectangle of free space enclosing dielectric
% slab)

cavities = 0;

wall = min([x_side y_side]) / 15;

% Determine dx and dy for walls
w_dx = round(wall / delta_x);
w_dy = round(wall / delta_y);

% Define helpful sizes
hN = floor(N/2);
tN = floor(N/3);
qN = floor(N/4);
fN = floor(N/5);
sN = floor(N/6);
eN = floor(N/8);
nN = floor(N/9);

hM = floor(M/2);
tM = floor(M/3);
qM = floor(M/4);
fM = floor(M/5);
sM = floor(M/6);
eM = floor(M/8);
nM = floor(M/9);

% Create slab
shape(1:M, hN:hN+2*w_dx) = 1;

if (cavities == 1)
    %Make cavities in walls
    shape(floor(nM/1.5):fM, hN-w_dx+round(eN/2):hN+w_dx-round(eN/2)) = 0;
    shape(fM+floor(nM/1.5):2*fM, hN-w_dx+round(eN/2):hN+w_dx-round(eN/2)) = 0;
    
    shape(hN-floor(nM/2):hN+floor(nM/2), hN-w_dx+round(eN/2):hN+w_dx-round(eN/2)) = 0;
    
    shape(N-floor(nM/1.5):-1:N-fM, hN-w_dx+round(eN/2):hN+w_dx-round(eN/2)) = 0;
    shape(N-(fM+floor(nM/1.5)):-1:N-(2*fM), hN-w_dx+round(eN/2):hN+w_dx-round(eN/2)) = 0;
end