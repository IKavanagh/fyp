% Define Building %
% Specify building materials
free_space = 0;
concrete = 1;
glass = 2;
wood = 3;

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

% Create complex building
building = zeros(M, N);

% Make all outer walls concrete
building(1:1+w_dy, 1:N) = 1;
building(M:-1:M-w_dy, 1:N) = 1;
building(1:M, 1:1+w_dx) = 1;
building(1:M, N:-1:N-w_dx) = 1;

% Create 2 column layout
building(1:M, hN:hN+w_dx) = 1;

% Create 2 rooms in left column and 3 in right
% Left column
building(hM:hM+w_dy, 1:hN) = 1;
% Right column
building(tM:tM+w_dy, hN:N) = 1;
building(2*tM:2*tM+w_dy, hN:N) = 1;

% Add metal cabinets
% building(hM+w_dy+1:hM+nM+w_dy, hN-nN:hN-1) = 4;

% Add doors and gaps
% Top left room, bottom right room
building(hM:hM+w_dy, nN:nN+d_dx) = 3;
% Top left room, top right room
building(M-eM:-1:M-eM-d_dy, hN:hN+w_dx) = 3;
% Bottom left room, bottom right room
building(round(sM/2):tM-round(sM/2), hN:hN+w_dx) = 3;
% Bottom left room, middle right room
building(tM+2*w_dx:hM-w_dx, hN:hN+w_dx) = 3;
% Top right room, middle right room
building(2*tM:2*tM+w_dy, hN+sN:hN+sN+d_dx) = 3;
% Middle right room, bottom right room
building(tM:tM+w_dy, hN+nN:N-nN) = 0;

% Add windows
% Top right
building(M:-1:M-w_dy, hN+sN:hN+sN+d_dx) = 2;
% Top left
building(M:-1:M-w_dy, eN:hN-eN) = 2;
% Bottom left
building(eM:hM-eM, 1:1+w_dx) = 2;
building(1:1+w_dy, sN:tN) = 2;
% Bottom right
building(1:1+w_dy, hN+eN:N-eN) = 2;

% Add Obstacles
% Pillars
% pillar_centre(1) = (-0.25 -0.25i)*x_side ;   
% pillar_radius(1) = 0.03*x_side;
% pillar_centre(2) = (0.25 - 0.25i)*x_side ;   
% pillar_radius(2) = 0.03*x_side;
% pillar_centre(3) = (-0.35 +0.25i)*x_side ;   
% pillar_radius(3) = 0.025*x_side;
% pillar_centre(4) = (0.35 +0.25i)*x_side ;    
% pillar_radius(4) = 0.025*x_side;
% pillar_centre(5) = (-0.05 +0.2i)*x_side ;    
% pillar_radius(5) = 0.035*x_side;
% 
% for ct1 = 1:N
%     for ct2 = 1:N     
%         
%         location = start_point + ((ct1 - 0.5)*delta_x + (ct2 - 0.5)*delta_y*1i);
%         
%         for ct3 = 1:5
%             distance = abs(location  - pillar_centre(ct3)); 
%             if(distance < pillar_radius(ct3))
%                 building(ct1,ct2) = 1; 
%             end
%         end
%     end
% end