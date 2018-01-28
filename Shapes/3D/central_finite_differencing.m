function Eirr = central_finite_differencing(A, position, N, M, O, dx, dy, dz)
    % Calculates the irrotational part of the electric field based on
    % Eirr = grad*div*A, where A = dv * convolution(G*D*E)

    start_point = position(1,:) - [0.5*dx, 0.5*dy, 0.5*dz];
    x_side = N*dx;
    y_side = M*dy;
    z_side = O*dz;
    
    problem_size = length(A(:, 1));
    
    Eirr = zeros(problem_size,3);
    
    for n = 1:problem_size
        % Eirr can't be calculated on the boundary
        if (sum(position(n, 1:3) > (start_point + [dx, dy, dz])) == 3) && ...
           (sum(position(n, 1:3)  < (start_point + [x_side, y_side, z_side] - [dx, dy, dz])) == 3)
            
        
            Eirr(n, 1) = (4 / (dx*dx))*(A(n-1,1)-2*A(n,1)+A(n+1,1)) + ... % Factor of 4 to cancel out division by 4 later on
                         (1 / (dx*dy))*(A(n-1-N,2)-A(n-1+N,2)-A(n+1-N,2)+A(n+1+N,2)) + ...
                         (1 / (dx*dz))*(A(n-1-N*M,3)-A(n-1+N*M,3)-A(n+1-N*M,3)+A(n+1+N*M,3));


            Eirr(n, 2) = (4 / (dy*dy))*(A(n-N,2)-2*A(n,2)+A(n+N,2)) + ... % Factor of 4 to cancel out division by 4 later on
                         (1 / (dx*dy))*(A(n-1-N,1)-A(n-1+N,1)-A(n+1-N,1)+A(n+1+N,1)) + ...
                         (1 / (dy*dz))*(A(n-N-N*M,3)-A(n-N+N*M,3)-A(n+N-N*M,3)+A(n+N+N*M,3));

            Eirr(n, 3) = (4 / (dz*dz))*(A(n-N*M,3)-2*A(n,3)+A(n+N*M,3)) + ... % Factor of 4 to cancel out division by 4 later on
                         (1 / (dx*dz))*(A(n-1-N*M,1)-A(n-1+N*M,1)-A(n+1-N*M,1)+A(n+1+N*M,1)) + ...
                         (1 / (dy*dz))*(A(n-N-N*M,2)-A(n-N+N*M,2)-A(n+N-N*M,2)+A(n+N+N*M,2));
        end
        
    end
    Eirr = (1 / 4) * Eirr; % Multiplication of jw and division by kb^2 not needed as A is calculated without this factor present
end