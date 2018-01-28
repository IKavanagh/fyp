function [A] = conv_fft_3D(G, E, N, M, O)

    % Calculates A = convolution(G*D*E) where G is a row of a 3D Toeplitz
    % matrix in the N direction and E = D*E is a 3D column vector, N are the
    % number of steps in the x direction, M the number of steps in the y
    % direction and O the number of steps in the z direction
    
    % G needs to be reshaped into a 2N x 2M x 2O matrix and padded with
    % carefully selected values just as for 2D case
    G = reshape(G, N, M, O);
    
    for counter = 1:O
        G(:, M+2:2*M, counter) = fliplr(G(:, 2:M, counter));
        G(N+2:2*N, :, counter) = flipud(G(2:N, :, counter));
    end
    G(:, :, O+2:2*O) = G(:, :, O:-1:2);

    % E needs to be reshaped into a 2N x 2M x 2O matrix for x, y and z
    % co-ordinates and padded with zeros
    Et = reshape(E(:, 1), N, M, O);
    Ex = zeros(2*N, 2*M, 2*O);
    Ex(1:N, 1:M, 1:O) = Et;
    
    Et = reshape(E(:, 2), N, M, O);
    Ey = zeros(2*N, 2*M, 2*O);
    Ey(1:N, 1:M, 1:O) = Et;
    
    Et = reshape(E(:, 3), N, M, O);
    Ez = zeros(2*N, 2*M, 2*O);
    Ez(1:N, 1:M, 1:O) = Et;

    % Take 3D fft of all matrices
    Gfft = fftn(G);
    Exfft = fftn(Ex);
    Eyfft = fftn(Ey);
    Ezfft = fftn(Ez);
    
    % Multiply Gfft by E(x, y, z) point by point and take 3D ifft
    Zix = ifftn(Gfft.*Exfft);
    Ziy = ifftn(Gfft.*Eyfft);
    Ziz = ifftn(Gfft.*Ezfft);

    Ax = zeros(N*M*O,1);
    Ay = zeros(N*M*O,1); 
    Az = zeros(N*M*O,1);
    
    Ax(:) = Zix(1:N,1:M,1:O);
    Ay(:) = Ziy(1:N,1:M,1:O);
    Az(:) = Ziz(1:N,1:M,1:O);
    
    A = [Ax, Ay, Az];
end