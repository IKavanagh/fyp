function Z = conv_fft_2d(X, Y, N)
    % Calculates X*Y where X is a row of a mirrored block toeplitz matrix
    % in the N direction and Y is a column vector
    % Size of problem is N x M
    
    M = length(Y) / N;
    
    % Reshape Y into M x N matrix
    Y = reshape(Y, N, M).';
    
    % Embed Y into a 2M x 2N matrix padded with zeros
    Y(M+1:2*M, 1:2*N) = 0;
    
    % Reshape X into M x N matrix
    X = reshape(X, N, M).';
    
    % Embed X into a 2M x 2N matrix of carefully selected values
    X(M+1, N+1) = 0;
    X(:, N+2:2*N) = fliplr(X(:, 2:N));
    X(M+2:2*M, :) = flipud(X(2:M, :));
    
    % Take fft of both matrices
    Yfft = fft2(Y);
    Xfft = fft2(X);

    % Multiply matrices point by point and take ifft
    Zi = ifft2(Yfft.*Xfft);
    
    % Take part of solution we need and convert back to column vector
    Z = reshape(Zi(1:M, 1:N).', M*N, 1);
end