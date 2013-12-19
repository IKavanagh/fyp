function Z = conv_fft_1d(X, Y)
    % Calculates X*Y where X is a row of a mirrored toeplitz matrix Y is a
    % column vector
    % Size of problem is N x N

    N = length(Y);

    % Extend Y to length 2N with zeros and convert to row vector
    Y(N+1:2*N, 1) = 0;
    Y = Y.';

    % Extend X to length 2N with a zero and then a reflection
    X(1, N+1) = 0;
    X(1, N+2:2*N) = fliplr(X(1, 2:N));

    % Take fft of both matrices
    Yfft = fft(Y);
    Xfft = fft(X);

    % Multiply matrices point by point and take ifft
    Zi = ifft(Yfft.*Xfft);

    % Take part of solution we need and convert back to column vector
    Z = Zi(1:N).';
end