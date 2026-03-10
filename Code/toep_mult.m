function y = toep_mult(T, x)
    % Efficient multiplication of Toeplitz matrix T with vector x
    %
    % T: Toeplitz matrix
    % x: vector to be multiplied
    
    % Get the first column and first row of the Toeplitz matrix
    c = T(:, 1);
    r = T(1, :);
    
    % Use FFT for efficient Toeplitz matrix-vector multiplication
    n = length(x);
    p = [c; 0; flipud(conj(r(2:end)'))];
    y = ifft(fft(p) .* fft([x; zeros(n, 1)]));
    y = y(1:n);
end