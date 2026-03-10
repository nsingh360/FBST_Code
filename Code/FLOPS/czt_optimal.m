
function Xk = czt_optimal(x, M, X, A)
    % CZT implementation based on Definition
    % Computes the Z-transform at points on a logarithmic spiral contour
    %
    % Inputs:
    %   x  - Input signal of length N
    %   M  - Number of frequency points
    %   A  - Starting point of contour
    %   X  - Ratio defining the logarithmic spiral
    %
    % Output:
    %   Xk - Computed Chirp Z-Transform values
    
    N = length(x);  % Length of input signal
    L = 2^nextpow2(N + M - 1);  % Optimal FFT length for efficient computation

    % Compute spiral contour points
    k = (0:M-1)';
    n = (0:N-1)';
    zk = A .* X.^(-k);  % Chirp-Z evaluation points

    % Compute chirp sequences
    Wn = X.^(n.^2 / 2);   % Input modulation
    Wk = X.^(-k.^2 / 2);  % Output modulation

    % Modulate input sequence
    x_mod = (x .* Wn);  

    % Construct chirp kernel for convolution
    chirp_padded = zeros(L, 1);
    chirp_padded(1:N) = conj(Wn);  
    chirp_padded(L-N+2:L) = conj(Wn(N:-1:2));  % Symmetric extension

    % Compute FFTs
    X_fft = fft(x_mod, L);
    C_fft = fft(chirp_padded, L);
    Y_fft = X_fft .* C_fft;

    % Compute inverse FFT and extract the first M points
    y = ifft(Y_fft, L);
    Xk = y(1:M) .* Wk;
end
