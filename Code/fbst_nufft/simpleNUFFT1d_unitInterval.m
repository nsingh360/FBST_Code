function F = simpleNUFFT1d_unitInterval(x, c, N, eps)
% 1D Type-1 NUFFT with Gaussian kernel, for x in [0,1], odd N, variable eps
%
% Inputs:
%   x: M points in [0,1]
%   c: M complex values
%   N: odd number of frequency bins
%   eps: desired accuracy
%
% Output:
%   F: N frequency components, k = -(N-1)/2 : (N-1)/2

    if mod(N,2) == 0
        error('N must be odd');
    end

    oversamp = 2;                   % oversampling factor
    Mgrid = oversamp * N;           % oversampled grid size

    % Map x to [0, 2*pi)
    x_scaled = mod(x,1) * 2*pi;

    % Uniform grid points in [0, 2pi)
    grid = (0:Mgrid-1) * (2*pi / Mgrid);

    % Kernel width (number of points)
    w = ceil(-log(eps));  % kernel half-width in grid points

    % Gaussian std dev for kernel
    sigma = w / (2*pi);

    % Precompute spreading radius in radians
    spread_radius = w * 2*pi / Mgrid;

    % Define Gaussian kernel
    kernel = @(d) exp(-(d./sigma).^2);

    % Initialize oversampled grid
    g = zeros(1, Mgrid);

    % Spread data onto grid with proper normalization
    for j = 1:length(x_scaled)
        % Distances with wrap-around in [-pi, pi)
        dist = grid - x_scaled(j);
        dist = mod(dist + pi, 2*pi) - pi;

        % Find indices within kernel support
        idx = find(abs(dist) <= spread_radius);

        % Compute kernel weights
        weights = kernel(dist(idx));

        % Normalize weights to sum to 1 for each point to conserve energy
        weights = weights / sum(weights);

        % Spread c(j) with normalized weights
        g(idx) = g(idx) + c(j) * weights;
    end

    % Compute FFT of oversampled grid (no fftshift here)
    G = fft(g);

    % Frequency vector (0 to Mgrid-1), shifted to [-Mgrid/2,...,Mgrid/2-1]
    kfreq = mod((0:Mgrid-1) + floor(Mgrid/2), Mgrid) - floor(Mgrid/2);

    % Fourier transform of Gaussian kernel at frequencies kfreq
    kernelFT = @(k) sqrt(pi)*sigma*exp(-(sigma*k*pi/Mgrid).^2);

    % Deconvolve kernel effect
    deconv = kernelFT(kfreq);
    F_oversamp = G ./ deconv;

    % Extract center N frequencies
    center = floor(Mgrid/2) + 1;
    halfN = (N-1)/2;
    idxExtract = (center - halfN):(center + halfN);
    idxExtract = mod(idxExtract - 1, Mgrid) + 1; % wrap indices 1-based

    % Return column vector
    F = F_oversamp(idxExtract).';
end
