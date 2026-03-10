
function y = dvm(N, alpha, z)
% DVM - Efficient Delay Vandermonde Matrix-Vector Product
% Implements Algorithm 3.10 from the Efficient DVM paper
%
% INPUTS:
%   N     - Size of the DVM (must be a power of 2)
%   alpha - Complex node (e.g., exp(-1j*omega*tau))
%   z     - Input vector (Nx1)
%
% OUTPUT:
%   y     - Output vector (Nx1)

if mod(N, 2) ~= 0
    error('N must be a power of 2');
end

if N == 2
    y = [1, alpha; 1, alpha^2] * z;
    return;
end

N1 = N / 2;

% Construct diagonal matrix DN
DN = diag(alpha.^(0:N-1).');

% Scale input
u = DN * z;

% Generate companion matrix coefficients using Algorithm 3.6
w = com(N, alpha);

% Companion matrix CN1
CN1 = diag(ones(N1-1,1), -1);
CN1(end, :) = -w(end:-1:1);

% Compute power CN1^N1 using Algorithm 3.7
C_power = CN1;
for k = 1:log2(N1)
    C_power = C_power * C_power;
end

% Construct C_tilde_N
D_tilde = diag(alpha.^(0:N1-1).');
C_tilde = [eye(N1), C_power; D_tilde, alpha^N1 * C_power * D_tilde];

% Apply transformation
v = C_tilde * u;

% Apply two recursive calls
v1 = dvm(N1, alpha^2, v(1:N1));
v2 = dvm(N1, alpha^2, v(N1+1:end));

% Permute
y = zeros(N,1);
y(1:2:end) = v1;
y(2:2:end) = v2;
end
