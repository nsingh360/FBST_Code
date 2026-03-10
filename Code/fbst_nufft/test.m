M = 1000;
N = 513; % odd
x = rand(M,1);
c = exp(1i*2*pi*x);
eps = 1e-6;

F = simpleNUFFT1d_unitInterval(x, c, N, eps);

k = -(N-1)/2 : (N-1)/2;
F_direct = zeros(1, N);
for idx = 1:N
    F_direct(idx) = sum(c .* exp(1i*2*pi*k(idx)*x));
end

fprintf('Relative max error = %.3e\n', max(abs(F(:) - F_direct(:))) / max(abs(F_direct)));
