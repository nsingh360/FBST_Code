function w = com(N, alpha)
% Compute coefficients of polynomial p(z) = (z - 1)(z - alpha^2)...(z - alpha^(N-2))

N1 = N / 2;
w = [1; zeros(N1-1, 1)];

for k = 1:N1-1
    alpha_pow = alpha^(2*k);
    w_new = zeros(N1, 1);
    w_new(1) = -alpha_pow * w(1);
    for j = 2:k+1
        w_new(j) = w(j-1) - alpha_pow * w(j);
    end
    w = w_new;
end

end
