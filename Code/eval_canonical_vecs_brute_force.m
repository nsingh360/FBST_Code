function [X1, X2, X3, X4] = eval_canonical_vecs_brute_force(T,N)
e1 = zeros(N,1);
e1(1) = 1;
en = zeros(N,1);
en(end) = 1;
alpha = 0.5;
v_alpha = zeros(N,1);
v_alpha(2:N) = -1*flip(T(1,2:N));
v_alpha(1) = alpha;

b = [e1 en];
gen_vecs = T\b;

u = gen_vecs(:,1);
v = gen_vecs(:,2);
u_c1 = u;
u_r1 = zeros(N,1);
u_r1(1) = u(1);
X1 = toeplitz(u_c1, u_r1);

% v_r1 = ones(N,1);
v_r1 = flip(v(1:N));
v_c1 = zeros(N,1);
v_c1(1) = v_r1(1);
X2 = toeplitz(v_c1, v_r1);

v_c1 = v;
v_r1 = zeros(N,1);
v_r1(1) = v(1);
X3 = toeplitz(v_c1, v_r1);

% u_r1 = zeros(N,1);
u_r1 = flip(u(1:N));
u_c1 = zeros(N,1);
u_c1(1) = u_r1(1);
X4 = toeplitz(u_c1, u_r1);
end