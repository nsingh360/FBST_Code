function [A_basis, V_nyq, K] = construct_slepian_basis(W, t, T, L, M, N, T_ap)

MN = M*N;
K = ceil(2*W*T_ap)+6;
V_u = dpss(MN+1,W*T_ap,K); % firt K uniformly sampled PSWF's
t_u = [0:(MN)]'; % uniform sampling points over the interval
t_nu = [M*N*(t-min(t(:)))/T_ap]; % non-uniform sample points, scaled to match interval
V_nu = interp1(t_u,V_u,t_nu,'pchip',0); % non-uniform

ls = reshape(linspace(-L,L-1,2*L),[2*L,1]);
F_L = ls*(1/2)*(1/(L*T));
exp_basis = exp(-1i*2*pi*(F_L*t.'));
A_basis = exp_basis*V_nu;

fs = 1/T;
t_nyq = [0:N-1]'/fs; % nyquist samples for testing
V_nyq = interp1(t_u,V_u,MN*(t_nyq-min(t))/T_ap,'pchip',0); % interpolate to nyquist

end