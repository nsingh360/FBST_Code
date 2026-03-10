function [U, C] = prolate_approx(t_nu,W,T,L,interp_factor)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
N = length(t_nu); % number of non-uniform samples
N_interp = N*interp_factor;
t_u = [0:N_interp]*N/N_interp;
K = ceil(2*W*T) + L; % truncation limit
[V,E] = dpss(N_interp+1,W*T,K); % generate uniformly sampled PSWF's
V_nu = interp1(t_u,V,t_nu,'pchip',0); % interpolate to non-uniform points
[U,S,W] = svd(V_nu,'econ'); % factor basis
C = S*W'*diag(E)*W*S';
C = C/trace(C);

end

