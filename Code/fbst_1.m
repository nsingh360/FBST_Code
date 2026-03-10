function [y_nyq, w_l, alpha] = fbst_1(y, L, B, A_eplitz, Psi, beam_idx, fc,Ws, T, delta, exp_fac)
%% fbst_1 finds coefficients using Chirp-z and the inverse problem using O(n^3)
%   Detailed explanation goes here
% A_eplitz : Toeplitz matrix for the step-2
% Psi : Basis matrix with uniformly sampled basis for testing

%% Adjusted for [-pi/2,pi/2]

%% exp_fac factor is used to go slightly over the bandlimite (\Omega + eta)


fc_tilde = (fc+Ws); % adjustment for spatial aliasing
fc_norm = fc/fc_tilde;
[M,N] = size(y);
%% chirp-z transform
F_transform = zeros(M,2*L+1);
A = exp(-1i*pi*exp_fac);%exp(-1i*2*pi*T*W);
W_czt = exp(-1i*pi*exp_fac/L);
for mm=1:M
    F_transform(mm,:) = czt(y(mm,:).',2*L+1,W_czt,A);
end
X_czt = zeros(2*L+1,2*B+1);
A= 1;
B_list = reshape(linspace(0,2*B,2*B+1), [1,2*B+1]);
for ll=1:2*L+1
    l_dash = ll-1;
    A = exp(1i*pi*(1/B)*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L))*B);
    W_czt = exp(1i*pi*(1/B)*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L)));
    coeff1 = exp(-1i*pi*(1/B)*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L))*(M/2-1/2)*B_list); % adjusting for array center
    coeff2 = exp(1i*pi*(1/B)*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L))*(M/2-1/2)*B); % adjusting for negative angles
    X_czt(ll,:) = czt(F_transform(:,ll),2*B+1,W_czt,A);
    X_czt(ll,:) = coeff2*X_czt(ll,:).*coeff1;
end
w_l = X_czt(:,beam_idx);

%% beamforming inverse problem
alpha = inv(A_eplitz + delta*eye(2*L+1))*w_l;
y_nyq = Psi*alpha;

end