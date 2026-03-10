function [y_nyq, w_l] = fbst_0(y, L, B, A_eplitz_inv, Psi, beam_idx, fc, T)
%% fbst_0 finds coefficients using Chirp-z and assumes the inverse is pre-evaluated O(n^2)
%   Detailed explanation goes here
% A_eplitz : Toeplitz matrix for the step-2
% Psi : Basis matrix with uniformly sampled basis for testing
[M,N] = size(y);
delta = 1e-5;
%% chirp-z transform
F_transform = zeros(M,2*L);
A = -1;%exp(-1i*2*pi*T*W);
W_czt = exp(-1i*pi/L);
tic;
for mm=1:M
    F_transform(mm,:) = czt(y(mm,:),2*L,W_czt,A);
end
X_czt = zeros(2*L,B);
A= 1;
for ll=1:2*L
    l_dash = ll-1;
    % A = exp(-1i*pi*(1 + 2*(W/fc)*(l_dash-L/2)/L));
    W_czt = exp(1i*pi*(1/B)*(1 + (1/(2*fc*T*L))*(l_dash-L)));
    X_czt(ll,:) = czt(F_transform(:,ll),B,W_czt,A);
end
w_l = X_czt(:,beam_idx);

%% beamforming inverse problem
alpha = A_eplitz_inv*w_l;
y_nyq = Psi*alpha;
end