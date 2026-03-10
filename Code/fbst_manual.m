function [y_nyq, w_l, alpha] = fbst_manual(y, L, B, A_eplitz, Psi, beam_idx, fc,Ws, T, delta, exp_fac)
%% fbst_manual finds coefficients using Chirp-z (implemented using ffts) and the inverse problem using O(n^3)
fc_tilde = (fc+Ws); % adjustment for spatial aliasing
fc_norm = fc/fc_tilde;
[M,N] = size(y);

%% Filters for chirp-z
%----- stage-1 filters
A_stage_1 = exp(-1i*pi*exp_fac);%exp(-1i*2*pi*T*W);
W_czt_stage_1 = exp(-1i*pi*exp_fac/L);
%------ filter design for chirp-z stage-1 (taken for each row of y)
m_len1 = N; % input signal length
k1_1 = 2*L+1; % chirp-z length
n_len = 1;
nfft_1 = 2^nextpow2(m_len1+k1_1-1);
%------- Premultiply data.
kk = ((-m_len1+1):max(k1_1-1,m_len1-1)).';
kk2 = (kk .^ 2) ./ 2;
ww_stage_1 = W_czt_stage_1 .^ (kk2);   % <----- Chirp filter is 1./ww
nn = (0:(m_len1-1))';
aa = A_stage_1 .^ ( -nn );
aa_stage_1 = aa.*ww_stage_1(m_len1+nn);
fv_stage_1 = fft( 1 ./ ww_stage_1(1:(k1_1-1+m_len1)), nfft_1 );   % <----- Chirp filter.

%------ stage-2 filters
m_len2 = M; % input signal length
k1_2 = 2*B+1; % chirp-z length
n_len = 1;
nfft_2 = 2^nextpow2(m_len2+k1_2-1);
%------- Premultiply data.
kk = ((-m_len2+1):max(k1_2-1,m_len2-1)).';
kk2 = (kk .^ 2) ./ 2;
nn = (0:(m_len2-1))';

A_coefs_stage_2 = zeros(2*L+1,1);
W_czt_coefs_stage_2 = zeros(2*L+1,1);
ww_stage2 = zeros(length(kk2),2*L+1);
aa_stage2 = zeros(m_len2,2*L+1);
FV_filters_stage_2 = zeros(nfft_2,2*L+1);

for ll=1:2*L+1
    l_dash = ll-1;
    A_coefs_stage_2(ll,1) = exp(1i*pi*(1/B)*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L))*B);
    aa = A_coefs_stage_2(ll,1) .^ ( -nn );
    W_czt_coefs_stage_2(ll,1) = exp(1i*pi*(1/B)*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L)));
    ww_stage2(:,ll) = W_czt_coefs_stage_2(ll,1) .^ (kk2);   % <----- Chirp filter is 1./ww
    aa_stage2(:,ll) = aa.*ww_stage2(m_len2+nn,ll);
    FV_filters_stage_2(:,ll) = fft( 1 ./ ww_stage2(1:(k1_2-1+m_len2)), nfft_2 );   % <----- Chirp filter.
end

F_transform = zeros(M,2*L+1);
% A = -1;%exp(-1i*2*pi*T*W);
% W_czt = exp(-1i*pi/L);
for mm=1:M
    y_in = y(mm,:).';
    y_in = y_in .* aa_stage_1(:,ones(1,n_len));
    fy = fft(  y_in, nfft_1 );
    fy = fy .* fv_stage_1(:,ones(1, n_len));
    g  = ifft( fy );
    g = g( m_len1:(m_len1+k1_1-1), :, :) .* ww_stage_1( m_len1:(m_len1+k1_1-1),ones(1, n_len) );
    F_transform(mm,:) = g;
end

X_czt = zeros(2*L+1,2*B+1);
% A= 1;
B_list = reshape(linspace(0,2*B,2*B+1), [1,2*B+1]);
for ll=1:2*L+1
    l_dash = ll-1;

    y_in = F_transform(:,ll);
    ww = W_czt_coefs_stage_2(ll,1) .^ (kk2);
    y_in = y_in .* aa_stage2(:,ll);
    fy = fft(  y_in, nfft_2 );
    fv = FV_filters_stage_2(:,ll);
    fy = fy .* fv(:,ones(1, n_len));
    g  = ifft( fy );
    g = g( m_len2:(m_len2+k1_2-1), :, :) .* ww_stage2( m_len2:(m_len2+k1_2-1),ll );
    % A = exp(1i*pi*(1/B)*(fc_norm + (1/(2*fc_tilde*T*L))*(l_dash-L))*B);
    % W_czt = exp(1i*pi*(1/B)*(fc_norm + (1/(2*fc_tilde*T*L))*(l_dash-L)));
    coeff1 = exp(-1i*pi*(1/B)*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L))*(M/2-1/2)*B_list); % adjusting for array center
    coeff2 = exp(1i*pi*(1/B)*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L))*(M/2-1/2)*B); % adjusting for negative angles
    X_czt(ll,:) = g;
    X_czt(ll,:) = coeff2*X_czt(ll,:).*coeff1;
end
w_l = X_czt(:,beam_idx);


% %% chirp-z transform
% F_transform = zeros(M,2*L+1);
% A = exp(-1i*pi*exp_fac);%exp(-1i*2*pi*T*W);
% W_czt = exp(-1i*pi*exp_fac/L);
% %------ filter design for chirp-z stage-1 (taken for each row of y)
% m = N; % input signal length
% k1 = 2*L+1; % chirp-z length
% n = 1;
% nfft = 2^nextpow2(m+k1-1);
% %------- Premultiply data.
% kk = ((-m+1):max(k1-1,m-1)).';
% kk2 = (kk .^ 2) ./ 2;
% ww = W_czt .^ (kk2);   % <----- Chirp filter is 1./ww
% nn = (0:(m-1))';
% aa = A .^ ( -nn );
% aa = aa.*ww(m+nn);
% 
% fv = fft( 1 ./ ww(1:(k1-1+m)), nfft );   % <----- Chirp filter.
% 
% for mm=1:M
%     y_in = y(mm,:).';
%     y_in = y_in .* aa(:,ones(1,n));
%     fy = fft(  y_in, nfft );
%     fy = fy .* fv(:,ones(1, n));
%     g  = ifft( fy );
%     g = g( m:(m+k1-1), :, :) .* ww( m:(m+k1-1),ones(1, n) );
%     F_transform(mm,:) = g;
% end
% 
% %---------- Filter design for stage-2 (taken for each column of F_transform)
% m = M; % input signal length
% k1 = 2*B+1; % chirp-z length
% n = 1;
% nfft = 2^nextpow2(m+k1-1);
% %------- Premultiply data.
% kk = ((-m+1):max(k1-1,m-1)).';
% kk2 = (kk .^ 2) ./ 2;
% nn = (0:(m-1))';
% 
% A_coefs = zeros(2*L+1,1);
% W_czt_coefs = zeros(2*L+1,1);
% FV_filters = zeros(nfft,2*L+1);
% 
% for ll=1:2*L+1
%     l_dash = ll-1;
%     A_coefs(ll,1) = exp(1i*pi*(1/B)*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L))*B);
%     W_czt_coefs(ll,1) = exp(1i*pi*(1/B)*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L)));
%     ww = W_czt_coefs(ll,1) .^ (kk2);   % <----- Chirp filter is 1./ww
%     FV_filters(:,ll) = fft( 1 ./ ww(1:(k1-1+m)), nfft );   % <----- Chirp filter.
% end
% 
% X_czt = zeros(2*L+1,2*B+1);
% B_list = reshape(linspace(0,2*B,2*B+1), [1,2*B+1]);
% for ll=1:2*L+1
%     l_dash = ll-1;
% 
%     y_in = F_transform(:,ll);
%     ww = W_czt_coefs(ll,1) .^ (kk2);
%     aa = A_coefs(ll,1) .^ ( -nn );
%     aa = aa.*ww(m+nn);
%     y_in = y_in .* aa(:,ones(1,n));
%     fy = fft(  y_in, nfft );
%     fv = FV_filters(:,ll);
%     fy = fy .* fv(:,ones(1, n));
%     g  = ifft( fy );
%     g = g( m:(m+k1-1), :, :) .* ww( m:(m+k1-1),ones(1, n) );
%     coeff1 = exp(-1i*pi*(1/B)*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L))*(M/2-1/2)*B_list); % adjusting for array center
%     coeff2 = exp(1i*pi*(1/B)*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L))*(M/2-1/2)*B); % adjusting for negative angles
%     X_czt(ll,:) = g;
%     X_czt(ll,:) = coeff2*X_czt(ll,:).*coeff1;
% end
% w_l = X_czt(:,beam_idx);

%% beamforming inverse problem
alpha = inv(A_eplitz + delta*eye(2*L+1))*w_l;
y_nyq = Psi*alpha;

end