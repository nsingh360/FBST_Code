clear
clc
close all

%% 1) Adjust for spatial aliasing (fc change) (Done)
%% 2) Adjust for negative indices (only for K) (-B to B) (Done)
%% 3) Adjust for negative elevation angle (TODO)

%% NOTE: K is associated with sin(phi) and J with cos(phi) 
%% therefore for sampling phi in [-pi/2,pi/2] J always need to be positive.

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)
exp_fac = 1.05;

M = 2^8; % \sqrt(M) x \sqrt(M) dimensional planar array
B = ceil((-3 + sqrt(9 + 8*(M-1)))/4); % no. of beams
N = 2^6; % no. of samples
L = ceil(N); % no. of samples in frequency is 2*L+1
MN = M*N;

S = log2(M)/2;
R = 2;
%% Angle sampling for final stage
r_lin = linspace(0,R^S-1,R^S);
a1 = 1 - (2*r_lin + 1)/R^S; % sin(theta)cos(phi)
b1 = a1; % sin(theta)sin(phi)
a_idx = 10;
b_idx = 15;

%% Sampling over Theta and Phi
% Theta = zeros(B^2,1);
% Phi = zeros(B^2,1);
J = linspace(0,B,B+1);
K = linspace(-B,B,2*B+1);
fc = 20e9;  % center frequency
c = physconst('LightSpeed');
Ws = 5e9; % to prevent spatial aliasing
fc_tilde = (fc+Ws); % adjustment for spatial aliasing
fc_norm = fc/fc_tilde;
lambda = c/(fc_tilde); % wavelngth
W = 5e9; % bandwidth
fs = 2*W; % sampling frequency
T = 1/fs;
% N = max(32,ceil(2/(2/1.05-1)*(sqrt(2)*(sqrt(M)-1)*W/(fc+W)-1)));
% L = ceil(N); % no. of samples in frequency is 2*L+1
% MN = M*N;
freq_resolution = fs/N; % frequency grid for sub-band beamformign

%% select a random theta, phi for beamforming
j_idx = 3;%B+1;%3;%randsample(1:B+1,1);
k_idx = 21;%B+1;%21;%randsample(1:2*B+1,1);
theta = asin(sqrt(J(j_idx)^2 + K(k_idx)^2)/(B));
phi = asin(K(k_idx)/(sqrt(J(j_idx)^2 + K(k_idx)^2)));
% phi = atan(b1(b_idx)/a1(a_idx));
% theta = asin(sqrt(a1(a_idx)^2 + b1(b_idx)^2));

Tau = zeros(sqrt(M),sqrt(M));
m = (-sqrt(M)/2+1/2):(sqrt(M)/2-1/2);


[X_mesh, Y_mesh] = meshgrid(m,m);
x_pos = [X_mesh(:),Y_mesh(:)]*lambda/2; % sensor positions
u_s = sin(theta)*[cos(phi);sin(phi)]; % normal vector for determining delays
tau = x_pos*u_s/c; % relative delays to phase center
eta = lambda/(2*c);

%% Signal specs that do not need to be redifined in loop
n_sinusoids = 100; % number of sinusoids in signal
n = [0:N-1]/fs; % temporal sample vectors
t_array = n-tau; % delays across array and time
t = t_array(:); %
T_aperture = max(t)-min(t);% temporal aperture of the array

demod_phase = repmat(exp(-1i*2*pi*fc*tau),N,1);
mod_phase = exp(1i*2*pi*fc*tau);

taps = 8;
edge_lim = taps + ceil(fs*(max(tau)-min(tau))/2);% clip signal to account for edge distortion
comp_idxs = (edge_lim+1):(N-edge_lim-1); % comparison indicies

taps2 = 16;
edge_lim2 = taps2 + ceil(fs*(max(tau)-min(tau))/2);% clip signal to account for edge distortion
comp_idxs2 = (edge_lim2+1):(N-edge_lim2-1); % comparison indicies

%% construct the toeplitz matrix
delta = 1e-5;
freq_samples = exp_fac/(2*T*L)*linspace(-L,L,2*L+1);
F_basis = exp(1i*2*pi*t*freq_samples).*demod_phase;
A_eplitz = F_basis'*F_basis;
toeplitz_inv = inv(A_eplitz + (delta)*eye(2*L+1));

%% get canonical vectors of Toeplitz inverse via brute force
[X1,X2,X3,X4] = eval_canonical_vecs_brute_force(A_eplitz + delta*eye(2*L+1),2*L+1);
x0 = X1(1,1); % first entry of first column of the inverse


%% equispaced sampled basis matrix
t_nyq = [0:N-1]/fs; % nyquist samples for testing
psi = exp(1i*2*pi*t_nyq'*freq_samples);

%% chirp filter initialization
%------- first stage chirp-z
A_stage1 = exp(-1i*pi*exp_fac);
W_czt_stage1 = exp(-1i*pi*exp_fac/L);
m_stage1 = N; % input signal length
k1_1 = 2*L+1; % chirp-z length
n_len = sqrt(M);
nfft = 2^nextpow2(m_stage1+k1_1-1);
%------- Premultiply data.
kk = ((-m_stage1+1):max(k1_1-1,m_stage1-1)).';
kk2 = (kk .^ 2) ./ 2;
ww = W_czt_stage1 .^ (kk2);   % <----- Chirp filter is 1./ww
nn = (0:(m_stage1-1))';
aa = A_stage1 .^ ( -nn );
aa_stage1 = aa.*ww(m_stage1+nn);

fv_stage1 = fft( 1 ./ ww(1:(k1_1-1+m_stage1)), nfft );   % <----- Chirp filter.

%--------- second stage chirp-z 
m_stage2 = sqrt(M); % input signal length
k1_2 = B+1; % chirp-z length
n_len2 = sqrt(M);
nfft_2 = 2^nextpow2(m_stage2+k1_2-1);

A_stage2 = 1;
kk = ((-m_stage2+1):max(k1_2-1,m_stage2-1)).';
kk2 = (kk .^ 2) ./ 2;
nn = (0:(m_stage2-1))';

W_czt_stage2_coeffs = zeros(2*L+1,1);
aa_stage2 = zeros(m_stage2,2*L+1);
FV_stage2 = zeros(nfft_2,2*L+1);
ww_stage2 = zeros(length(kk2),2*L+1);
for ll=1:2*L+1
    l_dash = ll - 1;
    W_czt_stage2_coeffs(ll,1) = exp(1i*pi*(1/(B))*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L)));
    ww_stage2(:,ll) = W_czt_stage2_coeffs(ll,1) .^ (kk2);   % <----- Chirp filter is 1./ww
    aa = A_stage2 .^ ( -nn );
    aa_stage2(:,ll) = aa.*ww_stage2(m_stage2+nn,ll);
    FV_stage2(:,ll) = fft( 1 ./ ww_stage2(1:(k1_2-1+m_stage2),ll), nfft_2 );   % <----- Chirp filter.
end

%------------ third stage chirp-z 
m_stage3 = sqrt(M); % input signal length
k1_3 = 2*B+1;
n_len3 = B+1;
nfft_3 = 2^nextpow2(m_stage3+k1_3-1);
kk = ((-m_stage3+1):max(k1_3-1,m_stage3-1)).';
kk2 = (kk .^ 2) ./ 2;
nn = (0:(m_stage3-1))';

A_coeff_stage3 = zeros(2*L+1,1);
W_czt_stage3_coeffs = zeros(2*L+1,1);
aa_stage3 = zeros(m_stage3,2*L+1);
FV_stage3 = zeros(nfft_3,2*L+1);
ww_stage3 = zeros(length(kk2),2*L+1);
for ll=1:2*L+1
    l_dash = ll - 1;
    A_coeff_stage3(ll,1) = exp(1i*pi*(1/(B))*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L))*B);
    W_czt_stage3_coeffs(ll,1) = exp(1i*pi*(1/(B))*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L)));
    ww_stage3(:,ll) = W_czt_stage3_coeffs(ll,1) .^ (kk2);   % <----- Chirp filter is 1./ww
    aa = A_coeff_stage3(ll,1) .^ ( -nn );
    aa_stage3(:,ll) = aa.*ww_stage3(m_stage3+nn,ll);
    FV_stage3(:,ll) = fft( 1 ./ ww_stage3(1:(k1_3-1+m_stage3),ll), nfft_3 );   % <----- Chirp filter.
end


%% simulation specifications
trials = 4;
SNRs = -30:6:30;
SNR_fbst = zeros(trials,length(SNRs));
SNR_ttd16 = zeros(trials, length(SNRs));
SNR_ttd32 = zeros(trials, length(SNRs));
SNR_fds16 = zeros(trials, length(SNRs));
SNR_fds32 = zeros(trials, length(SNRs));

SNR_DFT_beamformer = zeros(trials, length(SNRs));

for ii=1:trials
    for jj=1:length(SNRs)

        SNR_dB = SNRs(jj); % signal to noise ratio
        a = (randn(n_sinusoids,1) + 1i*randn(n_sinusoids,1))/sqrt(2); % sinusoid amplitudes
        f = W*(2*rand(n_sinusoids,1)-1); % sinusoid frequencies
        s = sig_gen(a,f,n_sinusoids,t);
        sigma = (norm(s)/sqrt(MN))*10^(-SNR_dB/20);
        noise = sigma*(randn(MN,1)+1i*randn(MN,1))/sqrt(2);
        SNR_check = db(norm(s)/norm(noise));
        fprintf('Trial: %d, Set SNR: %.2f, Measured SNR: %.2f\n',ii,SNR_dB,SNR_check)

        s_demod = s.*demod_phase + noise;
        y = s_demod;
        y_2d = reshape(y,[M,N]);
        y_3d = reshape(y_2d, [sqrt(M),sqrt(M),N]);

        %% chirp-z FBST
        F_transform = zeros(sqrt(M),sqrt(M),2*L+1); % transform w.r.t. temporal frequency
        % A = -1;%exp(-1i*2*pi*T*W);
        % W_czt = exp(-1i*pi/L);
        for mm1=1:sqrt(M)
            y_in = reshape(y_3d(mm1,:,:),n_len,N).';
            y_in = y_in .* aa_stage1(:,ones(1,n_len));
            fy = fft(  y_in, nfft );
            fy = fy .* fv_stage1(:,ones(1, n_len));
            g  = ifft( fy );
            g = g( m_stage1:(m_stage1+k1_1-1), :, :) .* ww( m_stage1:(m_stage1+k1_1-1),ones(1, n_len) );
            F_transform(mm1,:,:) = g.';
        end

        X2_transform = zeros(sqrt(M),B+1,2*L+1); % transform w.r.t. second spatial frequency (J)
        A= 1;
        B_list1 = reshape(linspace(0,2*B,2*B+1), [1,2*B+1]);
        B_list2 = reshape(linspace(0,B,B+1), [1,B+1]);
        for ll=1:2*L+1
                l_dash = ll-1;

                y_in = reshape(F_transform(:,:,ll),sqrt(M),sqrt(M)).';
                y_in = y_in .* reshape(aa_stage2(:,ll,ones(1,n_len2)),sqrt(M),n_len2);
                fy = fft(  y_in, nfft_2 );
                fy = fy .* reshape(FV_stage2(:,ll,ones(1, n_len2)),nfft_2,n_len2);
                g  = ifft( fy );
                g = g( m_stage2:(m_stage2+k1_2-1), :, :) .* reshape(ww_stage2( m_stage2:(m_stage2+k1_2-1),ll,ones(1, n_len2) ),k1_2,n_len2);

                % A = exp(1i*pi*(1/(B))*(fc_norm + (1/(2*fc_tilde*T*L))*(l_dash-L))*B);%exp(-1i*pi*(1 + 2*(W/fc)*(l_dash-L/2)/L));
                % W_czt = exp(1i*pi*(1/(B))*(fc_norm + (1/(2*fc_tilde*T*L))*(l_dash-L)));
                coeff = exp(-1i*pi*(1/(B))*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L))*(sqrt(M)/2-1/2)*B_list2); % adjusting for array center
                % coeff2 = exp(1i*pi*(1/B)*(fc_norm + (1/(2*fc_tilde*T*L))*(l_dash-L))*(sqrt(M)/2-1/2)*B); % adjusting for negative azimuthal
                X2_transform(:,:,ll) = g.';
                X2_transform(:,:,ll) = X2_transform(:,:,ll).*coeff;
        end

        X1_transform = zeros(2*B+1,B+1,2*L+1); % final transform matrix (after taking czt along the first spatial frequency) (K)
        for ll=1:2*L+1
                l_dash = ll-1;
                
                y_in = reshape(X2_transform(:,:,ll),sqrt(M),B+1);
                y_in = y_in .* reshape(aa_stage3(:,ll,ones(1,n_len3)),sqrt(M),n_len3);
                fy = fft(  y_in, nfft_3 );
                fy = fy .* reshape(FV_stage3(:,ll,ones(1, n_len3)),nfft_3,n_len3);
                g  = ifft( fy );
                g = g( m_stage3:(m_stage3+k1_3-1), :, :) .* reshape(ww_stage3( m_stage3:(m_stage3+k1_3-1),ll,ones(1, n_len3) ),k1_3,n_len3);

                % A = exp(1i*pi*(1/(B))*(fc_norm + (1/(2*fc_tilde*T*L))*(l_dash-L))*B); % check
                % W_czt = exp(1i*pi*(1/(B))*(fc_norm + (1/(2*fc_tilde*T*L))*(l_dash-L)));
                coeff = exp(-1i*pi*(1/(B))*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L))*(sqrt(M)/2-1/2)*B_list1); % adjusting for array center
                coeff2 = exp(1i*pi*(1/B)*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L))*(sqrt(M)/2-1/2)*B); % adjusting for negative azimuthal
                X1_transform(:,:,ll) = g;
                X1_transform(:,:,ll) = coeff2*X1_transform(:,:,ll).*coeff.';
        end

        w_l = reshape(X1_transform(k_idx,j_idx,:),[2*L+1,1]);


        % %% sanity check
        % ls = reshape(linspace(-L,L-1,2*L),[2*L,1]);
        % w1 = -pi*(fc_norm + (1/(2*fc_tilde*T*L))*ls)*sin(theta)*sin(phi);
        % w2 = -pi*(fc_norm + (1/(2*fc_tilde*T*L))*ls)*sin(theta)*cos(phi);
        % w3 = (pi/L)*ls;
        % F_L = ls*(1/2)*(1/(L*T));
        % w_l_check = zeros(2*L,1);
        % for ll=1:2*L
        %     sum_l = 0;
        %     for mm1=1:sqrt(M)
        %         for mm2 =1:sqrt(M) 
        %             for nn=1:N
        %                 sum_l = sum_l + y_3d(mm1,mm2,nn)*exp(-1i*w1(ll)*(m(mm1)))*exp(-1i*w2(ll)*(m(mm2)))*exp(-1i*w3(ll)*(nn-1));
        %             end
        %         end
        %     end
        %     w_l_check(ll) = sum_l;
        % end

        %% superfast toeplitz using canonical vectors
        b1 = toep_mult(X2,w_l);
        b2 = toep_mult(X4,w_l);
        b3 = toep_mult(X1,b1);
        b4 = toep_mult(X3,b2);
        alpha_superfast = (b3 - b4)/x0;
        y_nyq_superfast = psi*alpha_superfast;

        alpha = toeplitz_inv*w_l;
        y_nyq = psi*alpha;

        y_filter_ttd = TTD_Beamformer((y_2d.*mod_phase).',-tau, fs, @Psi_sinc, 'R',taps).';
        y_bf_ttd = sum(y_filter_ttd(:,1:N),1)/M;

        y_filter_ttd32 = TTD_Beamformer((y_2d.*mod_phase).',-tau, fs, @Psi_sinc, 'R',taps2).';
        y_bf_ttd32 = sum(y_filter_ttd32(:,1:N),1)/M;

        s_nyq = sig_gen(a,f,n_sinusoids,t_nyq');
        s_nyq_2 = sig_gen(a,f,n_sinusoids,t_nyq+1/fs);
        s_nyq_3 = sig_gen(a,f,n_sinusoids,t_nyq+2*S/fs);
        
        %---------R=16 Fast delay and sum
        y_multibeamformed_first = zeros(sqrt(M),N,1,R^S);
        y_3d_fds = reshape(y_2d.*mod_phase,[sqrt(M),sqrt(M),N]);
        %% first FDS along b1 
        for m_a1 = 1:sqrt(M)
            v_now = reshape(y_3d_fds(m_a1,:,:,1,1),[sqrt(M),N,1]);
            for stage=1:S
                v_prev = v_now;
                num_virtual_sensors = sqrt(M)/2^stage;
                num_sectors = R^stage;
                v_now = zeros(num_virtual_sensors,N,num_sectors); % init partial beamforming at stage s

                for mm=0:num_virtual_sensors-1
                    for rr = 0:num_sectors-1
                        m = mm+1;
                        r = floor((rr)/R)+1;
                        delays = ((-1)^(rr))*(1/(R^stage))*eta*[2^(stage-1)*(2*mm+0.5)-sqrt(M)/2; 2^(stage-1)*(2*mm+1.5)-sqrt(M)/2];
                        y_interp = reshape(v_prev([2*mm+1,2*mm+2],:,r),[2,N]);
                        v_prev_interped = TTD_Beamformer(y_interp.',-delays, fs, @Psi_sinc, 'R',taps).';
                        v_now(m,:,rr+1) = sum(v_prev_interped,1)./2;
                    end
                end
            end
            y_multibeamformed_first(m_a1,:,1,:) = v_now;
        end

        y_multibeamformed_total = zeros(N,R^S,R^S);
        %% second FDS along a1
        for r_b1=1:R^S
            v_now = reshape(y_multibeamformed_first(:,:,1,r_b1),[sqrt(M),N,1]);
            for stage=1:S
                v_prev = v_now;
                num_virtual_sensors = sqrt(M)/2^stage;
                num_sectors = R^stage;
                v_now = zeros(num_virtual_sensors,N,num_sectors); % init partial beamforming at stage s

                for mm=0:num_virtual_sensors-1
                    for rr = 0:num_sectors-1
                        m = mm+1;
                        r = floor((rr)/R)+1;
                        delays = ((-1)^(rr))*(1/(R^stage))*eta*[2^(stage-1)*(2*mm+0.5)-sqrt(M)/2; 2^(stage-1)*(2*mm+1.5)-sqrt(M)/2];
                        y_interp = reshape(v_prev([2*mm+1,2*mm+2],:,r),[2,N]);
                        v_prev_interped = TTD_Beamformer(y_interp.',-delays, fs, @Psi_sinc, 'R',taps).';
                        v_now(m,:,rr+1) = sum(v_prev_interped,1)./2;
                    end
                end
            end
            y_multibeamformed_total(:,:,r_b1) = v_now;
        end
        %--- adjust for positive a1 constraint
        if a_idx<=R^S/2 
            v_pred = y_multibeamformed_total(:,b_idx,a_idx);
        else
            v_pred = y_multibeamformed_total(:,R^S - b_idx+1,R^S - a_idx+1);
        end
        SNR_fds16(ii,jj) = norm(s_nyq_3(comp_idxs))/norm(s_nyq_3(comp_idxs) - v_pred(comp_idxs).');

        %---------R=32 Fast delay and sum
        y_multibeamformed_first = zeros(sqrt(M),N,1,R^S);
        y_3d_fds = reshape(y_2d.*mod_phase,[sqrt(M),sqrt(M),N]);
        %% first FDS along b1 
        for m_a1 = 1:sqrt(M)
            v_now = reshape(y_3d_fds(m_a1,:,:,1,1),[sqrt(M),N,1]);
            for stage=1:S
                v_prev = v_now;
                num_virtual_sensors = sqrt(M)/2^stage;
                num_sectors = R^stage;
                v_now = zeros(num_virtual_sensors,N,num_sectors); % init partial beamforming at stage s

                for mm=0:num_virtual_sensors-1
                    for rr = 0:num_sectors-1
                        m = mm+1;
                        r = floor((rr)/R)+1;
                        delays = ((-1)^(rr))*(1/(R^stage))*eta*[2^(stage-1)*(2*mm+0.5)-sqrt(M)/2; 2^(stage-1)*(2*mm+1.5)-sqrt(M)/2];
                        y_interp = reshape(v_prev([2*mm+1,2*mm+2],:,r),[2,N]);
                        v_prev_interped = TTD_Beamformer(y_interp.',-delays, fs, @Psi_sinc, 'R',taps2).';
                        v_now(m,:,rr+1) = sum(v_prev_interped,1)./2;
                    end
                end
            end
            y_multibeamformed_first(m_a1,:,1,:) = v_now;
        end

        y_multibeamformed_total = zeros(N,R^S,R^S);
        %% second FDS along a1
        for r_b1=1:R^S
            v_now = reshape(y_multibeamformed_first(:,:,1,r_b1),[sqrt(M),N,1]);
            for stage=1:S
                v_prev = v_now;
                num_virtual_sensors = sqrt(M)/2^stage;
                num_sectors = R^stage;
                v_now = zeros(num_virtual_sensors,N,num_sectors); % init partial beamforming at stage s

                for mm=0:num_virtual_sensors-1
                    for rr = 0:num_sectors-1
                        m = mm+1;
                        r = floor((rr)/R)+1;
                        delays = ((-1)^(rr))*(1/(R^stage))*eta*[2^(stage-1)*(2*mm+0.5)-sqrt(M)/2; 2^(stage-1)*(2*mm+1.5)-sqrt(M)/2];
                        y_interp = reshape(v_prev([2*mm+1,2*mm+2],:,r),[2,N]);
                        v_prev_interped = TTD_Beamformer(y_interp.',-delays, fs, @Psi_sinc, 'R',taps2).';
                        v_now(m,:,rr+1) = sum(v_prev_interped,1)./2;
                    end
                end
            end
            y_multibeamformed_total(:,:,r_b1) = v_now;
        end
        %--- adjust for positive a1 constraint
        if a_idx<=R^S/2 
            v_pred = y_multibeamformed_total(:,b_idx,a_idx);
        else
            v_pred = y_multibeamformed_total(:,R^S - b_idx+1,R^S - a_idx+1);
        end
        SNR_fds32(ii,jj) = norm(s_nyq_3(comp_idxs2))/norm(s_nyq_3(comp_idxs2) - v_pred(comp_idxs2).');



        %% sub-band processing
        % subband FFT across array elements 
        Y_f = fft(y_2d.').'; % mth column for mth array sample
        Y_f_3d = reshape(Y_f, [sqrt(M),sqrt(M),N]);

        X_transform_subband = zeros(2*B+1,B+1,N);
        beam_list2 = linspace(0,B,B+1);
        beam_list1 = linspace(0,2*B,2*B+1);
        Y_nb = zeros(N,1);
        f_set = zeros(N,1);

        % beams for each frequency bin sequentially
        for nn=0:N-1
            y_array_space = Y_f_3d(:,:,nn+1).';
            if nn<N/2
                f_curr = fc + nn*freq_resolution;
            elseif nn>=N/2
                f_curr = fc + (nn)*freq_resolution - N*freq_resolution;
            end
            f_set(nn+1) = f_curr;
            % take czt for each column of y_array space first
            A_czt2 = 1;%exp(1i*pi*f_curr/(fc+W));
            W_czt2 = exp(1i*pi*f_curr/((fc+W)*B));
            coeff2_1 = exp(-1i*pi*f_curr*((sqrt(M)/2-1/2)/(fc+W))*(beam_list2/B)).';
            temp_var = czt(y_array_space,B+1,W_czt2,A_czt2);
            temp_var = (temp_var.*coeff2_1).';

            % take czt for  each column of temp_var second
            A_czt1 = exp(1i*pi*f_curr/(fc+W));
            W_czt1 = exp(1i*pi*f_curr/((fc+W)*B));
            coeff1_1 = exp(-1i*pi*f_curr*((sqrt(M)/2-1/2)/(fc+W))*(beam_list1/B)).';
            coeff1_2 = exp(1i*pi*f_curr*(sqrt(M)/2-1/2)/(fc+W));
            temp_var2 = czt(temp_var,2*B+1,W_czt1,A_czt1);
            X_transform_subband(:,:,nn+1) = coeff1_2*(temp_var2.*coeff1_1);
        end

        % take ifft back across all beams
        Y_time_beamspace = ifft(X_transform_subband,[],3)/M;
        y_nyq_subband = squeeze(Y_time_beamspace(k_idx,j_idx,:));


        SNR_fbst(ii,jj) = norm(s_nyq)/norm(s_nyq - y_nyq_superfast);
        SNR_DFT_beamformer(ii,jj) = norm(s_nyq)/norm(s_nyq - y_nyq_subband);
        SNR_ttd16(ii,jj) = norm(s_nyq_2(comp_idxs))/norm(s_nyq_2(comp_idxs) - y_bf_ttd(comp_idxs));
        SNR_ttd32(ii,jj) = norm(s_nyq_2(comp_idxs2))/norm(s_nyq_2(comp_idxs2) - y_bf_ttd32(comp_idxs2));
    end
end

load('SNR_subband_UPA_N256.mat');

figure(1)
plot(SNRs,SNRs + db(M)/2,'LineWidth',1)
hold on
grid on
plot(SNRs,db(mean(SNR_fbst)),'-*','LineWidth',1)
hold on
grid on
plot(SNRs, db(mean(SNR_ttd16)),'-^','LineWidth',1)
hold on
grid on
plot(SNRs, db(mean(SNR_ttd32)),'-^','LineWidth',1)
hold on
grid on
plot(SNRs, db(mean(SNR_fds16)),'-^','LineWidth',1)
hold on
grid on
plot(SNRs, db(mean(SNR_fds32)),'-^','LineWidth',1)
hold on
grid on
plot(SNRs,db(mean(SNR_DFT_beamformer)),'-^','LineWidth',1);
hold on
grid on
plot(SNRs,db(mean(SNR_DFT_beamformer2)),'-^','LineWidth',1);
xlabel('Nominal SNR (dB)','Interpreter','latex','FontSize',12)
ylabel('Beamformed SNR (dB)','Interpreter','latex','FontSize',12)
legend({'Ideal', 'FBST', 'Delay and Sum (R=16)', 'Delay and Sum (R=32)', 'Fast Delay and Sum (R=16)', 'Fast Delay and Sum (R=32)','Sub-Band Processing (N=64)', 'Sub-Band Processing (N=256)'},'Location','northwest','Interpreter','latex','FontSize',10);
% exportgraphics(gcf, 'upa_ag_vs_snr_subband.pdf', 'ContentType', 'vector')
%% supporting functions
function [X] = sig_gen(a,f,n_sinusoids,t)
X = zeros(size(t));
for ii = 1:n_sinusoids % For each antenna
    sig = (a(ii)*exp(1i*2*pi*(f(ii))*(t)));
    X = X + sig;
end
end

function y = Psi_sinc(z)
y = sinc(z);
end

