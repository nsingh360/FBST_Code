clear
clc
close all

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

B = 2^6; % no. of beams
M = 2^7; % array size
N = 2^6; % no. of samples
L = ceil(N); % no. of samples in frequency is 2*L
exp_fac = 1.005; % factor to include freqs slightly out of band;
MN = M*N;

Theta = zeros(2*B+1,1);
idx = 121;%2*B;%randi(B,1);
for i=1:2*B+1
    Theta(i) = asin((i-1-B)/B);
end

%% parameters for fast delay and sum
S = log2(M);
R_fac = 2; % 2-radix fast delay and sum
L_fac = 2; % downsampling factor in spatial dimension

%% Angle sampling for final stage
r_lin = linspace(0,R_fac^S-1,R_fac^S);
Theta_fds = asin(1 - (2*r_lin + 1)/R_fac^S);
idx_fds = 8;

fc = 20e9;  % center frequency
Ws = 5e9; % to prevent spacial aliasing
c = physconst('LightSpeed');
lambda = c/(fc+Ws); % wavelngth
eta = lambda/(2*c);
W = 5e9; % bandwidth
fs = 2*W; % sampling frequency
T = 1/fs;
theta = Theta(idx);
m = [(-M/2+1/2):(M/2-1/2)]';
x_pos = m*lambda/2; % sensor positions
u_s = sin(theta); % normal vector for determining delays

freq_resolution = fs/N; % frequency grid for sub-band beamformign

%% Signal specs that do not need to be redifined in loop
n_sinusoids = 100; % number of sinusoids in signal
n = [0:N-1]/fs; % temporal sample vectors
tau = x_pos*u_s/c; % relative delays to phase center
t_array = n-tau; % delays across array and time
t = t_array(:); %
T_aperture = max(t)-min(t);% temporal aperture of the array
demod_phase = repmat(exp(-1i*2*pi*(fc)*tau),N,1);
mod_phase = exp(1i*2*pi*(fc)*tau);
% demod_phase = repmat(demod_phase,N,1);
taps = 8;
edge_lim = taps + ceil(fs*(max(tau)-min(tau))/2);% clip signal to account for edge distortion
comp_idxs = (edge_lim+1):(N-edge_lim-1); % comparison indicies

taps2 = 16;
edge_lim2 = taps2 + ceil(fs*(max(tau)-min(tau))/2);% clip signal to account for edge distortion
comp_idxs2 = (edge_lim2+1):(N-edge_lim2-1); % comparison indicies

edge_lim_dft = ceil(fs*(max(tau)-min(tau))/2);
comp_dft_idx = edge_lim_dft+1:(N-edge_lim_dft-1);

%% construct the toeplitz matrix
delta = 1e-6;
freq_samples = exp_fac/(2*T*L)*linspace(-L,L,2*L+1);
freq_samples2 = 1/(2*T*L)*linspace(-L,L,2*L+1);
F_basis = exp(1i*2*pi*t*freq_samples).*demod_phase;
F_basis2 = exp(1i*2*pi*t*freq_samples2).*demod_phase;
% for ll=1:2*L
%     for kk=1:2*L
%         f_l = (1/(2*T*L))*(ll-1-L);
%         f_k = (1/(2*T*L))*(kk-1-L);
%         A_eplitz(ll,kk) = sum(exp(1i*2*pi*(f_k - f_l)*t));
%     end
% end
A_eplitz = F_basis'*F_basis;
A_eplitz2 = F_basis2'*F_basis2;
toeplitz_inv = inv(A_eplitz + (delta)*eye(2*L+1));

%% get canonical vectors of Toeplitz inverse via brute force
[X1,X2,X3,X4] = eval_canonical_vecs_brute_force(A_eplitz + delta*eye(2*L+1),2*L+1);
x0 = X1(1,1); % first entry of first column of the inverse

% %% construct slepian matrix
% [A_basis, V_nyq, K] = construct_fast_slepian_basis(W, t, T, L, M, N, T_aperture);
% Phi = ((A_basis'*A_basis + delta*eye(K))^-1)*(A_basis');

%% slepian function for signal generation
K = ceil(2*W*T_aperture)+10;% subspace dimension
[V_u,Lambda] = dpss(MN+1,W*T_aperture,K); % firt K uniformly sampled PSWF's
t_u = [0:(MN)]'; % uniform sampling points over the interval
t_nu = [M*N*(t-min(t(:)))/T_aperture]; % non-uniform sample points, scaled to match interval
V_nu = interp1(t_u,V_u,t_nu,'pchip',0); % non-uniform
% Lambda = svd(V_nu);

%% equispaced sampled basis matrix
t_nyq = [0:N-1]/fs; % nyquist samples for testing
t_nyq2 = [0:N-1]/fs+1/fs;
t_nyq3 =  [0:N-1]/fs + S/fs;
V_nyq = interp1(t_u,V_u,MN*(t_nyq-min(t))/T_aperture,'pchip',0); % interpolate to nyquist
V_nyq2 = interp1(t_u,V_u,MN*(t_nyq2-min(t))/T_aperture,'pchip',0); % interpolate to nyquist
V_nyq3 = interp1(t_u,V_u,MN*(t_nyq3-min(t))/T_aperture,'pchip',0); % interpolate to nyquist
psi = exp(1i*2*pi*t_nyq'*freq_samples);
psi2 = exp(1i*2*pi*t_nyq'*freq_samples2);


%% slepian filter
K_filt = 2^4; % filter taps
N_filter = K_filt + N + 1; % filter length
dim = ceil(2*W*T_aperture)+8;% subspace dimension
K_mid = K_filt/2;

t_filt = [0:K_filt-1]/fs - tau;
t_filt = t_filt(:);
K_interp = 2*MN;
K_interp2 = 2*MN;
V_u_filt = dpss(K_interp+1,W*T_aperture,dim); % firt K uniformly sampled PSWF's
V_u_filt2 = dpss(K_interp2+1,W*T_aperture,dim); % firt K uniformly sampled PSWF's
t_u_filt = [0:(K_interp)]'; % uniform sampling points over the interval
t_nu_filt = [K_interp*(t_filt-min(t_filt(:)))/T_aperture]; % non-uniform sample points, scaled to match interval
V_nu_filt = interp1(t_u_filt,V_u_filt,t_nu_filt,'pchip',0); % non-uniform

t_filt_nyq = [0:K_filt-1]'/fs; % nyquist samples for testing
V_nyq_filt = interp1(t_u_filt,V_u_filt,K_interp*(t_filt_nyq-min(t_filt))/T_aperture,'pchip',0); % interpolate to nyquist
v_nyq = interp1(t_u_filt,V_u_filt,K_interp*(t_filt_nyq(K_mid)-min(t_filt))/T_aperture,'pchip',0);

%% Filter design
delta = 1e-5;
Phi_filt = ((V_nu_filt'*V_nu_filt + delta*eye(dim))^-1)*(V_nu_filt');
h = v_nyq*Phi_filt;
H = reshape(h,M,K_filt);
H_buff = circshift([H,zeros(M,N+1)],-K_mid+1,2);
Hf = fft(H_buff,[],2);
t_nyq = [0:N-1]/fs;
comp_idxs = (K_filt+1):(N-K_filt-1);


%% simulation specifications
trials = 25;
SNRs = -30:6:30;
SNR_fbst = zeros(trials,length(SNRs));
SNR_fbst_2 = zeros(trials, length(SNRs)); % fbst with no basis over the bandwidth
SNR_ttd32 = zeros(trials,length(SNRs));
SNR_ttd = zeros(trials, length(SNRs));
SNR_fbst_superfast_toeplitz = zeros(trials, length(SNRs));
SNR_fds = zeros(trials, length(SNRs));
SNR_fds32 = zeros(trials, length(SNRs));
SNR_DFT = zeros(trials, length(SNRs));

SNR_slepian_filt32 = zeros(trials, length(SNRs));

for ii = 1:trials
    for jj = 1:length(SNRs)

        SNR_dB = SNRs(jj); % signal to noise ratio
        a = (randn(n_sinusoids,K) + 1i*randn(n_sinusoids,K))/sqrt(2).*(Lambda.'); % sinusoid amplitudes
        % f = linspace(-W,W,n_sinusoids); % sinusoid frequencies
        s = slepian_sig_gen(V_nu,a,n_sinusoids,t);
        sigma = (norm(s)/sqrt(MN))*10^(-SNR_dB/20);
        noise = sigma*(randn(MN,1)+1i*randn(MN,1))/sqrt(2);
        SNR_check = db(norm(s)/norm(noise));
        fprintf('Trial: %d, Set SNR: %.2f, Measured SNR: %.2f\n',ii,SNR_dB,SNR_check)

        s_demod = s.*demod_phase + noise;

        y = s_demod;
        y = reshape(y, [M,N]);
        % F_transform = zeros(M,2*L);
        % A = -1;%exp(-1i*2*pi*T*W);
        % W_czt = exp(-1i*pi/L);
        % tic;
        % for mm=1:M
        %     F_transform(mm,:) = czt(y(mm,:),2*L,W_czt,A);
        % end
        % 
        % X_czt = zeros(2*L,B);
        % A= 1;
        % for ll=1:2*L
        %     l_dash = ll-1;
        %     % A = exp(-1i*pi*(1 + 2*(W/fc)*(l_dash-L/2)/L));
        %     W_czt = exp(1i*pi*(1/B)*(1 + (1/(2*fc*T*L))*(l_dash-L)));
        %     X_czt(ll,:) = czt(F_transform(:,ll),B,W_czt,A);
        % end
        % 
        % w_l = X_czt(:,idx);

        %% sanity check
        % y_check = zeros(M,N);
        % ls = reshape(linspace(-L,L-1,2*L),[2*L,1]);
        % w1 = -pi*(1 + (1/(2*fc*T*L))*ls)*sin(Theta(idx));
        % w2 = (pi/L)*ls;
        % F_L = ls*(1/2)*(1/(L*T));
        % for mm=1:M
        %     for nn=1:N
        %         y_check(mm,nn) = sum(w_l.*exp(1i*w1*mm).*exp(1i*w2*nn));
        %     end
        % end
        % w_l_check = zeros(2*L,1);
        % for ll=1:2*L
        %     sum_l = 0;
        %     for mm=1:M
        %         for nn=1:N
        %             sum_l = sum_l + y(mm,nn)*exp(-1i*w1(ll)*(m(mm)))*exp(-1i*w2(ll)*(nn-1));
        %         end
        %     end
        %     w_l_check(ll) = sum_l;
        % end
        % 
        % inner_prod_basis = zeros(2*L,MN);
        % for ll=1:2*L
        %     for mn=1:MN
        %         inner_prod_basis(ll,mn) = exp(-1i*2*pi*F_L(ll)*t(mn));
        %     end
        % end
        % w_l_check2 = inner_prod_basis*s;

        [y_nyq, w_l, alpha] = fbst_1(y, L, B, A_eplitz, psi, idx, fc,Ws, T, delta, exp_fac);
        % [y_nyq2, w_l2, alpha2] = fbst_1(y, L, B, A_eplitz2, psi2, idx, fc,Ws, T, delta, 1);

        % alpha_fbst = inv(A_eplitz + delta*eye(2*L))*F_basis'*y(:);
        % y_nyq = psi*alpha_fbst;

        % %% superfast toeplitz using canonical vectors
        % c_x1 = X1(:, 1);
        % r_x1 = X1(1, :);
        % p_x1 = [c_x1; 0; flipud(conj(r_x1(2:end)'))];
        % pf_x1 = fft(p_x1);
        % 
        % c_x2 = X2(:, 1);
        % r_x2 = X2(1, :);
        % p_x2 = [c_x2; 0; flipud(conj(r_x2(2:end)'))];
        % pf_x2 = fft(p_x2);
        % 
        % c_x3 = X3(:, 1);
        % r_x3 = X3(1, :);
        % p_x3 = [c_x3; 0; flipud(conj(r_x3(2:end)'))];
        % pf_x3 = fft(p_x3);
        % 
        % c_x4 = X4(:, 1);
        % r_x4 = X4(1, :);
        % p_x4 = [c_x4; 0; flipud(conj(r_x4(2:end)'))];
        % pf_x4 = fft(p_x4);
        % 
        % fft_wl = fft([w_l; zeros(length(w_l), 1)]);
        % b1 = ifft(pf_x2 .* fft_wl);
        % b1 = b1(1:length(w_l));
        % % b1 = toep_mult(X2,w_l);
        % b2 = ifft(pf_x4 .* fft_wl);
        % b2 = b2(1:length(w_l));
        % % b2 = toep_mult(X4,w_l);
        % b3 = ifft(pf_x1 .* fft([b1; zeros(length(b1), 1)]));
        % b3 = b3(1:length(b1));
        % % b3 = toep_mult(X1,b1);
        % b4 = ifft(pf_x3 .* fft([b2; zeros(length(b2), 1)]));
        % b4 = b4(1:length(b2));
        % alpha_superfast = (b3 - b4)/x0;
        % y_nyq_superfast = psi*alpha_superfast;
        % coeff_nyq = exp(-1i*pi*linspace(0,N-1,N)).';
        % y_nyq_transform = czt(alpha_superfast,N,exp(1i*pi/L),1);
        % y_nyq_transform = y_nyq_transform.*coeff_nyq;

        
        % tic;
        % for bb=1:B
        %     [y_nyq, w_l] = fbst_0(y, L, B, toeplitz_inv, psi, idx, fc, T);
        % end
        % runtime_fbst0 = toc;

        s_nyq = slepian_sig_gen(V_nyq,a,n_sinusoids,t_nyq');
        s_nyq_2 = slepian_sig_gen(V_nyq2,a,n_sinusoids,t_nyq+1/fs);
        s_nyq_3 = slepian_sig_gen(V_nyq3,a,n_sinusoids,t_nyq'+S/fs);

        y_filter_ttd = TTD_Beamformer((y.*mod_phase).',-tau, fs, @Psi_sinc, 'R',taps).';
        y_bf_ttd = (sum(y_filter_ttd(:,1:N),1)/M).';
        % 
        y_filter_ttd32 = TTD_Beamformer((y.*mod_phase).',-tau, fs, @Psi_sinc, 'R',taps2).';
        y_bf_ttd32 = (sum(y_filter_ttd32(:,1:N),1)/M).';
        % 
        %% fast delay and sum beamforming 16 taps
        v_now = y.*mod_phase;
        for stage=1:S
            v_prev = v_now;
            num_virtual_sensors = M/2^stage;
            num_sectors = R_fac^stage;
            v_now = zeros(num_virtual_sensors,N,num_sectors); % init partial beamforming at stage s
            for mm2=0:num_virtual_sensors-1
                for rr = 0:num_sectors-1
                    m2 = mm2+1;
                    r = floor((rr)/R_fac)+1;
                    delays = ((-1)^(rr))*(1/(R_fac^stage))*eta*[2^(stage-1)*(2*mm2+0.5)-M/2; 2^(stage-1)*(2*mm2+1.5)-M/2];
                    y_interp = reshape(v_prev([2*mm2+1,2*mm2+2],:,r),[2,N]);
                    v_prev_interped = TTD_Beamformer(y_interp.',-delays, fs, @Psi_sinc, 'R',taps).';
                    v_now(m2,:,rr+1) = sum(v_prev_interped,1)./2;
                end
            end
            % for mm2=0:num_virtual_sensors-1
            %     for rr = 0:num_sectors-1
            %         m2 = mm2+1;
            %         r = floor((rr)/R_fac)+1;
            %         delays = ((-1)^(rr))*(1/(R_fac^stage))*eta*[2^(stage-1)*(2*mm2+0.5)-M/2; 2^(stage-1)*(2*mm2+1.5)-M/2];
            %         v_prev_interped = TTD_Beamformer(y_interp.',-delays, fs, @Psi_sinc, 'R',taps2).';
            %         % y_interp = reshape(v_prev([2*mm2+1,2*mm2+2],:,r),[2,N]);
            %         % t_filt2 = [0:K_filt-1]/fs - delays;
            %         % t_filt2 = t_filt2(:);
            %         % t_u_filt = [0:(K_interp2)]'; % uniform sampling points over the interval
            %         % t_nu_filt2 = [K_interp2*(t_filt2-min(t_filt2(:)))/T_aperture]; % non-uniform sample points, scaled to match interval
            %         % V_nu_filt2 = interp1(t_u_filt,V_u_filt2,t_nu_filt2,'pchip',0); % non-uniform
            %         % 
            %         % V_nyq_filt2 = interp1(t_u_filt,V_u_filt2,K_interp2*(t_filt_nyq-min(t_filt2))/T_aperture,'pchip',0); % interpolate to nyquist
            %         % v_nyq2 = interp1(t_u_filt,V_u_filt2,K_interp2*(t_filt_nyq(K_mid)-min(t_filt2))/T_aperture,'pchip',0);
            %         % Phi_filt2 = ((V_nu_filt2'*V_nu_filt2 + delta*eye(dim))^-1)*(V_nu_filt2');
            %         % h2 = v_nyq2*Phi_filt2;
            %         % H2 = reshape(h2,2,K_filt);
            %         % H_buff2 = circshift([H2,zeros(2,N+1)],-K_mid+1,2);
            %         % Hf2 = fft(H_buff2,[],2);
            %         % 
            %         % Yf = fft(y_interp,N_filter,2);
            %         % 
            %         % v_prev_interped = ifft(Yf.*conj(Hf2),[],2);
            % 
            %         v_now(m2,:,rr+1) = sum(v_prev_interped(:,1:N),1).';
            %     end
            % end
        end
        v_pred = v_now(1,:,idx_fds).';
        SNR_fds(ii,jj) = norm(s_nyq_3(comp_idxs))/norm(v_pred(comp_idxs) - s_nyq_3(comp_idxs));
        % 
        % 
        %% fast delay and sum beamforming 32 taps
        v_now = reshape(y.*mod_phase, [M,N,1]);
        for stage=1:S
            v_prev = v_now;
            num_virtual_sensors = M/2^stage;
            num_sectors = R_fac^stage;
            v_now = zeros(num_virtual_sensors,N,num_sectors); % init partial beamforming at stage s

            for mm2=0:num_virtual_sensors-1
                for rr = 0:num_sectors-1
                    m2 = mm2+1;
                    r = floor((rr)/R_fac)+1;
                    delays = ((-1)^(rr))*(1/(R_fac^stage))*eta*[2^(stage-1)*(2*mm2+0.5)-M/2; 2^(stage-1)*(2*mm2+1.5)-M/2];
                    y_interp = reshape(v_prev([2*mm2+1,2*mm2+2],:,r),[2,N]);
                    v_prev_interped = TTD_Beamformer(y_interp.',-delays, fs, @Psi_sinc, 'R',taps2).';
                    v_now(m2,:,rr+1) = sum(v_prev_interped,1)./2;
                end
            end
        end
        v_pred = v_now(1,:,idx_fds).';
        SNR_fds32(ii,jj) = norm(s_nyq_3(comp_idxs2))/norm(v_pred(comp_idxs2) - s_nyq_3(comp_idxs2));
        
        % Yf = fft(y.*mod_phase,N_filter,2);
        % y_bf = sum(ifft(Yf.*conj(Hf),[],2),1);
        % y_comp = y_bf(comp_idxs).';

        %% sub-band beamforming
        % subband FFT across array elements 
        Y_f = fft(y.'); % mth column for mth array sample

        Y_nb = zeros(N,1);
        f_set = zeros(N,1);


        Y_freq_beamspace = zeros(N,2*B+1);
        beam_list = linspace(0,2*B,2*B+1);

        for nn=0:N-1
            ym = Y_f(nn+1,:).';
            if nn<N/2
                f_curr = fc + nn*freq_resolution;
            elseif nn>=N/2
                f_curr = fc + (nn)*freq_resolution - N*freq_resolution;
            end
            W_czt = exp(1i*pi*f_curr/((fc+W)*B));
            A_czt = exp(1i*pi*f_curr/(fc+W));
            coeff2 = exp(1i*pi*f_curr*(M/2-1/2)/(fc+W));
            coeff1 = exp(-1i*pi*f_curr*((M/2-1/2)/(fc+W))*(beam_list/B)).';
            temp_var = czt(ym,2*B+1,W_czt,A_czt);
            Y_freq_beamspace(nn+1,:) = coeff2*temp_var.*coeff1;
        end

        Y_time_beamspace = ifft(Y_freq_beamspace);
        Y_bf = Y_time_beamspace(:,idx)/M; % 


        SNR_ttd(ii,jj) = norm(s_nyq_2(comp_idxs))/norm(y_bf_ttd(comp_idxs) - s_nyq_2(comp_idxs));
        SNR_ttd32(ii,jj) = norm(s_nyq_2(comp_idxs2))/norm(y_bf_ttd32(comp_idxs2) - s_nyq_2(comp_idxs2));
        SNR_fbst(ii,jj) = norm(s_nyq)/norm(s_nyq - y_nyq);
        SNR_DFT(ii,jj) = norm(s_nyq)/norm(Y_bf - s_nyq);
        % SNR_fbst_2(ii,jj) = norm(s_nyq)/norm(s_nyq - y_nyq2);
        % SNR_fbst_superfast_toeplitz(ii,jj) = norm(s_nyq)/norm(s_nyq - y_nyq_superfast);
        % SNR_slepian_filt32(ii,jj) = norm(s_nyq(comp_idxs))/norm(s_nyq(comp_idxs) - y_comp);

    end
end

load('SNR_subband_N2048.mat');

figure(1)
plot(SNRs,SNRs + db(M)/2,'LineWidth',1)
hold on
grid on
plot(SNRs,db(mean(SNR_fbst)),'-*','LineWidth',1)
hold on
grid on
plot(SNRs,db(mean(SNR_ttd)),'-^','LineWidth',1);
hold on
grid on
plot(SNRs,db(mean(SNR_ttd32)),'-^','LineWidth',1);
hold on
grid on
plot(SNRs,db(mean(SNR_fds)),'-^','LineWidth',1);
hold on
grid on
plot(SNRs,db(mean(SNR_fds32)),'-^','LineWidth',1);
hold on
grid on
plot(SNRs,db(mean(SNR_DFT)),'-^','LineWidth', 1);
hold on
grid on
plot(SNRs, db(mean(SNR_subband)),'-^','LineWidth',1)
xlabel('Nominal SNR (dB)','Interpreter','latex','FontSize',12)
ylabel('Beamformed SNR (dB)','Interpreter','latex','FontSize',12)
legend({'Ideal','FBST', 'Delay and Sum (R=16)','Delay and Sum (R=32)','Fast Delay and Sum (R=16)',' Fast Delay and Sum (R=32)', 'Sub-Band Processing (N=64)', 'Sub-Band Processing (N=2048)'},'Location','northwest','Interpreter','latex','FontSize',11)
exportgraphics(gcf, 'ula_snr_with_sub_band.pdf', 'ContentType', 'vector');


%% supporting functions
function [X] = sig_gen(a,f,n_sinusoids,t)
X = zeros(size(t));
for ii = 1:n_sinusoids % For each antenna
    sig = (a(ii)*exp(1i*2*pi*(f(ii))*(t)));
    X = X + sig;
end
end

function y = Psi_b3(z)
y = (z.^3/6 + z.^2 + 2*z + 4/3).*(-2<=z).*(z<-1);
y = y + (-z.^3/2-z.^2 + 2/3).*(-1<=z).*(z<0);
y = y + (z.^3/2 - z.^2 + 2/3).*(0<=z).*(z<1);
y = y + (-z.^3/6 + z.^2 - 2*z + 4/3).*(1<=z).*(z<2);
end

function y = Psi_sinc(z)
y = sinc(z);
end

function [X] = slepian_sig_gen(V,a,n_slepians,t)
X = sum(a*V.').';
end