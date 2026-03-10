clear
clc
close all

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

B = 64; % no. of beams
M = 64; % array size
exp_fac=1.05;

Theta = zeros(2*B+1,1);
idx = 2*B;%2*B+1;%randi(B,1);
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
idx_fds = 9;

fc = 20e9;  % center frequency
Ws = 5e9; % to prevent spacial aliasing
c = physconst('LightSpeed');
lambda = c/(fc+Ws); % wavelngth
eta = lambda/(2*c);
W = 5e9; % bandwidth
N = 64;%ceil(1/(2/exp_fac-1)*((M-1)*W/(fc+W)-1)); % no. of samples
L = ceil(N); % no. of samples in frequency is 2*L
MN = M*N;
fs = 2*W; % sampling frequency
T = 1/fs;
theta = Theta(idx);
m = [(-M/2+1/2):(M/2-1/2)]';
x_pos = m*lambda/2; % sensor positions
u_s = sin(theta); % normal vector for determining delays

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

%% construct the toeplitz matrix
delta = 1e-5;
freq_samples = exp_fac/(2*T*L)*linspace(-L,L,2*L+1);
F_basis = exp(1i*2*pi*t*freq_samples).*demod_phase;
% for ll=1:2*L
%     for kk=1:2*L
%         f_l = (1/(2*T*L))*(ll-1-L);
%         f_k = (1/(2*T*L))*(kk-1-L);
%         A_eplitz(ll,kk) = sum(exp(1i*2*pi*(f_k - f_l)*t));
%     end
% end
A_eplitz = F_basis'*F_basis;
toeplitz_inv = inv(A_eplitz + (delta)*eye(2*L+1));
% Q = inv(A_eplitz + 1e-8*eye(2*L))*A_eplitz*inv(A_eplitz + 1e-8*eye(2*L));
% %% EVALUATE VARIANCE TERM
% sum_var = 0;
% for ll=1:2*L
%     for kk=1:2*L
%         f_l = (1/(2*T*L))*(ll-1-L);
%         f_k = (1/(2*T*L))*(kk-1-L);
%         if kk ~= ll
%             sum_var = sum_var + Q(ll,kk)*exp(1i*pi*(f_k - f_l)*T_aperture)*sin(pi*(f_k - f_l)*T_aperture)/(2*pi*(f_k-f_l));
%         else
%             sum_var = sum_var + Q(ll,kk)*exp(1i*pi*(f_k - f_l)*T_aperture)*T_aperture/2;
%         end
%     end
% end

%% get canonical vectors of Toeplitz inverse via brute force
[X1,X2,X3,X4] = eval_canonical_vecs_brute_force(A_eplitz + delta*eye(2*L+1),2*L+1);

%% equispaced sampled basis matrix
t_nyq = [0:N-1]/fs; % nyquist samples for testing
psi = exp(1i*2*pi*t_nyq'*freq_samples);

%% simulation specifications
trials = 2;
SNRs = -30:10:100;
SNR_fbst = zeros(trials,length(SNRs));
SNR_ttd32 = zeros(trials,length(SNRs));
SNR_ttd = zeros(trials, length(SNRs));
SNR_fbst_superfast_toeplitz = zeros(trials, length(SNRs));
SNR_fds = zeros(trials, length(SNRs));
SNR_fds32 = zeros(trials, length(SNRs));

for ii = 1:trials
    for jj = 1:length(SNRs)

        SNR_dB = SNRs(jj); % signal to noise ratio
        a = (randn(n_sinusoids,1) + 1i*randn(n_sinusoids,1))/sqrt(2); % sinusoid amplitudes
        f = W*(2*rand(n_sinusoids,1)-1); % sinusoid frequencies
        s = sig_gen(a,f,n_sinusoids,t);
        sigma = (norm(s)/sqrt(MN))*10^(-SNR_dB/20);
        noise = sigma*(randn(MN,1)+1i*randn(MN,1))/sqrt(2);
        SNR_check = db(norm(s)/norm(noise));
        fprintf('Trial: %d, Set SNR: %.2f, Measured SNR: %.2f\n',ii,SNR_dB,SNR_check)

        s_demod = s.*demod_phase;

        y = s_demod + noise;
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

        [y_nyq, w_l, alpha] = fbst_manual(y, L, B, A_eplitz, psi, idx, fc,Ws, T, delta,exp_fac);

        % alpha_fbst = inv(A_eplitz + delta*eye(2*L))*F_basis'*y(:);
        % y_nyq = psi*alpha_fbst;

        %% superfast toeplitz using canonical vectors
        b1 = toep_mult(X2,w_l);
        b2 = toep_mult(X4,w_l);
        b3 = toep_mult(X1,b1);
        b4 = toep_mult(X3,b2);
        alpha_superfast = b3 - b4;
        y_nyq_superfast = psi*alpha_superfast;
        % coeff_nyq = exp(-1i*pi*linspace(0,N-1,N)).';
        % y_nyq_transform = czt(alpha_superfast,N,exp(1i*pi/L),1);
        % y_nyq_transform = y_nyq_transform.*coeff_nyq;

        
        % tic;
        % for bb=1:B
        %     [y_nyq, w_l] = fbst_0(y, L, B, toeplitz_inv, psi, idx, fc, T);
        % end
        % runtime_fbst0 = toc;

        s_nyq = sig_gen(a,f,n_sinusoids,t_nyq');
        s_nyq_2 = sig_gen(a,f,n_sinusoids,t_nyq+1/fs);
        s_nyq_3 = sig_gen(a,f,n_sinusoids,t_nyq+S/fs);

        y_filter_ttd = TTD_Beamformer((y.*mod_phase).',-tau, fs, @Psi_sinc, 'R',taps).';
        y_bf_ttd = sum(y_filter_ttd(:,1:N),1)/M;

        y_filter_ttd32 = TTD_Beamformer((y.*mod_phase).',-tau, fs, @Psi_sinc, 'R',taps2).';
        y_bf_ttd32 = sum(y_filter_ttd32(:,1:N),1)/M;

        % %% fast delay and sum beamforming 16 taps
        % v_now = reshape(y.*mod_phase, [M,N,1]);
        % for stage=1:S
        %     v_prev = v_now;
        %     num_virtual_sensors = M/2^stage;
        %     num_sectors = R_fac^stage;
        %     v_now = zeros(num_virtual_sensors,N,num_sectors); % init partial beamforming at stage s
        % 
        %     for mm2=0:num_virtual_sensors-1
        %         for rr = 0:num_sectors-1
        %             m2 = mm2+1;
        %             r = floor((rr)/R_fac)+1;
        %             delays = ((-1)^(rr))*(1/(R_fac^stage))*eta*[2^(stage-1)*(2*mm2+0.5)-M/2; 2^(stage-1)*(2*mm2+1.5)-M/2];
        %             % v_prev_interped = TTD_Beamformer(y_interp.',-delays, fs, @Psi_sinc, 'R',taps2).';
        %             y_interp = reshape(v_prev([2*mm2+1,2*mm2+2],:,r),[2,N]);
        %             v_prev_interped = TTD_Beamformer(y_interp.',-delays, fs, @Psi_sinc, 'R',taps).';
        % 
        %             v_now(m2,:,rr+1) = sum(v_prev_interped,1)./2;
        %         end
        %     end
        % end
        % v_pred = v_now(1,:,idx_fds);
        % SNR_fds(ii,jj) = norm(s_nyq_3(comp_idxs))/norm(v_pred(comp_idxs) - s_nyq_3(comp_idxs));
        % 
        % 
        % %% fast delay and sum beamforming 32 taps
        % v_now = reshape(y.*mod_phase, [M,N,1]);
        % for stage=1:S
        %     v_prev = v_now;
        %     num_virtual_sensors = M/2^stage;
        %     num_sectors = R_fac^stage;
        %     v_now = zeros(num_virtual_sensors,N,num_sectors); % init partial beamforming at stage s
        % 
        %     for mm2=0:num_virtual_sensors-1
        %         for rr = 0:num_sectors-1
        %             m2 = mm2+1;
        %             r = floor((rr)/R_fac)+1;
        %             delays = ((-1)^(rr))*(1/(R_fac^stage))*eta*[2^(stage-1)*(2*mm2+0.5)-M/2; 2^(stage-1)*(2*mm2+1.5)-M/2];
        %             y_interp = reshape(v_prev([2*mm2+1,2*mm2+2],:,r),[2,N]);
        %             v_prev_interped = TTD_Beamformer(y_interp.',-delays, fs, @Psi_sinc, 'R',taps2).';
        %             v_now(m2,:,rr+1) = sum(v_prev_interped,1)./2;
        %         end
        %     end
        % end
        % v_pred = v_now(1,:,idx_fds);
        % SNR_fds32(ii,jj) = norm(s_nyq_3(comp_idxs))/norm(v_pred(comp_idxs) - s_nyq_3(comp_idxs));
        % 


        SNR_ttd(ii,jj) = norm(s_nyq_2(comp_idxs))/norm(y_bf_ttd(comp_idxs) - s_nyq_2(comp_idxs));
        SNR_ttd32(ii,jj) = norm(s_nyq_2(comp_idxs2))/norm(y_bf_ttd32(comp_idxs2) - s_nyq_2(comp_idxs2));
        SNR_fbst(ii,jj) = norm(s_nyq)/norm(s_nyq - y_nyq);
        SNR_fbst_superfast_toeplitz(ii,jj) = norm(s_nyq)/norm(s_nyq - y_nyq_superfast);

    end
end

figure(1)
plot(SNRs,SNRs + db(M)/2,'LineWidth',1)
hold on
grid on
plot(SNRs,db(mean(SNR_fbst)),'-*','LineWidth',1)
hold on
grid on
plot(SNRs, db(mean(SNR_ttd)),'-^','LineWidth',1)
hold on
grid on
plot(SNRs, db(mean(SNR_ttd32)),'-^','LineWidth',1)
hold on
grid on
plot(SNRs, db(mean(SNR_fds)),'-^','LineWidth',1)
hold on
grid on
plot(SNRs, db(mean(SNR_fds32)),'-^','LineWidth',1)
xlabel('Nominal SNR (dB)','Interpreter','latex','FontSize',12)
ylabel('Beamformed SNR (dB)','Interpreter','latex','FontSize',12)
legend({'Ideal','FBST', 'Delay and Sum (R=16)', 'Delay and Sum (R=32)', 'Fast Delay and Sum (R=16)','Fast Delay and Sum (R=32)'},'Location','northwest','Interpreter','latex','FontSize',12)
% exportgraphics(gcf, 'ula_snr.pdf', 'ContentType', 'vector');


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