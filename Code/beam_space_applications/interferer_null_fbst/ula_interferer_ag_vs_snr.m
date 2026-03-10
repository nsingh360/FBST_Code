clear
clc
close all

%% set random number generator seed
% digits(8);
seed = randi([1,10000]);
rng(seed)

B = 2^7; % no. of beams
M = 2^7; % array size
N = 2^6; % no. of samples
L = ceil(N); % no. of samples in frequency is 2*L
MN = M*N;
S = log2(M);
R= 2;

norm_fac = 2*L*M; % normalizing factor for topelitz matrix

Theta = zeros(2*B+1,1);
idx = 2*B+1;%120;%randi(B,1);
for i=1:2*B+1
    Theta(i) = asin((i-1-B)/B);
end

%% Angle sampling for final stage
r_lin = linspace(0,R^S-1,R^S);
Theta_fds = asin(1 - (2*r_lin + 1)/R^S);

fc = 20e9;  % center frequency
Ws = 5e9; % to prevent spatial aliasing
c = physconst('LightSpeed');
lambda = c/(fc+Ws); % wavelngth
W = 5e9; % bandwidth
fs = 2*W; % sampling frequency
T = 1/fs;
theta = Theta(idx);
m = [(-M/2+1/2):(M/2-1/2)]';
x_pos = m*lambda/2; % sensor positions
u_s = sin(theta); % normal vector for determining delays
fc_tilde = (fc+Ws); % adjustment for spatial aliasing
fc_norm = fc/fc_tilde;
eta = lambda/(2*c);

%% Signal specs that do not need to be redifined in loop
n_sinusoids = 100; % number of sinusoids in signal
n = [0:N-1]/fs; % temporal sample vectors
tau = x_pos*u_s/c; % relative delays to phase center
t_array = n-tau; % delays across array and time
t = t_array(:); %
T_aperture = max(t)-min(t);% temporal aperture of the array
demod_phase = repmat(exp(-1i*2*pi*fc*tau),N,1);
mod_phase = exp(1i*2*pi*(fc)*tau);

%% Interferer Specs
%% select angle to super resolve 
idx_1 = 239;%61;%239;
idx_2 = 240;
alpha = 0.5;
theta_i = asin(alpha*sin(Theta(idx_1)) + (1-alpha)*sin(Theta(idx_2)));
% idx_i = B;
% theta_i = Theta(idx_i);
u_i = sin(theta_i); % normal vector for determining delays
tau_i = x_pos*u_i/c;
t_i_array = n - tau_i;
t_i = t_i_array(:);
T_aperture_i = max(t_i) - min(t_i);
demod_phase_i  = repmat(exp(-1i*2*pi*fc*tau_i),N,1); % modulation vector after pre-steering
mod_phase_i = exp(1i*2*pi*(fc)*tau_i);

taps = 8;
edge_lim = taps + ceil(fs*(max(tau)-min(tau))/2);% clip signal to account for edge distortion
comp_idxs = (edge_lim+1):(N-edge_lim-1); % comparison indicies

taps2 = 16;
edge_lim2 = taps2 + ceil(fs*(max(tau)-min(tau))/2);% clip signal to account for edge distortion
comp_idxs2 = (edge_lim2+1):(N-edge_lim2-1); % comparison indicies

%% construct the toeplitz matrix
delta = 1e-5;
freq_samples = 1/(2*T*L)*linspace(-L,L-1,2*L);
F = exp(1i*2*pi*t*freq_samples);
A_eplitz = F'*F;
toeplitz_inv = inv(A_eplitz + (delta)*eye(2*L));
F = F.*demod_phase;

%% construct fourier and toeptliz matrix for interferer
F_i = exp(1i*2*pi*t_i*freq_samples);
A_eplitz_i = F_i'*F_i;
F_i = F_i.*demod_phase_i;

%% get canonical vectors of Toeplitz inverse via brute force
[X1,X2,X3,X4] = eval_canonical_vecs_brute_force(A_eplitz + delta*eye(2*L),2*L);

%% equispaced sampled basis matrix
t_nyq = min(t):T/1.1:max(t); % nyquist samples for testing
psi = exp(1i*2*pi*t_nyq'*freq_samples);
t_nyq3 = [0:N-1]/fs;
psi_nyq = exp(1i*2*pi*t_nyq3'*freq_samples);


%% simulation specifications
trials = 2;
SIRs = -10;
SNRs = 30;

SNR_fbst_nulled = zeros(trials,length(SNRs));
SNR_fbst_not_nulled = zeros(trials, length(SNRs));
SNR_fbst_no_i = zeros(trials, length(SNRs));
SNR_ttd = zeros(trials, length(SNRs));
SNR_ttd_nulled = zeros(trials, length(SNRs));
alpha_nmse = zeros(trials, length(SNRs));

for jj = 1:length(SNRs)
    for kk = 1:length(SIRs)
        SIR_dB = SIRs(kk); % signal to interferer ratio
        SNR_dB = SNRs(jj); % signal to noise ratio
        sigma = 10^(-SNR_dB/20);
    for ii = 1:trials
        SNR_dB = SNRs(jj); % signal to noise ratio
        fprintf("Trials %f SNR %f SIR %f:\n ", ii, SNR_dB, SIR_dB)
        a = (randn(n_sinusoids,1) + 1i*randn(n_sinusoids,1))/sqrt(2); % sinusoid amplitudes
        f = W*(2*rand(n_sinusoids,1)-1); % sinusoid frequencies
        s = sig_gen(a,f,n_sinusoids,t);
        s_scale = 1/rms(s);
        s_demod = s_scale*s.*demod_phase;

        a_i = (randn(n_sinusoids,1) + 1i*randn(n_sinusoids,1))/sqrt(2); % random sinusoid coefficients
        f_i = W*(2*rand(n_sinusoids,1)-1); % random sinusoid freqs
        s_i = demod_phase_i.*sig_gen(a_i,f_i,n_sinusoids,t_i);
        s_i_scale = 1/rms(s_i);
        s_i_demod = s_i_scale*s_i*10^(-SIR_dB/20);
        

        noise = sigma*(randn(MN,1)+1i*randn(MN,1))/sqrt(2);

        y = s_demod + s_i_demod + noise;
        y_no_i = s_demod+noise;
        % w_l_i_check = F_i'*y;
        y_array = reshape(y, [M,N]);
        y_array_no_i = reshape(y_no_i,[M,N]);

        %% interferer nulling
        % [~,w_l_s,~] = fbst_1(y_array, L, B, A_eplitz_i, psi, idx, fc,Ws, T, delta);
        % [~, w_l_i,~] = fbst_1(y_array, L, B, A_eplitz_i, psi, idx_i, fc,Ws, T, delta);

        %% measure run-times 
        y_i_nulled = y - F_i*inv(A_eplitz_i + delta*eye(2*L))*F_i'*y;
        P_null = eye(MN,MN) - F_i*inv(A_eplitz_i + delta*eye(2*L))*F_i';
        % ts_array_space = tic;
        % coeff_nulled = F'*(y - F_i*inv(A_eplitz_i + delta*eye(2*L))*F_i'*y);
        % runtime_array_space = toc(ts_array_space);
        % 
        % ts_beam_space = tic;
        % w_l_interp = interp1(rand(8,1),rand(64,8).',0.4,'spline').';
        % coeff_nulled_2 = w_l_s - F'*F_i*inv(A_eplitz_i + delta*eye(2*L))*w_l_i;
        % runtime_beam_space = toc(ts_beam_space);

        y_array_nulled = reshape(y_i_nulled, [M,N]); % used for ttd

        %% Nulling in beam space
        [M,N] = size(y_array);
        %--- chirp-z transform
        F_transform = zeros(M,2*L);
        A = -1;%exp(-1i*2*pi*T*W);
        W_czt = exp(-1i*pi/L);
        for mm=1:M
            F_transform(mm,:) = czt(y_array(mm,:).',2*L,W_czt,A);
        end
        X_czt = zeros(2*L,2*B+1);
        A= 1;
        B_list = reshape(linspace(0,2*B,2*B+1), [1,2*B+1]);
        for ll=1:2*L
            l_dash = ll-1;
            A = exp(1i*pi*(1/B)*(fc_norm + (1/(2*fc_tilde*T*L))*(l_dash-L))*B);
            W_czt = exp(1i*pi*(1/B)*(fc_norm + (1/(2*fc_tilde*T*L))*(l_dash-L)));
            coeff1 = exp(-1i*pi*(1/B)*(fc_norm + (1/(2*fc_tilde*T*L))*(l_dash-L))*(M/2-1/2)*B_list); % adjusting for array center
            coeff2 = exp(1i*pi*(1/B)*(fc_norm + (1/(2*fc_tilde*T*L))*(l_dash-L))*(M/2-1/2)*B); % adjusting for negative angles
            X_czt(ll,:) = czt(F_transform(:,ll),2*B+1,W_czt,A);
            X_czt(ll,:) = coeff2*X_czt(ll,:).*coeff1;
        end
        w_l = X_czt(:,idx);
        w_l_i = X_czt(:,idx_1);
        % % w_l_int = X_czt(:,idx_i);
        w_l_interp = interp1(Theta(idx_1-3:idx_2+3),X_czt(:,idx_1-3:idx_2+3).',theta_i,'spline').';
        % w_l_int2 = w_l_int - F_i'*F*inv(F'*F + delta*eye(2*L))*w_l;
        w_l_nulled = w_l - F'*F_i*inv(F_i'*F_i + delta*eye(2*L))*w_l_interp;

        alpha = toeplitz_inv*w_l_nulled;
        y_nyq = psi*alpha;

        beta_s = toeplitz_inv*w_l;
        beta_i = inv(F_i'*F_i + delta*eye(2*L))*w_l_interp;
        G = toeplitz_inv*F'*F_i;
        beta_tilde = beta_s - G*beta_i;
        G_ = psi*G*pinv(psi);
        % G_fds = psi_fds*G*pinv(psi_fds);
        y_nyq_time = psi*beta_tilde;

        %------ create interpolating vectors in time.
        theta_interps = Theta(idx_1-3:idx_2+3);
        y_nyqs_interps = zeros(length(theta_interps),length(t_nyq));
        count_indx = -3;
        for it=1:length(theta_interps)
            theta_interp = theta_interps(it);
            t_array_interp = n - x_pos*sin(theta_interp)/c;
            t_interp = t_array_interp(:);
            F_interp = exp(1i*2*pi*t_interp*freq_samples);
            y_nyqs_interps(it,:) = psi*inv(F_interp'*F_interp + delta*eye(2*L))*X_czt(:,idx_1+count_indx);
            count_indx = count_indx + 1;
        end

        y_nyq_s = psi*beta_s;
        y_nyq_i = psi*beta_i;
        y_nyq_i_interp = interp1(Theta(idx_1-3:idx_2+3),y_nyqs_interps,theta_i,'spline').';
        y_nyq_time2 = y_nyq_s - G_*y_nyq_i_interp;

        [y_nyq2, w_l2, alpha2] = fbst_1(y_array, L, B, A_eplitz, psi_nyq, idx, fc,Ws, T, delta);

        [y_nyq3, w_l3, alpha3] = fbst_1(y_array_no_i, L, B, A_eplitz, psi_nyq, idx, fc,Ws, T, delta);

        s_nyq = sig_gen(a,f,n_sinusoids,t_nyq')*s_scale;
        s_nyq3 = sig_gen(a,f,n_sinusoids,t_nyq3')*s_scale;

        % SNR_ttd(ii,jj) = norm(s_nyq_2(comp_idxs))/norm(y_bf_ttd(comp_idxs) - s_nyq_2(comp_idxs));
        % SNR_ttd32(ii,jj) = norm(s_nyq_2(comp_idxs2))/norm(y_bf_ttd32(comp_idxs2) - s_nyq_2(comp_idxs2));
        SNR_fbst_nulled(ii,jj) = norm(s_nyq(t_nyq<=(N-1)*T & t_nyq>=0))/norm(s_nyq(t_nyq<=(N-1)*T & t_nyq>=0) - y_nyq_time2(t_nyq<=(N-1)*T & t_nyq>=0));
        SNR_fbst_not_nulled(ii,jj) = norm(s_nyq3)/norm(s_nyq3 - y_nyq2);
        SNR_fbst_no_i(ii,jj) = norm(s_nyq3)/norm(s_nyq3 - y_nyq3);
        
        alpha_nmse(ii,jj) = norm(alpha3 - alpha)/norm(alpha);
    end
    end
end

figure(1)
plot(SNRs,SNRs + db(M)/2,'LineWidth',1)
hold on
grid on
plot(SNRs,db(mean(SNR_fbst_nulled)),'-^','LineWidth',1)
hold on
grid on
plot(SNRs, db(mean(SNR_fbst_not_nulled)),'-^','LineWidth',1)
hold on
grid on
plot(SNRs, db(mean(SNR_fbst_no_i)),'-^','LineWidth',1)
xlabel('Nominal SNR (dB)','Interpreter','latex')
ylabel('Beamformed SNR (dB)','Interpreter','latex')
legend({'Ideal', 'FBST Nulled', 'FBST not nulled', 'FBST pure'},'Location','northwest','Interpreter','latex')
exportgraphics(gcf, 'ula_snr_interference_fbst.pdf', 'ContentType', 'vector');
figure(2)
plot(SNRs, mean(alpha_nmse),'LineWidth',1)

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