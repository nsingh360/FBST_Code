clear
clc
close all

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

B = 2^6; % no. of beams
M = 2^6; % array size
N = 2^5; % no. of samples
L = ceil(N); % no. of samples in frequency is 2*L
exp_fac = 1.05; % factor to include freqs slightly out of band;
MN = M*N;

Theta = zeros(2*B+1,1);
idx = 25;%randi(B,1);
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
idx_fds = 12;

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
fc_tilde = (fc+Ws); % adjustment for spatial aliasing
fc_norm = fc/fc_tilde;

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

%% Interferer Specs
%% select angle to super resolve 
idx_1 = 90;%61;%239;
alpha = 1;
theta_i = Theta(idx_1);
% idx_i = B;
% theta_i = Theta(idx_i);
u_i = sin(theta_i); % normal vector for determining delays
tau_i = x_pos*u_i/c;
t_i_array = n - tau_i;
t_i = t_i_array(:);
T_aperture_i = max(t_i) - min(t_i);
demod_phase_i  = repmat(exp(-1i*2*pi*fc*tau_i),N,1); % modulation vector after pre-steering
mod_phase_i = exp(1i*2*pi*(fc)*tau_i);

%% create neighboring angles for beamspace interference
theta_i_low = theta_i - pi/80;
u_i_low = sin(theta_i_low); % normal vector for determining delays
tau_i_low = x_pos*u_i_low/c;
t_i_array_low = n - tau_i_low;
t_i_low = t_i_array_low(:);
T_aperture_i_low = max(t_i_low) - min(t_i_low);
demod_phase_i_low  = repmat(exp(-1i*2*pi*fc*tau_i_low),N,1); % modulation vector after pre-steering

theta_i_high = theta_i + pi/80;
u_i_high = sin(theta_i_high); % normal vector for determining delays
tau_i_high = x_pos*u_i_high/c;
t_i_array_high = n - tau_i_high;
t_i_high = t_i_array_high(:);
T_aperture_i_high = max(t_i_high) - min(t_i_high);
demod_phase_i_high  = repmat(exp(-1i*2*pi*fc*tau_i_high),N,1); % modulation vector after pre-steering

%% correlation matrix
n_sig = [-2*N:4*N]/(2*W);
e_s  = repmat(exp(1i*2*pi*fc*(-tau)),N,1);
Bs = bsxfun(@times,e_s,sinc(2*W*(t-n_sig)));
Bs_scale = 1/sqrt(trace(Bs'*Bs));
Bs = Bs*Bs_scale;
Rs = (Bs*Bs');

e_i  = repmat(exp(1i*2*pi*fc*(-tau_i)),N,1);
Bi = bsxfun(@times,e_i,sinc(2*W*(t_i-n_sig)));
Bi_scale = 1/sqrt(trace(Bi'*Bi));
Bi = Bi*Bi_scale;

%% construct the toeplitz matrix
delta = 1e-5;
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

F_basis_i = exp(1i*2*pi*t_i*freq_samples).*demod_phase_i;

F_basis_i_low = exp(1i*2*pi*t_i_low*freq_samples).*demod_phase_i_low;
F_basis_i_high = exp(1i*2*pi*t_i_high*freq_samples).*demod_phase_i_high;
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
V_nyq = interp1(t_u,V_u,MN*(t_nyq-min(t))/T_aperture,'pchip',0); % interpolate to nyquist
psi = exp(1i*2*pi*t_nyq'*freq_samples);
psi2 = exp(1i*2*pi*t_nyq'*freq_samples2);


%% simulation specifications
trials = 2;
SIRs = -30;
SNRs = -30:6:30;

SNR_fbst_mvdr = zeros(trials,length(SNRs));
SNR_fbst_mvdr_sparse = zeros(trials, length(SNRs));
SNR_fbst_mvdr_sparse_source = zeros(trials, length(SNRs));

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
        y_array = reshape(y,[M,N]);
        sigma_i = 10^(-SIR_dB/20);
        Ri = (sigma_i^2)*(Bi*Bi');
        R = Rs + Ri + sigma^2/MN*eye(MN);
        Rinv = R^(-1);

        % %% FBST MVDR
        % F_transform = zeros(M,2*L+1);
        % A = exp(-1i*pi*exp_fac);%exp(-1i*2*pi*T*W);
        % W_czt = exp(-1i*pi*exp_fac/L);
        % for mm=1:M
        %     F_transform(mm,:) = czt(y_array(mm,:).',2*L+1,W_czt,A);
        % end
        % X_czt = zeros(2*L+1,2*B+1);
        % A= 1;
        % B_list = reshape(linspace(0,2*B,2*B+1), [1,2*B+1]);
        % for ll=1:2*L+1
        %     l_dash = ll-1;
        %     A = exp(1i*pi*(1/B)*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L))*B);
        %     W_czt = exp(1i*pi*(1/B)*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L)));
        %     coeff1 = exp(-1i*pi*(1/B)*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L))*(M/2-1/2)*B_list); % adjusting for array center
        %     coeff2 = exp(1i*pi*(1/B)*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L))*(M/2-1/2)*B); % adjusting for negative angles
        %     X_czt(ll,:) = czt(F_transform(:,ll),2*B+1,W_czt,A);
        %     X_czt(ll,:) = coeff2*X_czt(ll,:).*coeff1;
        % end
        % w_s = X_czt(:,idx);
        beam_space_map = [F_basis'; F_basis_i_low'; F_basis_i_high'];
        beam_space_map2 = [F_basis'];
        beam_space_map = ((beam_space_map*beam_space_map')^(-0.5))*beam_space_map;
        beam_space_map2 = ((beam_space_map2*beam_space_map2')^(-0.5))*beam_space_map2;
        F_square = beam_space_map*F_basis;
        F_square2 = beam_space_map2*F_basis;
        w_s = beam_space_map*y;
        w_s2 = beam_space_map2*y;

        Rsource = beam_space_map*R*beam_space_map';
        Rs_inv = inv(Rsource + 1e-5*eye(195));

        Rsource2 = beam_space_map2*R*beam_space_map2';
        Rs_inv2 = inv(Rsource2 + 1e-5*eye(2*L+1));

        s_nyq = sig_gen(a,f,n_sinusoids,t_nyq')*s_scale;
        inv_m = inv(F_basis'*Rinv*F_basis + 1e-4*eye(2*L+1));
        alpha = inv_m*F_basis'*Rinv*y;
        y_nyq = psi*alpha;

        alpha_m = inv(F_square'*Rs_inv*F_square + 1e-4*eye(2*L+1))*F_square'*Rs_inv*w_s;
        y_nyq_m = psi*alpha_m;

        alpha_m2 = inv(F_square2'*Rs_inv2*F_square2 + 1e-4*eye(2*L+1))*F_square2'*Rs_inv2*w_s2;
        y_nyq_m2 = psi*alpha_m2;

        SNR_fbst_mvdr(ii,jj) = norm(s_nyq)/norm(s_nyq - y_nyq);
        SNR_fbst_mvdr_sparse(ii,jj) = norm(s_nyq)/norm(s_nyq - y_nyq_m);
        SNR_fbst_mvdr_sparse_source(ii,jj) = norm(s_nyq)/norm(s_nyq - y_nyq_m2);
    end
    end
end

figure(1)
plot(SNRs,SNRs + db(M)/2,'LineWidth',1)
hold on
grid on
plot(SNRs,db(mean(SNR_fbst_mvdr)),'-*','LineWidth',1)
hold on
grid on
plot(SNRs,db(mean(SNR_fbst_mvdr_sparse)),'-*','LineWidth',1)
hold on
grid on
plot(SNRs, db(mean(SNR_fbst_mvdr_sparse_source)),'-*','LineWidth',1);
xlabel('Nominal SNR (dB)','Interpreter','latex','FontSize',12)
ylabel('Beamformed SNR (dB)','Interpreter','latex','FontSize',12)
legend({'Ideal','FBST MVDR','FBST MVDR Sparse', 'FBST MVDR Sparse source'},'Location','northwest','Interpreter','latex','FontSize',12)
% exportgraphics(gcf, 'ula_snr_new.pdf', 'ContentType', 'vector');


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