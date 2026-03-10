clear
clc
close all

%% set random number generator seed
% digits(8);
seed = randi([1,10000]);
rng(seed)

M = 2^7; % array size
N = 2^6; % no. of samples
L = ceil(N);
MN = M*N;
S = log2(M);
R= 2;


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
idx = 15;
theta = Theta_fds(idx);
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
idx_1 = 30;%61;%239;
idx_2 = 240;
alpha = 0.5;
theta_i = Theta_fds(idx_1);%asin(alpha*sin(Theta(idx_1)) + (1-alpha)*sin(Theta(idx_2)));
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
F = F.*demod_phase;
A_eplitz = F'*F;
toeplitz_inv = inv(A_eplitz + (delta)*eye(2*L));

%% construct fourier and toeptliz matrix for interferer
F_i = exp(1i*2*pi*t_i*freq_samples);
A_eplitz_i = F_i'*F_i;
F_i = F_i.*demod_phase_i;

t_nyq = [0:N-1]/fs;
t_nyq2 = t_nyq + 1/fs;
t_nyq3 = t_nyq + S/(fs);
psi_ds = exp(1i*2*pi*t_nyq2'*freq_samples);
psi_fds = exp(1i*2*pi*t_nyq3'*freq_samples);

%% simulation specifications
trials = 5;
SIRs = -20;
SNRs = -30:6:30;

SNR_ttd = zeros(trials, length(SNRs));
SNR_ttd_nulled = zeros(trials, length(SNRs));
SNR_ttd_nulled_beamspace = zeros(trials, length(SNRs));
SNR_fds = zeros(trials, length(SNRs));
SNR_fds_nulled = zeros(trials, length(SNRs));
SNR_fds_nulled_beamspace = zeros(trials, length(SNRs));

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
        y_array = reshape(y, [M,N]);

        y_i_nulled = y - F_i*inv(A_eplitz_i + delta*eye(2*L))*F_i'*y;
        % P_null = eye(MN,MN) - F_i*inv(A_eplitz_i + delta*eye(2*L))*F_i';

        y_array_nulled = reshape(y_i_nulled, [M,N]); % used for ttd

        s_nyq2 = sig_gen(a,f,n_sinusoids,t_nyq2')*s_scale;
        s_nyq3 = sig_gen(a,f,n_sinusoids,t_nyq3')*s_scale;

        %% TTD
        y_filter_ttd = TTD_Beamformer((y_array.*mod_phase).',-tau, fs, @Psi_sinc, 'R',taps).';
        y_bf_ttd = sum(y_filter_ttd(:,1:N),1)/M;
        SNR_ttd(ii,jj) = norm(s_nyq2(comp_idxs).')/norm(y_bf_ttd(comp_idxs) - s_nyq2(comp_idxs).');

        % %% TTD interferer
        % y_filter_ttd_i = TTD_Beamformer((y_array.*mod_phase_i).',-tau_i, fs, @Psi_sinc, 'R',taps).';
        % y_bf_ttd_i = sum(y_filter_ttd_i(:,1:N),1)/M;

        y_filter_ttd_nulled = TTD_Beamformer((y_array_nulled.*mod_phase).',-tau, fs, @Psi_sinc, 'R',taps).';
        y_bf_ttd_nulled = sum(y_filter_ttd_nulled(:,1:N),1)/M;
        SNR_ttd_nulled(ii,jj) = norm(s_nyq2(comp_idxs).')/norm(y_bf_ttd_nulled(comp_idxs) - s_nyq2(comp_idxs).');

        %% Fast delay and sum (for source nulled)
        v_now = y_array_nulled.*mod_phase;
        for stage=1:S
            v_prev = v_now;
            num_virtual_sensors = M/2^stage;
            num_sectors = R^stage;
            v_now = zeros(num_virtual_sensors,N,num_sectors); % init partial beamforming at stage s

            for mm=0:num_virtual_sensors-1
                for rr = 0:num_sectors-1
                    m = mm+1;
                    r = floor((rr)/R)+1;
                    delays = ((-1)^(rr))*(1/(R^stage))*eta*[2^(stage-1)*(2*mm+0.5)-M/2; 2^(stage-1)*(2*mm+1.5)-M/2];
                    y_interp = reshape(v_prev([2*mm+1,2*mm+2],:,r),[2,N]);
                    v_prev_interped = TTD_Beamformer(y_interp.',-delays, fs, @Psi_sinc, 'R',taps).';
                    v_now(m,:,rr+1) = sum(v_prev_interped,1)./2;
                end
            end
        end
        % v_now = squeeze(v_now);
        y_fds_s_nulled = v_now(1,:,idx).';

        %% fast delay and sum (for source beam)
        v_now = y_array.*mod_phase;
        for stage=1:S
            v_prev = v_now;
            num_virtual_sensors = M/2^stage;
            num_sectors = R^stage;
            v_now = zeros(num_virtual_sensors,N,num_sectors); % init partial beamforming at stage s

            for mm=0:num_virtual_sensors-1
                for rr = 0:num_sectors-1
                    m = mm+1;
                    r = floor((rr)/R)+1;
                    delays = ((-1)^(rr))*(1/(R^stage))*eta*[2^(stage-1)*(2*mm+0.5)-M/2; 2^(stage-1)*(2*mm+1.5)-M/2];
                    y_interp = reshape(v_prev([2*mm+1,2*mm+2],:,r),[2,N]);
                    v_prev_interped = TTD_Beamformer(y_interp.',-delays, fs, @Psi_sinc, 'R',taps).';
                    v_now(m,:,rr+1) = sum(v_prev_interped,1)./2;
                end
            end
        end
        % v_now = squeeze(v_now);
        y_fds_s = v_now(1,:,idx).';

        %% fast delay and sum (for interferer beam)
        v_now = y_array.*mod_phase_i;
        for stage=1:S
            v_prev = v_now;
            num_virtual_sensors = M/2^stage;
            num_sectors = R^stage;
            v_now = zeros(num_virtual_sensors,N,num_sectors); % init partial beamforming at stage s

            for mm=0:num_virtual_sensors-1
                for rr = 0:num_sectors-1
                    m = mm+1;
                    r = floor((rr)/R)+1;
                    delays = ((-1)^(rr))*(1/(R^stage))*eta*[2^(stage-1)*(2*mm+0.5)-M/2; 2^(stage-1)*(2*mm+1.5)-M/2];
                    y_interp = reshape(v_prev([2*mm+1,2*mm+2],:,r),[2,N]);
                    v_prev_interped = TTD_Beamformer(y_interp.',-delays, fs, @Psi_sinc, 'R',taps).';
                    v_now(m,:,rr+1) = sum(v_prev_interped,1)./2;
                end
            end
        end
        % v_now = squeeze(v_now);
        y_fds_i = v_now(1,:,idx_1).';

        y_filter_ttd_nulled = TTD_Beamformer((y_array.*mod_phase).',-tau, fs, @Psi_sinc, 'R',taps).';
        y_ds_s = sum(y_filter_ttd_nulled(:,1:N),1)/M;

        y_filter_ttd_nulled = TTD_Beamformer((y_array.*mod_phase_i).',-tau_i, fs, @Psi_sinc, 'R',taps).';
        y_ds_i = sum(y_filter_ttd_nulled(:,1:N),1)/M;

        

        G = toeplitz_inv*F'*F_i;
        G_fds = psi_fds*G*pinv(psi_fds);
        G_ds = psi_ds*G*pinv(psi_ds);
        y_ds_nyq = y_ds_s.' - G_ds*(y_ds_i.');

        y_fds_nyq = y_fds_s - G_fds*(y_fds_i); 
        SNR_fds(ii,jj) = norm(s_nyq3(comp_idxs))/norm(y_fds_s(comp_idxs) - s_nyq3(comp_idxs));
        SNR_fds_nulled(ii,jj) = norm(s_nyq3(comp_idxs))/norm(y_fds_s_nulled(comp_idxs) - s_nyq3(comp_idxs));
        SNR_fds_nulled_beamspace(ii,jj) = norm(s_nyq3(comp_idxs))/norm(y_fds_nyq(comp_idxs) - s_nyq3(comp_idxs));
        SNR_ttd_nulled_beamspace(ii,jj) = norm(s_nyq2(comp_idxs))/norm(y_ds_nyq(comp_idxs) - s_nyq2(comp_idxs));

    end
    end
end


figure(1)
plot(SNRs,SNRs + db(M)/2,'LineWidth',1)
hold on
grid on
plot(SNRs,db(mean(SNR_ttd)),'-^','LineWidth',1)
hold on
grid on
% plot(SNRs, db(mean(SNR_ttd_nulled)),'-^','LineWidth',1)
% hold on
% grid on
plot(SNRs, db(mean(SNR_ttd_nulled_beamspace)),'-^','LineWidth',1)
hold on
grid on
plot(SNRs,db(mean(SNR_fds)),'-^','LineWidth',1)
hold on
grid on
% plot(SNRs, db(mean(SNR_fds_nulled)),'-^','LineWidth',1)
% hold on
% grid on
plot(SNRs, db(mean(SNR_fds_nulled_beamspace)),'-^','LineWidth',1)
xlabel('Nominal SNR (dB)','Interpreter','latex')
ylabel('Beamformed SNR (dB)','Interpreter','latex')
legend({'Ideal', 'TTD','TTD Nulled beamspace','FDS', 'FDS Nulled beamspace'},'Location','northwest','Interpreter','latex')
% exportgraphics(gcf, 'ula_snr_interference.pdf', 'ContentType', 'vector');

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