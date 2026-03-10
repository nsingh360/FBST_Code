clear
clc
close all

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

B = 2^7; % no. of beams
M = 2^7; % array size
N = 2^6; % no. of samples
L = ceil(N); % no. of samples in frequency is 2*L
MN = M*N;

norm_fac = 2*L*M; % normalizing factor for topelitz matrix

Theta = zeros(2*B+1,1);
for i=1:2*B+1
    Theta(i) = asin((i-1-B)/B);
end

%% select angle to super resolve 
idx_1 = ceil(1.5*B);
idx_2 = ceil(1.5*B)+1;
alpha = 0.5;
theta = asin(alpha*sin(Theta(idx_1)) + (1-alpha)*sin(Theta(idx_2)));

fc = 20e9;  % center frequency
Ws = 5e9; % to prevent spatial aliasing
c = physconst('LightSpeed');
lambda = c/(fc+Ws); % wavelngth
W = 5e9; % bandwidth
fs = 2*W; % sampling frequency
T = 1/fs;
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

%% construct the toeplitz matrix
delta = 1e-5;
frequency_samples = (1/(2*T*L))*linspace(-L,L-1,2*L);
F = exp(1i*2*pi*t*frequency_samples);
F = F.*demod_phase;
A_eplitz = F'*F;
toeplitz_inv = inv(A_eplitz + (delta)*eye(2*L));

%% equispaced sampled basis matrix
t_nyq = [0:N-1]/fs; % nyquist samples for testing
psi = exp(1i*2*pi*t_nyq.'*frequency_samples);

%% simulation specification
trials = 5;
SNRs = -30:6:30;

nmse_beam = zeros(trials,length(SNRs));
SNR_true_beam = zeros(trials, length(SNRs));
SNR_interp_beam = zeros(trials, length(SNRs));

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
        y_array = reshape(y, [M,N]);

        %% chirp-z transform
        F_transform = zeros(M,2*L);
        A = -1;%exp(-1i*2*pi*T*W);
        W_czt = exp(-1i*pi/L);
        for mm=1:M
            F_transform(mm,:) = czt(y_array(mm,:),2*L,W_czt,A);
        end
        X_czt = zeros(2*L,2*B+1);
        A= 1;
        B_list = reshape(linspace(0,2*B,2*B+1), [1,2*B+1]);
        for ll=1:2*L
            l_dash = ll-1;
            A = exp(1i*pi*(1/B)*(fc/(fc+Ws) + (1/(2*(fc+Ws)*T*L))*(l_dash-L))*B);
            W_czt = exp(1i*pi*(1/B)*(fc/(fc+Ws) + (1/(2*(fc+Ws)*T*L))*(l_dash-L)));
            coeff1 = exp(-1i*pi*(1/B)*(fc/(fc+Ws) + (1/(2*(fc+Ws)*T*L))*(l_dash-L))*(M/2-1/2)*B_list); % adjusting for array center
            coeff2 = exp(1i*pi*(1/B)*(fc/(fc+Ws) + (1/(2*(fc+Ws)*T*L))*(l_dash-L))*(M/2-1/2)*B); % adjusting for negative angles
            X_czt(ll,:) = czt(F_transform(:,ll),2*B+1,W_czt,A);
            X_czt(ll,:) = coeff2*X_czt(ll,:).*coeff1;
        end

        w_l_interp = interp1(Theta(idx_1-3:idx_2+3),X_czt(:,idx_1-3:idx_2+3).',theta,'spline').';
        %% extract the true beam in theta direction
        w_l_true = F'*y;

        nmse_beam(ii,jj) = norm(w_l_true - w_l_interp)/norm(w_l_true);

        s_nyq = sig_gen(a,f,n_sinusoids,t_nyq');
        y_nyq_true = psi*toeplitz_inv*w_l_true;
        y_nyq_interp = psi*toeplitz_inv*w_l_interp;

        SNR_true_beam(ii,jj) = norm(s_nyq)/norm(y_nyq_true - s_nyq);
        SNR_interp_beam(ii,jj) = norm(s_nyq)/norm(y_nyq_interp - s_nyq);

    end
end

figure(1)
plot(SNRs,mean(nmse_beam),'LineWidth',1)
xlabel('Nominal SNR (dB)','Interpreter','latex')
ylabel('NMSE (dB)','Interpreter','latex')

figure(2)
plot(SNRs,SNRs + db(M)/2,'LineWidth',1)
hold on
grid on
plot(SNRs,db(mean(SNR_true_beam)),'-^','LineWidth',1)
hold on
grid on
plot(SNRs, db(mean(SNR_interp_beam)),'-^','LineWidth',1)
xlabel('Nominal SNR (dB)','Interpreter','latex')
ylabel('Beamformed SNR (dB)','Interpreter','latex')
legend({'Ideal','True Beam','Interpolated Beam'},'Location','northwest','Interpreter','latex')

%% supporting functions
function [X] = sig_gen(a,f,n_sinusoids,t)
X = zeros(size(t));
for ii = 1:n_sinusoids % For each antenna
    sig = (a(ii)*exp(1i*2*pi*(f(ii))*(t)));
    X = X + sig;
end
end
