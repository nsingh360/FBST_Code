clear
clc
close all

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

B = 2^7; % no. of beams
M = 2^6; % array size
N = 2^8; % no. of samples per batch
K = 2^4; % no. of batches
L = ceil(N); % no. of samples in frequency is 2*L
eta = 0.001; % relative transition width
MN = M*N;

Theta = zeros(B,1);
idx = 90;%randi(B,1);
for i=1:B
    Theta(i) = asin((i-1)/B);
end

fc = 20e9;  % center frequency
c = physconst('LightSpeed');
lambda = c/fc; % wavelngth
W = 5e9; % bandwidth
fs = 2.5*W; % sampling frequency
T = 1/fs;
theta = Theta(idx);
m = [(-M/2+1/2):(M/2-1/2)]';
x_pos = m*lambda/2; % sensor positions
u_s = sin(theta); % normal vector for determining delays

%% Signal specs that do not need to be redifined in loop
n_sinusoids = 20; % number of sinusoids in signal
n = [0:N-1]/fs; % temporal sample vectors
tau = x_pos*u_s/c; % relative delays to phase center
t_array = n-tau; % delays across array and time
t = t_array(:); %
T_aperture = max(t)-min(t);% temporal aperture of the array


batch_samples = cell(K,1);
Fourier_ext_basis = cell(K,1);
nyquist_samples = cell(K,1); % nyquist samples for testing_
N_per_batch = N/K; % samples per batch
L_per_batck = ceil(N_per_batch); %no. of basis is 2*L
F_L = (1/(2*T*L_per_batck))*linspace(-L_per_batck,L_per_batck-1,2*L_per_batck);
for k = 1:K
    t_k = nonzeros(t.*(t>=((k-1)*T_aperture/K+min(t))).*(t<=(k*T_aperture/K+min(t))));
    nyquist_samples{k} = linspace((k-1)*T_aperture/K+min(t), k*T_aperture/K+min(t), N_per_batch);%(k-1)*T_aperture/K + min(t) : 1/fs : k*T_aperture/K+min(t);
    batch_samples{k} = t_k;

    Fourier_ext_basis{k} = exp(1i*2*pi*t_k*F_L);
end


trials = 10;
SNRs = -30:6:30;
SNR_per_batch = zeros(trials, length(SNRs),K);

for ii=1:trials
    for jj=1:length(SNRs)
        SNR_dB = SNRs(jj); % signal to noise ratio
        fprintf('Trial: %d, Set SNR: %.2f\n',ii,SNR_dB)
        %% create K batches of observed sample
        Y = cell(K,1);
        S = cell(K,1);
        S_nyq = cell(K,1);
        a = (randn(n_sinusoids,1) + 1i*randn(n_sinusoids,1))/sqrt(2); % sinusoid amplitudes
        f = W*(2*rand(n_sinusoids,1)-1); % sinusoid frequencies
        for k=1:K
            samples_K = length(batch_samples{k});
            s = sig_gen(a,f,n_sinusoids,batch_samples{k});
            sigma = (norm(s)/sqrt(samples_K))*10^(-SNR_dB/20);
            noise = sigma*(randn(samples_K,1)+1i*randn(samples_K,1))/sqrt(2);
    
            s_demod = s + noise;
    
            y_k = s_demod;
            Y{k} = y_k;
            S{k} = s;
            S_nyq{k} = sig_gen(a,f,n_sinusoids,nyquist_samples{k});
        end

        %% K batch reconstruction
        w_per_batch = cell(K,1);
        for k=1:K
            A = Fourier_ext_basis{k};
            wk = inv(A'*A + 1e-5*eye(2*L_per_batck))*A'*Y{k};
            w_per_batch{k} = wk;
            SNR_per_batch(ii,jj,k) = norm(S_nyq{k}.')/norm(S_nyq{k}.' - exp(1i*2*pi*nyquist_samples{k}.'*F_L)*wk);
        end

    end
end

labels = cell(1,K+1);
labels{1} = 'ideal';
figure(1)
plot(SNRs,SNRs + db(M)/2,'LineWidth',1)
hold on
grid on
for k=1:K
    SNR_batch_K = SNR_per_batch(:,:,k);
    plot(SNRs,db(mean(SNR_batch_K)),'-^','LineWidth',1);
    hold on
    grid on
    labels{k+1} = ['Batch ' num2str(k)];
end
xlabel('Nominal SNR (dB)','Interpreter','latex')
ylabel('Beamformed SNR (dB)','Interpreter','latex')
legend(labels,'Location','northwest','Interpreter','latex')


%% supporting function
function [X] = sig_gen(a,f,n_sinusoids,t)
X = zeros(size(t));
for ii = 1:n_sinusoids % For each antenna
    sig = (a(ii)*exp(1i*2*pi*(f(ii))*(t)));
    X = X + sig;
end
end
