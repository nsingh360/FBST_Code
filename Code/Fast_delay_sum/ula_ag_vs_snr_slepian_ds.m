clear
clc
close all

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

%% array spatial and temporal specificiations
M = 2^7; % array size
N = 2^6; % temporal samples
K = 2^4; % filter taps
MN = M*N;
fc = 20e9;  % center frequency
c = physconst('LightSpeed');
W = 5e9; % bandwidth
lambda = c/(fc+W); % wavelngth
fs = 2*W; % sampling frequency
phi = pi/4; % azimuth
theta = pi/2; % elevation
m = [(-M/2+1/2):(M/2-1/2)]';
x_pos = m*lambda/2; % sensor positions
u_s = sin(theta); % normal vector for determining delays


%% Signal specs that do not need to be redifined in loop
n_sinusoids = 500; % number of sinusoids in signal
N_filter = K + N + 1; % filter length
n = [0:N-1]/fs; % temporal sample vectors
tau = x_pos*u_s/c; % relative delays to phase center
t_array = n-tau; % delays across array and time
t = t_array(:); %
T = max(t)-min(t);% temporal aperture of the array
L = 8;
dim = ceil(2*W*T)+L;% subspace dimension

%% PSWF interpolation
t_filt = [0:K-1]/fs - tau;
t_filt = t_filt(:);
K_interp = 4*MN;
V_u = dpss(K_interp+1,W*T,dim); % firt K uniformly sampled PSWF's
t_u = [0:(K_interp)]'; % uniform sampling points over the interval
t_nu = [K_interp*(t_filt-min(t_filt(:)))/T]; % non-uniform sample points, scaled to match interval
V_nu = interp1(t_u,V_u,t_nu,'pchip',0); % non-uniform

t_filt_nyq = [0:K-1]'/fs; % nyquist samples for testing
V_nyq = interp1(t_u,V_u,K_interp*(t_filt_nyq-min(t_filt))/T,'pchip',0); % interpolate to nyquist
v_nyq = interp1(t_u,V_u,K_interp*(t_filt_nyq(K/2)-min(t_filt))/T,'pchip',0);

%% Filter design
delta = 1e-5;
Phi = ((V_nu'*V_nu + delta*eye(dim))^-1)*(V_nu');
h = v_nyq*Phi;
H = reshape(h,M,K);
H_buff = circshift([H,zeros(M,N+1)],-K/2+1,2);
Hf = fft(H_buff,[],2);
t_nyq = [0:N-1]/fs;


%% simulation specifications
trials = 5;
SNRs = -40:6:100;
SNR_recov = zeros(trials,length(SNRs));
SNR_h_recov = zeros(trials,length(SNRs));
comp_idxs = (K+1):(N-K-1);


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
        
        y = s + noise;
        y_array = reshape(y,M,N);
        Yf = fft(y_array,N_filter,2);
        y_bf = sum(ifft(Yf.*conj(Hf),[],2),1);
        y_comp = y_bf(comp_idxs);

        
        s_nyq = sig_gen(a,f,n_sinusoids,t_nyq);
        s_comp = s_nyq(comp_idxs);
        SNR_rec = (norm(s_comp)/norm(y_comp-s_comp));
        SNR_recov(ii,jj) = SNR_rec;
        fprintf('Recovery SNR: %.2f, AG: %.2f, Lower bound AG: %.2f, Ideal AG: %.2f \n',db(SNR_rec),db(SNR_rec)-SNR_dB,db(MN/dim)/2,db(M)/2)
    end
end

figure(1)
plot(SNRs,SNRs + db(M)/2,'LineWidth',1)
hold on
grid on
plot(SNRs,db(mean(SNR_recov)),'--','LineWidth',1)
plot(SNRs,db(mean(SNR_h_recov)),'--','LineWidth',1)
xlabel('Nominal SNR (dB)','Interpreter','latex')
ylabel('Beamformed SNR (dB)','Interpreter','latex')
legend({'Ideal', 'Slepian Embedding'},'Location','northwest','Interpreter','latex')



% 
% tit = ['data/slepian_ula_ag_v_snr_L_' num2str(L) '_N_' num2str(N) '.mat'];
% save(tit,'SNR_recov','L','SNRs','trials')


%% supporting functions
function [X] = sig_gen(a,f,n_sinusoids,t)
X = zeros(size(t));
for ii = 1:n_sinusoids % For each antenna
    sig = (a(ii)*exp(1i*2*pi*f(ii)*(t)));
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