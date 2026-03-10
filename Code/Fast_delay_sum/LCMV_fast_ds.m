clear
clc
close all

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

%% array spatial and temporal specificiations
M = 2^6; % array size
N = 2^7; % temporal samples
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
n_sinusoids = 100; % number of sinusoids in signal
n = [0:N-1]/fs; % temporal sample vectors
tau = x_pos*u_s/c; % relative delays to phase center
t_array = n - tau; % delays across array and time
t = t_array(:); %
T = max(t)-min(t);% temporal aperture of the array

%% Interferer specs that do not need to be redfined in loop
phi_i = -pi/4;
theta_i = pi/3;
u_i = sin(theta_i); % normal vector for determining delays
tau_i = x_pos*u_i/c;
t_i_array = n - tau_i;
t_i = t_i_array(:);
T_i = max(t_i) - min(t_i);
e_i  = repmat(exp(1i*2*pi*fc*(-tau_i+tau)),1,N);

%% MVDR filter design
K =16;
KM = M*K;
N_filter = K + N + 1; % filter length
t_mvdr = t_array(:,1:K);
t_mvdr = t_mvdr(:);
T_mvdr = max(t_mvdr)-min(t_mvdr);
t_i_mvdr = t_i_array(:,1:K); % temporal points over which interferer is seen
t_i_mvdr = t_i_mvdr(:);
T_i_mvdr = max(t_i_mvdr)-min(t_i_mvdr);
L_s = 6;
L_i = 15;

% signal covariance matrix
[Us, Cs] = prolate_approx(KM*(t_mvdr-min(t_mvdr))/T_mvdr,W,T_mvdr,L_s,1);
Cs_inv = Cs^-1;

% interferer covariance matrix
e_i_mvdr  = repmat(exp(1i*2*pi*fc*(-tau_i+tau)),K,1);
[Ui, Ci] = prolate_approx(KM*(t_i_mvdr-min(t_i_mvdr))/T_i_mvdr,W,T_i_mvdr,L_i,1);
Ci_inv = Ci^-1;
Ui = e_i_mvdr.*Ui;

U_full = [Us,Ui];
U_gram = U_full'*U_full;

C = kron(eye(K),ones(M,1)); % constraint matrix
F = eye(K,1); % response vector

% constraint design
J = ceil(2*W*(max(t_mvdr)-min(t_mvdr)))+9;
J2 = ceil(2*W*(max(t_i_mvdr) - min(t_i_mvdr))) + 9;
J_t = J+J2;
f_J = linspace(-W,W,10*J);
f_J2 = linspace(-W,W,10*J2);
Z = exp(1i*2*pi*t_mvdr.*f_J);
Z2 = exp(1i*2*pi*t_i_mvdr.*f_J2);
[U, E, V] = svd([Z,Z2],'econ');

C = U(:,1:J_t); % constraint matrix
F = (E(1:J_t,1:J_t)^-1)*V(:,1:J_t)'*[ones(length(f_J),1);0.1*ones(length(f_J2),1)]; % response vector

C_con = Z;
C_con2 = [Z,Z2];
F_con = [ones(length(f_J),1)];
F_con2 = [ones(length(f_J),1);zeros(length(f_J2),1)];

delta = 1e-5;
w_q = C_con*inv(C_con'*C_con + delta*eye(370))*F_con;

B = null(C_con');

%% simulation specifications
trials = 5;
SIRs = -30;
SNRs = -30:6:30;
SNR_recov = zeros(trials,length(SNRs),length(SIRs));
t_nyq = [0:N-1]'/fs; % nyquist samples for testing
comp_idxs = (K+1):(N-K-1);


for jj = 1:length(SNRs)
    for kk = 1:length(SIRs)
        SIR_dB = SIRs(kk); % signal to interferer ratio
        SNR_dB = SNRs(jj); % signal to noise ratio
        sigma = 10^(-SNR_dB/20);
        sigma_h = (sigma^2)/MN;
        Z = [Cs_inv, zeros(size(Cs,1),size(Ci,1));zeros(size(Ci,1),size(Cs,1)), (10^(SIR_dB/20))*Ci_inv];
        G = ((Z + (1/sigma_h)*(U_gram))^-1)/sigma_h^2;
        RC = C/sigma_h - U_full*(G*(U_full'*C));
        Rxx = inv(eye(M*K)./sigma_h-U_full*(G*(U_full')));
        h = w_q - B*inv(B'*Rxx*B)*B'*Rxx*w_q;
        % h = RC*(((C'*RC)^-1)*F);
        H = reshape(h,M,K); % filter bank form
        Hf = fft(H,N_filter,2);
        
        for ii = 1:trials
            % signal generation
            a = (randn(n_sinusoids,1) + 1i*randn(n_sinusoids,1))/sqrt(2); % sinusoid amplitudes
            f = W*(2*rand(n_sinusoids,1)-1); % sinusoid frequencies
            s_array = sig_gen(a,f,n_sinusoids,t_array);
            s_scale = 1/rms(s_array(:));
            s_array = s_array*s_scale;
            
            % interferer generation
            a_i = (randn(n_sinusoids,1) + 1i*randn(n_sinusoids,1))/sqrt(2); % random sinusoid coefficients
            f_i = W*(2*rand(n_sinusoids,1)-1); % random sinusoid freqs
            s_i_array = e_i.*sig_gen(a_i,f_i,n_sinusoids,t_i_array);
            s_i_scale = 1/rms(s_i_array(:));
            s_i_array = s_i_scale*s_i_array*10^(-SIR_dB/20);
            SIR_check = db(norm(s_array(:))/norm(s_i_array(:)));
            fprintf('Trial: %d, Set SIR: %.2f, Measured SIR: %.2f\n',ii,SIR_dB,SIR_check)
            
            
            % noise generation
            noise_array = sigma*(randn(M,N)+1i*randn(M,N))/sqrt(2);
            SNR_check = db(norm(s_array(:))/norm(noise_array(:)));
            fprintf('Set SNR: %.2f, Measured SNR: %.2f\n',SNR_dB,SNR_check)
            
            y = s_array + s_i_array + noise_array;
            yf = fft(y ,N_filter,2);
            y_filter = ifft(conj(Hf).*yf,[],2); % perform filtering
            y_bf = sum(y_filter(:,1:N),1); % combine channels
            
            s_nyq = sig_gen(a,f,n_sinusoids,n)*s_scale;
            SNR_rec = (norm(s_nyq(comp_idxs))/norm(y_bf(comp_idxs)-s_nyq(comp_idxs)));
            SNR_recov(ii,jj,kk) = SNR_rec;
            fprintf('Recovery SNR: %.2f, AG: %.2f, Lower bound AG: %.2f, Ideal AG: %.2f \n',db(SNR_rec),db(SNR_rec)-SNR_dB,db(MN/K)/2,db(M)/2)
        end
    end
end

figure(1)
plot(SNRs,SNRs + db(M)/2,'LineWidth',1)
hold on
grid on
plot(SNRs,db(mean(SNR_recov)),'--','LineWidth',1)
xlabel('Nominal SNR (dB)','Interpreter','latex')
ylabel('Beamformed SNR (dB)','Interpreter','latex')
legend({'Ideal', 'Slepian Embedding'},'Location','northwest','Interpreter','latex')




tit = ['data/lcmv_ula_ag_v_snr_sinr_K_' num2str(K) '.mat'];
save(tit,'SNR_recov','K','SNRs','trials')


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