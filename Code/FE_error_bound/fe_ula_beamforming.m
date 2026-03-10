clear
clc
close all

%% NOTE : SOLVE for the orignial time interval mapping to [-1,1] creates wrong covariance matrix

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

gamma = 2.5; % oversampling for nyquist criteria
beta = 2; % oversampling for fourier extension
B = 2^7; % no. of beams
M = 2^7; % array size
N = 2^6; % no. of samples
L = ceil(beta/2 * N); % no. of samples in frequency is 2*L
MN = M*N;

Theta = zeros(B,1);
idx = B;%randi(B,1);
for i=1:B
    Theta(i) = asin((i-1)/B);
end

fc = 20e9;  % center frequency
c = physconst('LightSpeed');
lambda = c/fc; % wavelngth
W = 5e9; % bandwidth
fs = gamma*W; % sampling frequency
T = 1/fs;
theta = Theta(idx);
m = [(-M/2+1/2):(M/2-1/2)]';
x_pos = m*lambda/2; % sensor positions
u_s = sin(theta); % normal vector for determining delays

%% Signal specs that do not need to be redifined in loop
n_sinusoids = 10; % number of sinusoids in signal
n = [0:N-1]/fs; % temporal sampple vectors
tau = x_pos*u_s/c; % relative delays to phase center
t_array = n-tau; % delays across array and time
t = t_array(:); %
T_aperture = max(t)-min(t);% temporal aperture of the array

%% Map the array aperture to [-1,1]
t_st = (2/T_aperture)*t - (max(t) + min(t))/T_aperture;
% t_st = unifrnd(-1,1,MN,1); % define uniform time sampling
fs_st = T_aperture/(2*T);
W_st = fs_st/gamma;
fc_st = fc*(T_aperture/2);
tau_st = m*u_s/(2*fc_st);
demod_phase = repmat(exp(-1i*2*pi*fc_st*tau_st),N,1);
mod_phase = exp(1i*2*pi*fc_st*tau_st);

%% define the fourier extension basis and covariance matrix
A_basis = zeros(MN,2*L);
covariance = zeros(2*L,2*L);
delta = 1e-5;
for ll=1:2*L
    f_l = (fs_st/(2*L))*(ll-1-L);
    A_basis(:,ll) = exp(1i*2*pi*(f_l)*t_st);
    for kk=1:2*L
        f_k = (fs_st/(2*L))*(kk-1-L);
        covariance(ll,kk) = 0.5*integral(@(t)exp(2*pi*1i*(f_k-f_l)*t),-1,1);
    end
end
A_demoded = A_basis.*demod_phase;
Phi = (A_demoded'*A_demoded + delta*eye(2*L))\A_demoded';
sample_covariance = A_basis'*(A_basis)/MN;
[V,D] = eig(sample_covariance);
eig_vals = diag(D);

%% bounded statistical leverage condition
numerator = zeros(MN,1);
for mn=1:MN
    numerator(mn,1) = norm(inv(covariance + delta*eye(2*L))*(A_basis(mn,:)'));
end
denominator = sqrt(trace(inv(covariance + delta*eye(2*L))*covariance));

delta = 1e-5;
A_eplitz = zeros(2*L,2*L);
for ll=1:2*L
    for kk=1:2*L
        f_l = (fs_st/(2*L))*(ll-1-L);
        f_k = (fs_st/(2*L))*(kk-1-L);
        A_eplitz(ll,kk) = sum(exp(1i*2*pi*(f_k - f_l)*t_st));
    end
end

%% equispaced sampled basis matrix
psi = zeros(N,2*L);
t_nyq = [-N/2:N/2-1]/fs_st; % nyquist samples for testing
for nn=1:N
    for ll=1:2*L
        f_l = (fs_st/(2*L))*(ll-1-L);
        psi(nn,ll) = exp(1i*2*pi*f_l*t_nyq(nn));
    end
end

%% create a finely sampled grid Tau
% Psi_0 = A_basis.';
% samples = 10000;
% t_grid = linspace(-1,1,samples);
% F_L = (fs_st/(2*L))*linspace(-L,L-1,2*L);
% Psi_grid = exp(1i*2*pi*(F_L.')*t_grid);
% coeff_grid = zeros(samples, MN);
% err_grid = zeros(samples,2*L);
% for i=1:samples
%     coeff_grid(i,:) = inv(Psi_0'*Psi_0 + 1e-8*eye(MN))*Psi_0'*Psi_grid(:,i);
%     err_grid(i,:) = Psi_grid(:,i) - Psi_0*coeff_grid(i,:).';
% end
% err_term = norm(err_grid,'fro');
% coeff_term = norm(coeff_grid,'fro');

trials = 2;
SNRs = -30:6:30;
noise_var_db = zeros(length(SNRs),1);
SNR_fe = zeros(trials,length(SNRs));
SNR_fe_chirp = zeros(trials, length(SNRs));
error_variance = zeros(trials, length(SNRs));
error_variance_bound = zeros(trials, length(SNRs));
error_variance_bound2 = zeros(trials, length(SNRs));

error_bias = zeros(trials, length(SNRs));
error_bias_bound = zeros(trials, length(SNRs));
error_exact = zeros(trials, length(SNRs));

error_bound_discrete = zeros(trials, length(SNRs));

for ii=1:trials
    for jj=1:length(SNRs)

        SNR_dB = SNRs(jj); % signal to noise ratio
        a = (randn(n_sinusoids,1) + 1i*randn(n_sinusoids,1))/sqrt(2); % sinusoid amplitudes
        f = W_st*(2*rand(n_sinusoids,1)-1); % sinusoid frequencies
        s = sig_gen(a,f,n_sinusoids,t_st);
        sigma = (norm(s)/sqrt(MN))*10^(-SNR_dB/20);
        noise = sigma*(randn(MN,1)+1i*randn(MN,1))/sqrt(2);
        SNR_check = db(norm(s)/norm(noise));
        fprintf('Trial: %d, Set SNR: %.2f, Measured SNR: %.2f\n',ii,SNR_dB,SNR_check)

        s_demod = s.*demod_phase;

        y = s_demod + noise;
        alpha = Phi*y;

        %% calculate the convergence error between continuous FE and the signal (The bias)
        T_extended = (beta/2)*(1/fs_st)*N; 
        [alpha_true,error,error_bound] = fe_bias_error(@sig_gen, @fe_continuous, a, f, n_sinusoids, T_extended,fs_st, W_st, L,t_st,0);
        error_variance(ii,jj) = abs((alpha_true - alpha)'*sample_covariance*(alpha_true-alpha));
        error_variance_bound(ii,jj) = beta*sigma^2/M;

        %% calculate variance bound assuming fixed design
        sum_1 = 0;
        sum_2 = 0;
        for ll=1:2*L
            sum_1 = sum_1 + eig_vals(ll)/((eig_vals(ll)/delta + 1)^2)*(abs(V(:,ll)'*alpha_true))^2;
            sum_2 = sum_2 + (eig_vals(ll)/(eig_vals(ll) + delta))^2;
        end
        error_variance_bound2(ii,jj) = 2*(sum_1 + sigma^2/MN*sum_2);

        error_exact(ii,jj) = sqrt(integral(@(t)(abs(sig_gen(a,f,n_sinusoids,t) - fe_continuous(alpha,fs_st,L,t))).^2,-1,1));
        total_En = integral(@(t_int)(abs(sig_gen(a,f,n_sinusoids,t_int))).^2,-1,1);

        error_bias(ii,jj) = error;
        error_bias_bound(ii,jj) = error_bound;

        error_bound_discrete(ii,jj) = error_bound + coeff_term*sqrt(error_variance_bound(ii,jj)) + err_term*norm(alpha_true-alpha);

        y_nyq = psi*alpha;
        s_nyq = sig_gen(a,f,n_sinusoids,t_nyq');
        SNR_fe(ii,jj) = norm(s_nyq)/norm(s_nyq - y_nyq);

        %% chirp-z 
        y_z = reshape(y, [M,N]);
        F_transform = zeros(M,2*L);
        A = -1;%exp(-1i*2*pi*T*W);
        W_czt = exp(-1i*pi/L);
        l_pts = reshape(linspace(-L,L-1,2*L), [1,2*L]);
        offset_n = (max(t)+min(t))/(2*T);
        for mm=1:M
            F_transform(mm,:) = czt(y_z(mm,:),2*L,W_czt,A);
            offset_coeff = exp(1i*pi/L*offset_n*l_pts);
            F_transform(mm,:) = F_transform(mm,:).*offset_coeff; % Adjust for linear mapping
        end
        X_czt = zeros(2*L,B);
        A= 1;
        B_list = reshape(linspace(0,B-1,B), [1,B]);
        for ll=1:2*L
            l_dash = ll-1;
            % A = exp(-1i*pi*(1 + 2*(W/fc)*(l_dash-L/2)/L));
            W_czt = exp(1i*pi*(1/B)*(1 + (1/(2*fc_st*1/fs_st*L))*(l_dash-L)));
            coeff = exp(-1i*pi*(1/B)*(1 + (1/(2*fc_st*1/fs_st*L))*(l_dash-L))*(M/2-1/2)*B_list);
            X_czt(ll,:) = czt(F_transform(:,ll),B,W_czt,A);
            X_czt(ll,:) = X_czt(ll,:).*coeff; % Adjust for array center
        end
        w_l = X_czt(:,idx);

        alpha_z = inv(A_eplitz + delta*eye(2*L))*w_l;
        y_nyq_z = psi*alpha_z;
        SNR_fe_chirp(ii,jj) = norm(s_nyq)/norm(s_nyq - y_nyq_z);

    end
end

%% plot beamformed SNR
figure(1)
plot(SNRs,SNRs + db(M)/2,'LineWidth',1)
hold on
grid on
plot(SNRs,db(mean(SNR_fe)),'-^','LineWidth',1)
hold on
grid on
plot(SNRs, db(mean(SNR_fe_chirp)),'-^','LineWidth',1)
xlabel('Nominal SNR (dB)','Interpreter','latex')
ylabel('Beamformed SNR (dB)','Interpreter','latex')
legend({'Ideal','FE','FE Chirp-z'},'Location','northwest','Interpreter','latex')

%% plot bias and bias bound
figure(2)
semilogy(flip(SNRs), flip(mean(error_bias.^2)),'-o','LineWidth',1)
hold on
grid on
semilogy(flip(SNRs), flip(mean(error_bias_bound.^2)),'-^','LineWidth',1)
xlabel('Noise variance (dB)','Interpreter','latex')
ylabel('Bias','Interpreter','latex')
legend({'bias $\vert\vert f - f_{N}\vert\vert$','bias bound'},'Location','northwest','Interpreter','latex')

%% plot variance and variance bound
figure(3)
semilogy(flip(SNRs), flip(mean(error_variance)),'-o','LineWidth',1)
hold on
grid on
semilogy(flip(SNRs), flip(mean(error_variance_bound)),'-^','LineWidth',1)
hold on
grid on
xlabel('Noise variance (dB)','Interpreter','latex')
ylabel('Variance','Interpreter','latex')
legend({'variance $\vert\vert \hat{f}_{N} - f_{N}\vert\vert$','variance bound'},'Location','northwest','Interpreter','latex')

%% plot total error and total error bound
figure(4)
semilogy(SNRs, mean(error_exact),'-o','LineWidth',1)
hold on
grid on
semilogy(SNRs, mean(error_bias_bound + sqrt(2*error_variance_bound)),'-^','LineWidth',1)
hold on
grid on
semilogy(SNRs, mean(error_bound_discrete),'-^','LineWidth',1)
xlabel('Nominal SNR (dB)','Interpreter','latex')
ylabel('Total error','Interpreter','latex')
legend({'error exact','error bound (bias variance)','error bound (discrete)'},'Location','northwest','Interpreter','latex')



%% supporting functions
function [X] = fe_continuous(alpha,fs,L,t)
X=0;
for l=1:2*L
    f_l = (fs/(2*L))*(l-1-L);
    X = X+alpha(l)*exp(1i*2*pi*f_l*t);
end
end

function [X] = sig_gen(a,f,n_sinusoids,t)
X = zeros(size(t));
for ii = 1:n_sinusoids % For each antenna
    sig = (a(ii)*exp(1i*2*pi*(f(ii))*(t)));
    X = X + sig;
end
end

function [d] = effective_dimension(eigen_val, delta, p)
d = 0;
for i=1:length(eigen_val)
    d = d + (eigen_val(i,1)/(eigen_val(i,1) + delta))^p;
end
end