%% Bound over the original time interval (not mapped to [-1,1])
clear
clc
close all
%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

gamma = 2.5; % oversampling for nyquist criteria
beta = 2; % oversampling for fourier extension
B = 2^6; % no. of beams
M = 2^6; % array size
N = 2^4; % no. of samples
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
% t = sort(t);
% t = unifrnd(min(t),max(t),MN,1);
T_aperture = max(t)-min(t);% temporal aperture of the array

demod_phase = repmat(exp(-1i*2*pi*fc*tau),N,1);
%% define the fourier extension basis and covariance matrix
covariance = zeros(2*L,2*L);
delta = 1e-5;
F_L = linspace(-L,L-1,2*L)*fs/(2*L);
A_basis = exp(1i*2*pi*t*F_L);
for ll=1:2*L
    f_l = (fs/(2*L))*(ll-1-L);
    for kk=1:2*L
        f_k = (fs/(2*L))*(kk-1-L);
        covariance(ll,kk) = integral(@(t)exp(2*pi*1i*(f_k-f_l)*t),min(t),max(t),'RelTol', 1e-12, 'AbsTol', 1e-14)/T_aperture;
    end
end
A_demoded = A_basis.*demod_phase;
Phi = (A_demoded'*A_demoded + delta*eye(2*L))\A_demoded';
sample_covariance = A_basis'*(A_basis)./MN;
[V,D] = eig(sample_covariance);
eig_vals = diag(D);
A_eplitz = A_basis'*A_basis;

%% equispaced sampled basis matrix
t_nyq = [0:N-1]/fs; % nyquist samples for testing
psi = exp(1i*2*pi*t_nyq.'*F_L);

% %% check exact covariance matrix
% F_hat = (F_L.' - F_L).';
% cov_check = exp(-1i*2*pi*F_hat*(M-1)/2*sin(theta)/(2*fc)).*(exp(1i*2*pi*F_hat*T_aperture)-1)./(T_aperture*1i*2*pi*F_hat);
% for i=1:2*L
%     cov_check(i,i) = 1;
% end
% 
% %% check sample_covariance matrix 
% C1 = exp(-1i*2*pi*(F_L.'-F_L).'*(M-1)/2*sin(theta)/(2*fc));
% C2 = 1./(1 - exp(-1i*2*pi*(F_L.'-F_L).'*sin(theta)/(2*fc)));
% C3 = 1./(1 - exp(-1i*2*pi*(F_L.'-F_L).'*T));
% C4 = exp(1i*2*pi*(F_L.' - F_L).'*(T_aperture)) + exp(-1i*2*pi*(F_L.' - F_L).'*(T + sin(theta)/(2*fc))) - exp(1i*2*pi*(F_L.' - F_L).'*((M-1)*sin(theta)/(2*fc) - T)) ...
%     - exp(1i*2*pi*(F_L.' - F_L).'*((N-1)*T - sin(theta)/(2*fc)));
% sample_cov_check = C1.*C2.*C3.*C4./MN;
% for i=1:2*L
%     sample_cov_check(i,i) = 1;
% end

%% check matrix rank
R = zeros(2*L,2*L);
A1 = exp(-1i*2*pi*(F_L.'-F_L).'*(M-1)/2*sin(theta)/(2*fc));
A2 = 1./(1 - exp(-1i*2*pi*(F_L.'-F_L).'*sin(theta)/(2*fc)));
A3 = 1./(1 - exp(-1i*2*pi*(F_L.'-F_L).'*T));
A4 = 1 + exp(-1i*2*pi*(F_L.' - F_L).'*(T + sin(theta)/(2*fc))) - exp(1i*2*pi*(F_L.' - F_L).'*((M-1)*sin(theta)/(2*fc) - T)) ...
    - exp(1i*2*pi*(F_L.' - F_L).'*((N-1)*T - sin(theta)/(2*fc)));
R = A1.*A2.*A3.*A4;
for i=1:2*L
    R(i,i) =  0;%MN*(T + tau(2)-tau(1))/T_aperture-1;
end
% 
Psi_0 = A_basis.';
samples = 5;
t_grid = linspace(0,(N-1)*T,samples);
Psi_grid = exp(1i*2*pi*(F_L.')*t_grid);
coeff_grid = zeros(samples, MN);
coeff_grid_tsvd = zeros(samples,MN);
err_grid = zeros(samples,2*L);
err_grid_tsvd = zeros(samples, 2*L);
t_sorted = sort(t);
err_grid_bn = zeros(samples,1);
err_norm = zeros(samples,1);
delta_max = max(diff(sort(t)));
for i=1:samples
    t_pt = t_grid(i);
    lower_idx = find(t_sorted <= t_pt, 1, 'last'); % Last element less than or equal to sample
    upper_idx = find(t_sorted >= t_pt, 1, 'first');
    delta_t = t_sorted(upper_idx) - t_sorted(lower_idx);
    % coeff_grid(i,:) = l1_ls(Psi_0, Psi_grid(:,i),10);
    coeff_grid(i,:) = inv(Psi_0'*Psi_0 + 1e-5*eye(MN))*Psi_0'*Psi_grid(:,i);
    coeff_grid_tsvd(i,:) = Psi_0\Psi_grid(:,i);
    err_grid(i,:) = Psi_grid(:,i) - Psi_0*coeff_grid(i,:).';
    err_grid_tsvd(i,:) = Psi_grid(:,i) - Psi_0*coeff_grid_tsvd(i,:).';
    err_norm(i) = norm(err_grid(i,:));
    err_grid_bn(i,1) = sum((1 - cos(2*pi*F_L*delta_max/2)).^2);
end

%% cubic lagrange interpolation using 4 points
t_grid = linspace(2*T,(N-1)*T,samples);

for i=1:length(t_grid)
    
end

err_term = norm(err_grid,'fro');
coeff_term = norm(coeff_grid,'fro');

err_term_tsvd = norm(err_grid_tsvd,'fro');
coeff_term_tsvd = norm(coeff_grid_tsvd,'fro');

trials = 5;
SNRs = -30:6:30;
noise_var_db = zeros(length(SNRs),1);
SNR_fe = zeros(trials,length(SNRs));
SNR_fe_chirp = zeros(trials, length(SNRs));

err_non_unif = zeros(trials, length(SNRs));
err_non_unif_bound = zeros(trials, length(SNRs));
err_non_unif_bound2 = zeros(trials, length(SNRs)); % bounding the mismatch bias

for ii=1:trials
    for jj=1:length(SNRs)

        SNR_dB = SNRs(jj); % signal to noise ratio
        a = (randn(n_sinusoids,1) + 1i*randn(n_sinusoids,1))/sqrt(2); % sinusoid amplitudes
        f = W*(2*rand(n_sinusoids,1)-1); % sinusoid frequencies
        s = sig_gen(a,f,n_sinusoids,t);
        sigma = (norm(s)/sqrt(MN))*10^(-SNR_dB/20);
        noise = sigma*(randn(MN,1)+1i*randn(MN,1))/sqrt(2);
        SNR_check = db(norm(s)/norm(noise));
        fprintf('Trial: %d, Set SNR: %.2f, Measured SNR: %.2f\n',ii,SNR_dB,SNR_check)

        s_demod = s.*demod_phase;

        y = s_demod+noise;
        alpha = Phi*y;

        %% calculate the convergence error between continuous FE and the signal (The bias)
        T_extended = (beta)*(1/fs)*N/T_aperture; 
        [alpha_true,error,error_bound,max_val,rho] = fe_bias_error(@sig_gen, @fe_continuous, a, f, n_sinusoids, T_extended,fs, W, L,t,1);
        
        % 
        %% calculate variance bound assuming fixed design
        sum_1 = 0;
        sum_2 = 0;
        for ll=1:2*L
            sum_1 = sum_1 + eig_vals(ll)/((eig_vals(ll)/delta + 1)^2)*(abs(V(:,ll)'*alpha_true))^2;
            sum_2 = sum_2 + (eig_vals(ll)/(eig_vals(ll) + delta))^2;
        end
        % error_bias_bound(ii,jj) = error_bound;

        % error_bound_discrete(ii,jj) = error_bound + coeff_term*sqrt(error_variance_bound(ii,jj)) + err_term*norm(alpha_true-alpha);

        %% try non-uniform bound
        s_grid = sig_gen(a,f,n_sinusoids,t_grid);
        err_non_unif(ii,jj) = norm(s_grid.' - Psi_grid.'*alpha)./samples;
        term1 =  coeff_term*sqrt(MN*beta*sigma^2/M);
        term2 = err_term*(norm(alpha)+norm(alpha_true));
        term2_bnd = err_term_tsvd*(norm(alpha)+norm(alpha_true));
        term1_bnd = coeff_term_tsvd*sqrt(MN*beta*sigma^2/M);
        err_non_unif_bound(ii,jj) = error+(term1+term2)./samples;
        err_non_unif_bound2(ii,jj) = error + (term1_bnd + term2_bnd)./samples;
        % 
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
            % offset_coeff = exp(1i*pi/L*offset_n*l_pts);
            % F_transform(mm,:) = F_transform(mm,:).*offset_coeff; % Adjust for linear mapping
        end
        X_czt = zeros(2*L,B);
        A= 1;
        B_list = reshape(linspace(0,B-1,B), [1,B]);
        for ll=1:2*L
            l_dash = ll-1;
            % A = exp(-1i*pi*(1 + 2*(W/fc)*(l_dash-L/2)/L));
            W_czt = exp(1i*pi*(1/B)*(1 + (1/(2*fc*1/fs*L))*(l_dash-L)));
            coeff = exp(-1i*pi*(1/B)*(1 + (1/(2*fc*1/fs*L))*(l_dash-L))*(M/2-1/2)*B_list);
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

%% plot total error and total error bound
figure(4)
semilogy(SNRs, mean(err_non_unif),'-o','LineWidth',1)
hold on
grid on
semilogy(SNRs, mean(err_non_unif_bound),'-^','LineWidth',1)
hold on
grid on
semilogy(SNRs, mean(err_non_unif_bound2),'-^','LineWidth',1)
xlabel('Nominal SNR (dB)','Interpreter','latex')
ylabel('Total error','Interpreter','latex')
legend({'Error Exact','Error bound', 'Error bound2'},'Location','northwest','Interpreter','latex')

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

function [X] = raise_n(F_L, tau, T, M, N, n)

Fm = (F_L.' - F_L).^n;
fctiral = factorial(n);
X = ((1i)^n)*((2*pi).^n)*Fm*((N*T).^2 + (M*tau).^n)./fctiral;

end
