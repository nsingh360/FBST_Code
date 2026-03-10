clear
clc
close all

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

B = 2^7; % no. of beams
M = 2^7; % array size
N = 2^5; % no. of samples per batch
K = 2^4; % no. of batches
L = ceil(N); % no. of samples in frequency is 2*L
eta = 0.4; % relative transition width
MN = M*N;

norm_fac = 2*L*M; % normalizing factor for topelitz matrix

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
n_sinusoids = 100; % number of sinusoids in signal
n = [0:N-1]/fs; % temporal sample vectors
tau = x_pos*u_s/c; % relative delays to phase center
t_array = n-tau; % delays across array and time
t = t_array(:); %
T_aperture = max(t)-min(t);% temporal aperture of the array

%% map array aperture to K batches
T_K = zeros(MN,K);
fs_K = T_aperture/T;
W_K = fs_K/2.5;
fc_K = fc*T_aperture;
tau_K = m*u_s/(2*fc_K);
demod_phase = repmat(exp(-1i*2*pi*fc_K*tau_K),N,1);
for k=1:K
    T_K(:,k) = unifrnd(k-1-eta,k-eta,MN,1);%t/T_aperture + ((k-1-eta)*T_aperture - min(t))/T_aperture;% %
end
f_l = zeros(2*L,1);
for ll=1:2*L
    f_l(ll) = (fs_K/(2*L))*(ll-1-L);
end
% T_K = T_K(:);

% %% construct the toeplitz matrix
delta = 1e-5;
% A_eplitz = zeros(2*L,2*L);
% for ll=1:2*L
%     for kk=1:2*L
%         f_l = (1/(2*T*L))*(ll-1-L);
%         f_k = (1/(2*T*L))*(kk-1-L);
%         A_eplitz(ll,kk) = sum(exp(1i*2*pi*(f_k - f_l)*T_K));
%     end
% end
% toeplitz_inv = inv(A_eplitz + (delta)*eye(2*L));

%% simulation specifications
trials = 2;
SNRs = 30;
SNR_fbst_entire = zeros(trials,length(SNRs));
% SNR_fbst_per_batch = zeros(trials, length(SNRs));
SNR_fbst_per_batch = zeros(trials, length(SNRs),K);

for ii = 1:trials
    for jj = 1:length(SNRs)

        SNR_dB = SNRs(jj); % signal to noise ratio
        %% create K batches of observed sample
        Y = zeros(MN,K);
        S = zeros(MN,K);
        a = (randn(n_sinusoids,1) + 1i*randn(n_sinusoids,1))/sqrt(2); % sinusoid amplitudes
        f = W_K*(2*rand(n_sinusoids,1)-1); % sinusoid frequencies
        for k=1:K
            s = sig_gen(a,f,n_sinusoids,T_K(:,k));
            sigma = (norm(s)/sqrt(MN))*10^(-SNR_dB/20);
            noise = sigma*(randn(MN,1)+1i*randn(MN,1))/sqrt(2);
    
            s_demod = s.*demod_phase + noise;
    
            y_k = s_demod;
            Y(:,k) = y_k;
            S(:,k) = s;
        end

        %% K-batch reconstruction
        A_mat = zeros(MN,2*L,K);
        B_mat = zeros(MN,2*L,K-1);
        D_mat = zeros(2*L,2*L,K-1);
        DK_ = zeros(2*L,2*L);
        E_mat = zeros(2*L,2*L,K-1);
        Q_mat = zeros(2*L,2*L,K-1);
        U_mat = zeros(2*L,2*L,K-1);
        QK_ = zeros(2*L,2*L);
        V_vec = zeros(2*L,K-1);
        W_vec = zeros(2*L,K-1);
        wK_ = zeros(2*L,1);
        vK_ = zeros(2*L,1);
        for k=1:K
            tk = T_K(:,k);
            Ak = demod_phase.*g(tk-(k-1),eta).*exp(2*1i*pi*(tk-(k-1))*f_l');
            A_mat(:,:,k) = g(tk-(k-1),eta).*exp(2*1i*pi*(tk-(k-1))*f_l');
            if k<=K-1
                tk_ = T_K(:,k+1);
                Bk = demod_phase.*g(tk_-(k-1),eta).*exp(2*1i*pi*(tk_-(k-1))*f_l');
                B_mat(:,:,k) = g(tk_-(k-1),eta).*exp(2*1i*pi*(tk_-(k-1))*f_l');
    
                Ak_ = demod_phase.*g(tk_-k,eta).*exp(2*1i*pi*(tk_-k)*f_l');
    
                D_mat(:,:,k) = Ak'*Ak + Bk'*Bk + delta*eye(2*L);
                E_mat(:,:,k) = Ak_'*Bk;
    
                W_vec(:,k) = Ak'*Y(:,k) + Bk'*Y(:,k+1);
            else
                DK_ = Ak'*Ak + delta*eye(2*L);
                wK_ = Ak'*Y(:,k);
            end

            if k==1
                Q_mat(:,:,k) = D_mat(:,:,k);
            elseif k>1 && k<=K-1
                U_mat(:,:,k-1) = inv(Q_mat(:,:,k-1))*E_mat(:,:,k-1)';
                Q_mat(:,:,k) = D_mat(:,:,k) - E_mat(:,:,k-1)*U_mat(:,:,k-1);
            else
                U_mat(:,:,k-1) = inv(Q_mat(:,:,k-1))*E_mat(:,:,k-1)';
                QK_ = DK_ - E_mat(:,:,k-1)*U_mat(:,:,k-1);
            end

            if k==1
                V_vec(:,k) = inv(Q_mat(:,:,k))*W_vec(:,k);
            elseif k>1 && k<=K-1
                V_vec(:,k) = inv(Q_mat(:,:,k))*(W_vec(:,k) - E_mat(:,:,k-1)*V_vec(:,k-1));
            else 
                vK_ = inv(QK_)*(wK_ - E_mat(:,:,k-1)*V_vec(:,k-1));
            end
        end

        alphaK = zeros(2*L,K);
        alphaK(:,K) = vK_;
        for k=K-1:-1:1
            alphaK(:,k) = V_vec(:,k) - U_mat(:,:,k)*alphaK(:,k+1);
        end

        %% plot reconstruction
        for k=1:K
            t_nyq = linspace(k-1-eta,k-eta,N)';
            s_nyq = sig_gen(a,f,n_sinusoids,t_nyq);
            % tk = t_nyq.*(k-1-eta<=t_nyq).*(t_nyq < k-eta);
            % s_nyq_k = s_nyq(tk~=0);
            % tk = tk(tk~=0);
            Ak = g(t_nyq-(k-1),eta).*exp(2*1i*pi*(t_nyq-(k-1))*f_l');
            Bk = g(t_nyq-(k-2),eta).*exp(2*1i*pi*(t_nyq-(k-2))*f_l');
            if k==1
                s_recon = Ak*alphaK(:,k);
            else
                s_recon = Ak*alphaK(:,k) + Bk*alphaK(:,k-1);
            end
            SNR_fbst_per_batch(ii,jj,k) = norm(s_nyq)/norm(s_nyq - s_recon);
        end

        %% streaming least squares
        % alpha = zeros(2*L,K);
        % U_mat = zeros(2*L,2*L,K);
        % V_mat = zeros(2*L,K);
        % 
        % Ak = demod_phase.*g(T_K(:,1),eta).*exp(2*1i*pi*(T_K(:,1))*f_l.');
        % Bk = zeros(MN,2*L); %% initialize 0 batch matrices
        % Ek = zeros(2*L,2*L);
        % Uk = zeros(2*L,2*L);
        % Dk_ = Ak'*Ak + delta*eye(2*L);
        % Qk_ = Dk_;
        % wk_ = Ak'*Y(:,1);
        % vk = inv(Qk_)*wk_;
        % U_mat(:,:,1) = Uk;
        % V_mat(:,1) = vk;
        % for iter=1:K-1
        %     k = iter;
        %     tk = T_K(:,k+1);
        %     Ak_ = demod_phase.*g(tk-k,eta).*exp(2*1i*pi*(tk-k)*f_l.');
        %     Bk_ = demod_phase.*g(tk-(k-1),eta).*exp(2*1i*pi*(tk-(k-1))*f_l.');
        % 
        %     Dk = Dk_ + Bk_'*Bk_;
        %     Qk = Dk - Ek*Uk;
        % 
        %     Ek_ = Ak_'*Bk_;
        %     Uk = inv(Qk)*Ek_';
        %     Dk_ = Ak_'*Ak_ + delta*eye(2*L);
        %     Qk_ = Dk_ - Ek_*Uk;
        % 
        %     wk = wk_ + Bk_'*Y(:,k+1);
        %     wk_ = Ak_'*Y(:,k+1);
        % 
        %     vk = inv(Qk)*(wk - Ek*vk);
        %     vk_ = inv(Qk_)*(wk_ - Ek_*vk);
        % 
        %     U_mat(:,:,k+1) = Uk;
        %     V_mat(:,k+1) = vk;
        % 
        %     alpha(:,k+1) = vk_;
        %     for kk = k:-1:1
        %         alpha(:,k) = V_mat(:,k) - U_mat(:,:,k)*alpha(:,k+1);
        %     end
        % 
        % end


        % for k=1:K
        %     tk = T_K(:,k);
        %     yk = Y(:,k);
        %     if k==1
        %         Ak = demod_phase.*g(tk-k-1,eta).*exp(2*1i*pi*(tk-k-1)*f_l');
        %         Dk_ = Ak'*Ak + delta*eye(MN);
        %         Qk_ = Dk_;
        % 
        %         wk_ = Ak'*yk;
        %     else
        %         Ak = demod_phase.*g(tk-k-1,eta).*exp(2*1i*pi*(tk-k-1)*f_l');
        %         Bk = demod_phase.*g(tk-k-2,eta).*exp(2*1i*pi*(tk-k-2)*f_l');
        % 
        %         Dk = Dk_ + Bk'*Bk;
        %         if k==2
        %             Qk = Dk;
        % 
        %         else
        %             Qk = Dk - Ek*Uk;
        %         end
        %         Ek = Bk'*Ak;
        %         Uk = inv(Qk)*Ek';
        % 
        %         Dk_ = Ak'*Ak + delta*eye(MN);
        %         Qk_ = Dk_ - Ek*Uk;
        % 
        % 
        %     end
        % end

        
    end
end


labels = cell(1,K+1);
labels{1} = 'ideal';
figure(1)
plot(SNRs,SNRs + db(M)/2,'LineWidth',1)
hold on
grid on
for k=1:K
    SNR_batch_K = SNR_fbst_per_batch(:,:,k);
    plot(SNRs,db(mean(SNR_batch_K)),'-^','LineWidth',1);
    hold on
    grid on
    labels{k+1} = ['Batch ' num2str(k)];
end
xlabel('Nominal SNR (dB)','Interpreter','latex')
ylabel('Beamformed SNR (dB)','Interpreter','latex')
legend(labels,'Location','northwest','Interpreter','latex')



%% supporting functions
function [X] = sig_gen(a,f,n_sinusoids,t)
X = zeros(size(t));
for ii = 1:n_sinusoids % For each antenna
    sig = (a(ii)*exp(1i*2*pi*(f(ii))*(t)));
    X = X + sig;
end
end

%% windowing function
function out = beta(t)
out = sin((pi/2)*(sin(t*pi/2).^2));
end

function out = g(t,eta)
out = beta((t+eta)/(2*eta)).*(-eta<=t).*(t<eta) + (eta<=t).*(t<1-eta) + ...
    beta((-t+1+eta)/(2*eta)).*(1-eta<=t).*(t<=1+eta);
end