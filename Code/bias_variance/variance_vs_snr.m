clear
clc
close all

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

B = 2^6; % no. of beams
M = 2^7; % array size
N = 2^6; % no. of samples
L = ceil(N); % no. of samples in frequency is 2*L
MN = M*N;

norm_fac = 2*L*M; % normalizing factor for topelitz matrix

% Theta = zeros(2*B,1);
% idx = B;%randi(B,1);
% for i=1:2*B
%     Theta(i) = asin((i-1-B)/B);
% end

fc = 20e9;  % center frequency
c = physconst('LightSpeed');
lambda = c/(fc+5e9); % wavelngth
W = 5e9; % bandwidth
fs = 2.5*W; % sampling frequency
T = 1/fs;
theta = pi/3;%Theta(idx);
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
freqs = 1/(2*T*L)*linspace(-L,L-1,2*L);
F_basis = exp(1i*2*pi*t*freqs).*demod_phase;
A_eplitz = F_basis'*F_basis;

toeplitz_inv = inv(A_eplitz + (delta)*eye(2*L));
Q = toeplitz_inv*A_eplitz*toeplitz_inv;
% %% EVALUATE VARIANCE TERM
% sum_var = 0;
% for ll=1:2*L
%     for kk=1:2*L
%         f_l = (1/(2*T*L))*(ll-1-L);
%         f_k = (1/(2*T*L))*(kk-1-L);
%         if kk ~= ll
%             sum_var = sum_var + Q(ll,kk)*exp(1i*pi*(f_k - f_l)*T_aperture)*sin(pi*(f_k - f_l)*T_aperture)/(2*pi*(f_k-f_l));
%         else
%             sum_var = sum_var + Q(ll,kk)*exp(1i*pi*(f_k - f_l)*T_aperture)*T_aperture/2;
%         end
%     end
% end
exp_f = exp(1i*2*pi*T_aperture/2*(freqs' - freqs));
sincs = (T_aperture/2)*sinc(T_aperture*(freqs' - freqs));
variance_factor2 = 2*sum(Q.*exp_f.*sincs,'all');

variance_factor = 0;
for kk=1:2*L
    for ll=1:2*L
        k = kk-1-L;
        l = ll-1-L;
        f_l = (1/(2*T*L))*(l);
        f_k = (1/(2*T*L))*(k);
        w_tilde = 2*pi*(f_k - f_l);
        if k~=l 
            variance_factor = variance_factor + Q(kk,ll)*(exp(1i*2*pi*(f_k - f_l)*T_aperture/2))*(sin(2*pi*(f_k-f_l)*T_aperture/2)/(2*pi*(f_k-f_l)));
            % variance_factor = variance_factor + Q(kk,ll) * exp(1i*w_tilde*min(t)) * (exp(1i*w_tilde*T_aperture) - 1) / (1i*w_tilde);
        else
            variance_factor = variance_factor + Q(kk,ll)*T_aperture/2;
            % variance_factor = variance_factor + Q(kk,ll)*T_aperture;
        end
    end
end
variance_factor = 2*variance_factor;

%% simulation specifications
trials = 5;
SNRs = -30:6:30;
variance_exp = zeros(trials, length(SNRs));
variance_cal = zeros(trials, length(SNRs));
variance_obs = zeros(trials, length(SNRs));
variance_stoc = zeros(trials, length(SNRs));

samples = 64;
t_eq = linspace(min(t),max(t),samples)';
noise_P = exp(1i*2*pi*t_eq*freqs)*toeplitz_inv*F_basis';

for ii=1:trials
    for jj=1:length(SNRs)
        SNR_dB = SNRs(jj); % signal to noise ratio
        a = (randn(n_sinusoids,1) + 1i*randn(n_sinusoids,1))/sqrt(2); % sinusoid amplitudes
        f = W*(2*rand(n_sinusoids,1)-1); % sinusoid frequencies
        s = sig_gen(a,f,n_sinusoids,t);
        sigma = (norm(s)/sqrt(MN))*10^(-SNR_dB/20);
        noise = sigma*(randn(MN,1) + 1i*randn(MN,1))/sqrt(2);
        SNR_check = db(norm(s)/norm(noise));

        y = s.*demod_phase + noise;
        fprintf('Trial: %d, Set SNR: %.2f, Measured SNR: %.2f\n',ii,SNR_dB,SNR_check)

        % sum_var = 0;
        % for ll=1:2*L
        %     for kk=1:2*L
        %         f_l = (1/(2*T*L))*(ll-1-L);
        %         f_k = (1/(2*T*L))*(kk-1-L);
        %         if kk ~= ll
        %             sum_var = sum_var + Q(ll,kk)*exp(1i*pi*(f_k - f_l)*T_aperture)*sin(pi*(f_k - f_l)*T_aperture)/(2*pi*(f_k-f_l));
        %         else
        %             sum_var = sum_var + Q(ll,kk)*exp(1i*pi*(f_k - f_l)*T_aperture)*T_aperture;
        %         end
        %     end
        % end

        noise_variance = (norm(noise)/norm(s)).^2;
        variance_exp(ii,jj) = noise_variance/(M);

        %% measure observed variance
        s_true = sig_gen(a,f,n_sinusoids,t_eq);
        energy = integral(@(t)(abs(sig_gen(a,f,n_sinusoids,t)).^2),0,(N-1)*T,'RelTol', 1e-12, 'AbsTol', 1e-14);
        s_hat = exp(1i*2*pi*t_eq*freqs)*toeplitz_inv*(exp(1i*2*pi*t*freqs).*demod_phase)'*y;
        variance_obs(ii,jj) = ((norm(s_true - s_hat)/norm(s_true)).^2);
        norm_noise = sigma*sqrt(trace(noise_P*noise_P'));
        variance_stoc(ii,jj) = (norm_noise/norm(s)).^2;

        variance_cal(ii,jj) = (sigma^2)*variance_factor*(2*W*T_aperture)*samples*(1/(norm(s).^2));
    end
end

figure(1)
plot(SNRs,db(mean(variance_exp))/2,'LineWidth',1)
hold on
grid on
plot(SNRs,db(mean(real(variance_cal)))/2,'-^','LineWidth',1)
hold on
grid on
plot(SNRs, db(mean(variance_obs))/2, '-*','LineWidth',1)
hold on
grid on
plot(SNRs, db(mean(variance_stoc))/2,'-*','LineWidth',1)
xlabel('Nominal SNR (dB)','Interpreter','latex','FontSize',12)
ylabel('Variance ','Interpreter','latex','FontSize',12)
legend({'Ideal', 'Calculated', 'Observed', 'Stoc'},'Location','northeast','Interpreter','latex','FontSize',12)
% exportgraphics(gcf, 'variance_vs_snr.pdf', 'ContentType', 'vector');

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


