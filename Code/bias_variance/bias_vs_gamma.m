clear
clc
close all

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

M = 2^7; % array size
N = 2^6; % no. of samples
gamma_list = linspace(0.5,2.5,20); % oversampling factor for L

fc = 20e9;  % center frequency
c = physconst('LightSpeed');
Ws = 5e9;
lambda = c/(fc+Ws); % wavelngth
W = 5e9; % bandwidth
fs = 2*W; % sampling frequency
T = 1/fs;
theta = pi/3;%Theta(idx);
m = [(-M/2+1/2):(M/2-1/2)]';
x_pos = m*lambda/2; % sensor positions
u_s = sin(theta); % normal vector for determining delays

MN = M*N;
n_sinusoids = 100; % number of sinusoids in signal
n = [0:N-1]/fs; % temporal sample vectors
tau = x_pos*u_s/c; % relative delays to phase center
t_array = n-tau; % delays across array and time
t = t_array(:); %
T_aperture = max(t)-min(t);% temporal aperture of the array

ext = 1.05;
gamma_l = T_aperture*fs/N;
gamma_ext = ext*gamma_l;

%% BIAS EVALUATION
samples = 2000;
% Sinc_matrix = 2*W*sinc(2*W*(t-t'));
% Sinc_matrix = 2*W*T_aperture*Sinc_matrix/trace(Sinc_matrix);

t_eq = linspace(min(t),max(t),samples);
t_eq2 = linspace(min(t),max(t),samples);%linspace(0,(N-1)*T,samples);
% t_eq = [1000*(t_eq-min(t_eq(:)))/T_aperture]; % non-uniform sample points, scaled to match interval

% Sinc_matrix_eq = 2*W*sinc(2*W*(t_eq-t_eq'));
% Sinc_matrix_eq = 2*W*T_aperture*Sinc_matrix_eq/trace(Sinc_matrix_eq);
% [V_eq,D_eq] = eigs(Sinc_matrix_eq,K);
K = ceil(2*W*T_aperture) + 8;
[V_u,Lambda] = dpss(samples,W*T_aperture,K);
V_nu = interp1(t_eq,V_u,t,'pchip',1); % non-uniform

SNR_dB = 10;
bias_fourier = zeros(length(gamma_list),1);
bias = zeros(length(gamma_list),1);
bias_ext = zeros(length(gamma_list),1);
trials = 10;


for jj=1:length(gamma_list)
        gamma = gamma_list(jj);
        if mod(ceil(gamma*N),2)==0
            L = ceil(gamma*N)/2;
            
        else
            L = (ceil(gamma*N)+1)/2;
        end

        fprintf('Snapshots: %d, L: %d \n',N, 2*ceil(gamma/2*N));
        L_ext = L;
    
        % demod_phase = repmat(exp(-1i*2*pi*fc*tau),N,1);
        % mod_phase = exp(1i*2*pi*fc*tau);
        freqs_fourier = 1/(T_aperture)*linspace(-L,L-1,2*L);
        freqs = 1/(2*T*L)*linspace(-L,L-1,2*L);
        freqs_ext = ext/(2*T*L_ext)*linspace(-L_ext,L_ext-1,2*L_ext);
        F_basis = exp(1i*2*pi*t*freqs);
        A_eplitz = F_basis'*F_basis;

        F_basis_ext = exp(1i*2*pi*t*freqs_ext);
        A_eplitz_ext = F_basis_ext'*F_basis_ext;

        F_basis_fourier = exp(1i*2*pi*t*freqs_fourier);
        A_eplitz_fourier = F_basis_fourier'*F_basis_fourier;
        
        %% construct the toeplitz matrix
        delta = 1e-5;
        toeplitz_inv = inv(A_eplitz + delta*eye(2*L));
        toeplitz_inv_ext = inv(A_eplitz_ext + delta*eye(2*L_ext));
        toeplitz_inv_fourier = inv(A_eplitz_fourier + delta*eye(2*L_ext));
        Q = toeplitz_inv*A_eplitz*toeplitz_inv;
        % Q = Q./trace(Q);
        P = toeplitz_inv*F_basis';
        P_ext = toeplitz_inv_ext*F_basis_ext';
        P_fourier = toeplitz_inv_fourier*F_basis_fourier';
        
        % %% construct slepian matrix
        % K = ceil(2*W*T_aperture) + 4; % number of diagonal entries to consider
        % V_u = dpss(MN+1,W*T_aperture,K); % firt K uniformly sampled PSWF's
        % t_u = [0:(MN)]'; % uniform sampling points over the interval
        % t_nu = [M*N*(t-min(t(:)))/T_aperture]; % non-uniform sample points, scaled to match interval
        % V_nu = interp1(t_u,V_u,t_nu,'pchip',1); % non-uniform
        
        % %% equispaced sampled basis matrix
        % psi = zeros(N,2*L);
        % t_nyq = [0:N-1]/fs; % nyquist samples for testing
        % for ll=1:2*L
        %     f_l = (1/(2*T*L))*(ll-1-L);
        %     psi(:,ll) = exp(1i*2*pi*f_l*t_nyq);
        % end
        
        F_eq = zeros(samples,2*L);
        F_eq_ext = zeros(samples,2*L_ext);
        F_eq_fourier = zeros(samples, 2*L);
        for ll=1:2*L_ext
            f_l_ext = (ext/(2*T*L_ext))*(ll-1-L_ext);
            F_eq_ext(:,ll) = exp(1i*2*pi*f_l_ext*t_eq);
        end
        for ll=1:2*L
            f_l = (1/(2*T*L))*(ll-1-L);
            F_eq(:,ll) = exp(1i*2*pi*f_l*t_eq);
        end
        for ll=1:2*L
            f_l = (1/(T_aperture))*(ll-1-L);
            F_eq_fourier(:,ll) = exp(1i*2*pi*f_l*t_eq);
        end
        F_prod = F_eq*P*V_nu;
        F_prod_ext = F_eq_ext*P_ext*V_nu;
        F_prod_fourier = F_eq_fourier*P_fourier*V_nu;
        bias_fourier(jj,1) = trace(((F_prod_fourier - V_u)'*(F_prod_fourier - V_u))*diag(Lambda/(2*W*T_aperture)));
        bias(jj,1) = trace(((F_prod - V_u)'*(F_prod - V_u))*diag(Lambda/(2*W*T_aperture)));
        bias_ext(jj,1) = trace(((F_prod_ext - V_u)'*(F_prod_ext - V_u))*diag(Lambda/(2*W*T_aperture)));
end


figure(1)
semilogy(gamma_list,bias_fourier,'LineWidth',1);
hold on
grid on
semilogy(gamma_list,bias_ext,'LineWidth',1)
xline(gamma_ext,'Color',[0.8500 0.3250 0.0980],'LineStyle','--','LineWidth',2)
% hold on
% xline(gamma_ext,'Color',[0.9290, 0.6940, 0.1250],'LineStyle','--')
xlabel('$\gamma$','Interpreter','latex','FontSize',12)
ylabel('Bias','Interpreter','latex','FontSize',12)
legend({'Fourier series', 'Fourier extension'},'Location','northeast','Interpreter','latex', 'FontSize',12)
exportgraphics(gcf, 'bias_vs_gamma.pdf', 'ContentType', 'vector');