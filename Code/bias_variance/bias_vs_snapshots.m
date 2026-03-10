clear
clc
close all

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

M = 2^7; % array size
N_list = [1 2 4 8 16 32 64 72 80 92 100]; % no. of samples
beta_list = [1,1.4,1.5,1.6,2]; % oversampling factor for L

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


SNR_dB = 10;
bias = zeros(length(N_list),length(beta_list));
bias_ext = zeros(length(N_list),length(beta_list));
ext = 1.01;
bias_lower_bound = zeros(length(N_list));
trials = 10;

gamma_lower = zeros(length(N_list),1);

for ii=1:length(N_list)
    %% signal specs
    N = N_list(ii);
    MN = M*N;
    n_sinusoids = 100; % number of sinusoids in signal
    n = [0:N-1]/fs; % temporal sample vectors
    tau = x_pos*u_s/c; % relative delays to phase center
    t_array = n-tau; % delays across array and time
    t = t_array(:); %
    T_aperture = max(t)-min(t);% temporal aperture of the array
    gamma_lower(ii) = ((N-1)*T + (M-1)*sin(theta)/(2*(fc+Ws)))/(N*T);
    gamma = gamma_lower(ii);
    gamma_ext = ext*gamma_lower(ii);

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

    for jj=1:length(beta_list)
        sigma = 10^(-SNR_dB/20);
        beta = beta_list(jj);
        beta = gamma*beta;
        beta_ext = gamma_ext*beta;
        if mod(ceil(beta*N),2)==0
            L = ceil(beta*N)/2;
            L_ext = ceil(beta_ext*N)/2;
        else
            L = (ceil(beta*N)+1)/2;
            L_ext = (ceil(beta_ext*N)+1)/2;
        end
        fprintf('Snapshots: %d, L: %d\n',N, 2*ceil(beta/2*N));
    
        % demod_phase = repmat(exp(-1i*2*pi*fc*tau),N,1);
        % mod_phase = exp(1i*2*pi*fc*tau);
        freqs = 1/(2*T*L)*linspace(-L,L-1,2*L);
        freqs_ext = ext/(2*T*L_ext)*linspace(-L_ext,L_ext-1,2*L_ext);
        F_basis = exp(1i*2*pi*t*freqs);
        A_eplitz = F_basis'*F_basis;

        F_basis_ext = exp(1i*2*pi*t*freqs_ext);
        A_eplitz_ext = F_basis_ext'*F_basis_ext;
        
        %% construct the toeplitz matrix
        delta = 1e-5;
        toeplitz_inv = inv(A_eplitz + delta*eye(2*L));
        toeplitz_inv_ext = inv(A_eplitz_ext + delta*eye(2*L_ext));
        Q = toeplitz_inv*A_eplitz*toeplitz_inv;
        % Q = Q./trace(Q);
        P = toeplitz_inv*F_basis';
        P_ext = toeplitz_inv_ext*F_basis_ext';
        
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
        for ll=1:2*L_ext
            f_l_ext = (ext/(2*T*L_ext))*(ll-1-L_ext);
            F_eq_ext(:,ll) = exp(1i*2*pi*f_l_ext*t_eq);
        end
        for ll=1:2*L
            f_l = (1/(2*T*L))*(ll-1-L);
            F_eq(:,ll) = exp(1i*2*pi*f_l*t_eq);
        end
        F_prod = F_eq*P*V_nu;
        F_prod_ext = F_eq_ext*P_ext*V_nu;
        bias(ii,jj) = trace(((F_prod - V_u)'*(F_prod - V_u))*diag(Lambda/(2*W*T_aperture)));
        bias_ext(ii,jj) = trace(((F_prod_ext - V_u)'*(F_prod_ext - V_u))*diag(Lambda/(2*W*T_aperture)));
    end
end

figure(1)
Legend = cell(length(beta_list),1);
for jj=1:length(beta_list)
    semilogy(N_list,bias(:,jj),'LineWidth',1);
    if (jj==1)
        Legend{jj} = sprintf('L = $\\lceil\\gamma_l N \\rceil$');
    elseif (jj==length(beta_list))
        Legend{jj} = sprintf('L = $\\lceil %d \\gamma_l N \\rceil$',beta_list(jj));
    else
        Legend{jj} = sprintf('L = $\\lceil %.1f \\gamma_l N \\rceil$',beta_list(jj));
    end
    hold on
    grid on
end
xlim([min(N_list),max(N_list)])
set(gca, 'XScale', 'log'); % Set x-axis to logarithmic scale
xticks([1 2 4 8 16 32 64]); % Set specific log2 tick marks
xticklabels({'1' '2^1' '2^2' '2^3' '2^4'  '2^5' '2^6' }); % Label as powers of 2
xlabel('Snapshots (N)','Interpreter','latex','FontSize',12)
%ylabel('Beamformed SNR (dB)','Interpreter','latex')
legend(Legend,'Location','west','Interpreter','latex','FontSize',14)
lll = legend;
lll.Position = [0.25 0.58 0.15 0.2];  % Moves the legend to bottom-right-ish
% exportgraphics(gcf, 'bias_vs_N.pdf', 'ContentType', 'vector');

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