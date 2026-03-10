clear
clc
close all

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

M = 2^7; % array size
N_list = 2.^(0:9);[2,4,8,16,32,64,128,256, 512];%[2^1:8:2^6]; % no. of samples
beta_list = [1, 1.4,1.5,1.6, 2]; % oversampling factor for L

fc = 20e9;  % center frequency
c = physconst('LightSpeed');
W = 5e9; % bandwidth
fs = 2.0*W; % sampling frequency
Ts = 1/fs;
lambda = c/(fc+W); % wavelngth
theta = pi/3;%Theta(idx);
m = [(-M/2+1/2):(M/2-1/2)]';
x_pos = m*lambda/2; % sensor positions
u_s = sin(theta); % normal vector for determining delays
tau = x_pos*u_s/c;

delta = 1e-6; % regularization constant

% predefined saved values
bias = zeros(length(N_list),length(beta_list));
gamma_lower = zeros(length(N_list),1);
variance_exp = zeros(length(N_list)); % true array gain
variance_cal = zeros(length(N_list),length(beta_list));
variance_fac = zeros(length(N_list),1);
variance_gamma_lower = zeros(length(N_list),1);

for ii=1:length(N_list)
    
    N = N_list(ii);
    MN = M*N;
    n = [0:N-1]*Ts;
    t_array = n-tau;
    t = t_array(:);
    T = max(t)-min(t);
    t = t-min(t);
    gamma_lower(ii) = ((N)*Ts + (M-1)*(tau(2)-tau(1)))/((N)*Ts);
    gamma = gamma_lower(ii);%max(max(gamma_lower(ii),beta_list));
    variance_exp(ii) = (2*W*T/MN);
    
    for jj=1:length(beta_list)
        beta = beta_list(jj)*gamma;
        if mod(ceil(beta*N),2)==0
            L = ceil(beta*N)/2;
        else
            L = (ceil(beta*N)+1)/2;
        end
        fprintf('Snapshots: %d, L: %d\n',N, 2*L);
        
        ell = [-L:(L-1)];
        f_ell = ell/(2*L*Ts);
        F = exp(1i*2*pi*(t)*f_ell);
        
        A = F'*F;
        A = (A + A')/2;
        
        A_inv = inv(A + delta*eye(2*L));
        A_inv = (A_inv + A_inv')/2;
        
        Q = (A_inv*A)*A_inv;
        Q = (Q+Q')/2;
        
        T_op = (N-1)*Ts - min(t);
        % variance calc
        B = (T_op/2)*sinc((f_ell'-f_ell)*T_op);
        E = exp(1i*pi*(f_ell'-f_ell)*(T_op));
        R = Q.*(B.*E);
        
        variance_cal(ii,jj) = (2*real(sum(R(:))))/T_op;
        
        
    end
end

Legend = cell(length(beta_list)+1,1);
for jj=1:length(beta_list)
    if (jj==1)
        Legend{jj} = sprintf('L = $\\lceil \\gamma_l N \\rceil$');
    elseif (jj==length(beta_list))
        Legend{jj} = sprintf('L = $\\lceil %d \\gamma_l N \\rceil$',beta_list(jj));
    else
        Legend{jj} = sprintf('L = $\\lceil %.1f \\gamma_l N \\rceil$',beta_list(jj));
    end
end
Legend{end} = sprintf('Ideal array gain $(\\frac{\\sigma^2}{M}) $');

plot(N_list,variance_cal,'LineWidth',1)
hold on
grid on
plot(N_list, ones(size(N_list))*(1/M))
xlim([min(N_list),max(N_list)])
ylim([0,0.1])
% set(gca, 'XScale', 'log'); % Set x-axis to logarithmic scale
xticks([1 64 128 256 512 1024]); % Set specific log2 tick marks
xticklabels({'1' '2^6' '2^7' '2^8' '2^9' '2^{10}'}); % Label as powers of 2
xlabel('Snapshots (N)','Interpreter','latex','FontSize',12)
legend(Legend,'Location','northeast','Interpreter','latex','FontSize',14)
exportgraphics(gcf, 'variance_vs_N.pdf', 'ContentType', 'vector');