clear
clc
close all

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

%% array spatial and temporal specificiations
M = 2^3; % array size
N = 2^5; % temporal samples
MN = M*N;
fc = 6e9;  % center frequency
c = physconst('LightSpeed');
W = 200e6; % bandwidth
lambda = c/(fc + W); % wavelngth
fs = (2*W)*1.5; % sampling frequency
phi = pi/4; % azimuth
theta = pi/3; % elevation
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
theta_i = pi/8;
u_i = sin(theta_i); % normal vector for determining delays
tau_i = x_pos*u_i/c;
t_i_array = n - tau_i;
t_i = t_i_array(:);
T_i = max(t_i) - min(t_i);
e_i  = repmat(exp(1i*2*pi*fc*(-tau_i+tau)),1,N);

%% MVDR filter design
Rs = exp(-1i*2*pi*fc*(tau - tau')).*sinc(2*W*(tau - tau'));
Rs = Rs/trace(Rs);

Ri = exp(-1i*2*pi*fc*(tau_i - tau_i')).*sinc(2*W*(tau_i - tau_i'));
Ri = Ri/trace(Ri);

Rn = eye(length(tau));

%% Slepian basis generation
T_1 = max(tau)-min(tau);
L = 0;
K =ceil(2*W*T_1)+L;
V_nu = dpss(M,W*T_1,K);
V_mod = bsxfun(@times,exp(-1i*2*pi*fc*tau),V_nu);

SIR_dB = -30;
SNR_dB = -10;
sigma_n = (1/sqrt(M))*10^(-SNR_dB/20);
sigma_i = 10^(-SIR_dB/20);
Rt = Rs + Ri*(sigma_i^2) + Rn*(sigma_n^2);
Rt_inv = Rt^-1;
Phi_full = ((V_mod'*Rt_inv*V_mod)^(-1))*(V_mod'*Rt_inv);

%% measurement design
V = orth(Phi_full');
var_nom = real(trace(((V_mod'*Rt_inv*V_mod)^(-1))) - trace(Rs))
D = 1000;
bits = 1;
thresh = 1e-9;
var_thresh_dB = 1; % vary this to get better or worse quantized measurments

phi = zeros(M,10*K);
phi_orth = zeros(M,10*K);
P = 1;
    
var_prev = inf;
    for ii = 1:100*K
        U = (randn(K,D) + 1i*randn(K,D))/sqrt(2);
        Q = V*U;%-phi_orth*(phi_orth'*V*U);
        kappa = max(sqrt(sum(abs(V).^2,2)))/sqrt(2);
        Q = (Q)./kappa;
        Q = (quantizer_fixed(Q,bits,1,0));
        
        Q_scale = Q./sqrt(sum(abs(Q).^2,1));
        
        [H,S,Wii] = svd(phi_orth(:,1:(P-1))'*V);
        sig_max = max(max(diag(S)));
        sig_min = min(min(diag(S)));
        Ht = H(:,diag(S)>thresh*sig_max);
        Wt = Wii(:,diag(S)>thresh*sig_max);
        Vi = V*Wt;
        G = V'*Q_scale - Wt*((Vi'*Q_scale));
        G_norm = sqrt(sum(abs(G).^2,1));
%         G_norm = G_norm.^2 + (1-min(min(diag(S))))*(sum(abs(V'*Q_scale).^2,1));
        [~,idx] = max(G_norm);
        
        if P>=K
            phi_ii = [phi(:,1:(P-1)) Q(:, idx(1))];
            Re = ((phi_ii'*V_mod)'*(((phi_ii'*Rt)*phi_ii)^-1)*(phi_ii'*V_mod))^-1;
            var_ii = real(trace(Re)-trace(Rs));
            
            if db((var_ii)/var_nom)/2<=var_thresh_dB
                fprintf('Number of %d-bit measurements = %d\n',bits,P)
                phi = [phi(:,1:(P-1)) Q(:, idx(1))];
                break
            end
            
            if var_ii<var_prev        
                phi(:,P) = Q(:, idx(1));
                phi_orth(:,P) = phi(:,P) - phi_orth(:,1:P-1)*(phi_orth(:,1:P-1)'*phi(:,P));
                phi_orth(:,P) = phi_orth(:,P)/norm(phi_orth(:,P));
                P = P+1;
                var_prev = var_ii;
            end
    
        else
            phi(:,P) = Q(:, idx(1));
            if ii==1
                phi_orth(:,P) = phi(:,P)/norm(phi(:,P));
            else
                phi_orth(:,P) = phi(:,P) - phi_orth(:,1:P-1)*(phi_orth(:,1:P-1)'*phi(:,P));
                phi_orth(:,P) = phi_orth(:,P)/norm(phi_orth(:,P));
            end
            P = P+1;
        end
        
    end

Psi = phi'*V_mod;
Rw = (phi'*Rt*phi)^(-1);
Gamma = ((Psi'*Rw*Psi)^-1)*Psi'*Rw;

%% Filter generation
R = 2;
edge_lim = R + ceil(fs*(max(tau)-min(tau))/2);% clip signal to account for edge distortion
comp_idxs = (edge_lim+1):(N-edge_lim-1); % comparison indicies


N_filter = N+2*R+1;
wrapN = @(x, n) (1 + mod(x-1, n));
[fb] = Kernel_Filter_Gen(fs, @Psi_sinc,'R',R);
[h, h_idx] = fb.Filter_Gen(-tau);
h = squeeze(h);
h_idx = squeeze(h_idx);
filter_idx = [min(h_idx(:)):max(h_idx(:))];

H = zeros(M,max( h_idx(:) - min(h_idx(:))+1));
for ii=1:M
    H(ii,h_idx(ii,:) - min(h_idx(:))+1) = h(ii,:);
end

H_buffer = zeros(N_filter,M);
for ii=1:M
    H_buffer(wrapN(filter_idx,N_filter),ii) = H(ii,:);
end



tap_lim = size(H,2);
h = H(:)/M;
t_tap = ones(M,1)*[0:tap_lim-1]/fs;

%% compute beampattern

f_beam = linspace(-W,1*W,10);
theta_beam = linspace(-pi/2,pi/2,2000);
beam_pattern = zeros(length(f_beam),length(theta_beam));
e_tau = repmat(exp(1i*2*pi*(fc)*(tau)),tap_lim,1);
e_mod_tau = exp(-1i*2*pi*fc*(t_tap));
e_mod_tau = e_mod_tau(:);

%% THESE ARE THE EMBEDDED FILTERS
H_embed = (H'*(V_nu*Gamma))'/M;% This is the filter bank, it is P times # of taps
h_embed = H_embed(:);% this is the vecotrized form

%% calculate pattern


for ii = 1:length(f_beam)
    for jj = 1:length(theta_beam)
        tau_jj = -repmat(x_pos*(sin(theta_beam(jj)))/c,1,tap_lim);
        t_jj = (t_tap-tau_jj);
%         t_jj = tau_jj(:);

        e_jj = reshape(exp(1i*2*pi*(f_beam(ii))*t_jj).*exp(1i*2*pi*fc*tau_jj),M,tap_lim);
        e_le_jj = (phi'*e_jj);
        beam_pattern(ii,jj) = (abs(h_embed'*e_le_jj(:)));
    end
end

figure(1)
plot(rad2deg(theta_beam),zeros(length(theta_beam),1),'r--','linewidth',1)
hold on
plot(rad2deg(theta_beam),db(beam_pattern),'k','linewidth',1)
plot(rad2deg([theta,theta]),[-90,10],'g--','linewidth',1)
plot(rad2deg([theta_i,theta_i]),[-90,10],'g--','linewidth',1)
grid on
xlabel('$\theta$ (degrees)','Interpreter','latex')
ylabel('Response(dB)','Interpreter','latex')
legend({'Distortionless response'},'location','northwest','Interpreter','latex')
xlim([-90,90])
ylim([-45,.5])
% saveas(gca,'figs/le_mvdr_ula_beampattern.png')

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