clear
clc
close all

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

B = 2^7; % no. of beams
M = 2^6; % array size
N = 2^5; % no. of samples
L = ceil(N); % no. of samples in frequency is 2*L
exp_fac = 1.01; % factor to include freqs slightly out of band;
MN = M*N;
SNR_dB = -10; % signal to noise ratio
SIR_dB = [-30,-20]; % signal to interferer ratio

Theta = zeros(2*B+1,1);
idx = 25;%randi(B,1);
for i=1:2*B+1
    Theta(i) = asin((i-1-B)/B);
end

Theta_vals = [0.5*Theta(ceil(B/4.5))+0.5*Theta(ceil(B/4.5)+1),0.5*Theta(B)+0.5*Theta(B+1),0.5*Theta(ceil(3*B/2+B/4))+0.5*Theta(ceil(3*B/2+B/4)+1)];

fc = 20e9;  % center frequency
Ws = 5e9; % to prevent spacial aliasing
c = physconst('LightSpeed');
lambda = c/(fc+Ws); % wavelngth
eta = lambda/(2*c);
W = 5e9; % bandwidth
fs = 2*W; % sampling frequency
T = 1/fs;
theta = Theta_vals(1);
m = [(-M/2+1/2):(M/2-1/2)]';
x_pos = m*lambda/2; % sensor positions
u_s = sin(theta); % normal vector for determining delays
fc_tilde = (fc+Ws); % adjustment for spatial aliasing
fc_norm = fc/fc_tilde;

%% Signal specs that do not need to be redifined in loop
n_sinusoids = 100; % number of sinusoids in signal
n = [0:N-1]/fs; % temporal sample vectors
tau = x_pos*u_s/c; % relative delays to phase center
t_array = n-tau; % delays across array and time
t = t_array(:); %
T_aperture = max(t)-min(t);% temporal aperture of the array
demod_phase = repmat(exp(-1i*2*pi*(fc)*tau),N,1);

%% Interferer Specs
%% select angle to super resolve 
idx_1 = 90;%61;%239;
theta_i = Theta_vals(2);
% idx_i = B;
% theta_i = Theta(idx_i);
u_i = sin(theta_i); % normal vector for determining delays
tau_i = x_pos*u_i/c;
t_i_array = n - tau_i;
t_i = t_i_array(:);
T_aperture_i = max(t_i) - min(t_i);
demod_phase_i  = repmat(exp(-1i*2*pi*fc*tau_i),N,1); % modulation vector after pre-steering

theta_i_low = Theta_vals(2) - pi/120;
% idx_i = B;
% theta_i = Theta(idx_i);
u_i_low = sin(theta_i_low); % normal vector for determining delays
tau_i_low = x_pos*u_i_low/c;
t_i_array_low = n - tau_i_low;
t_i_low = t_i_array_low(:);
T_aperture_i_low = max(t_i_low) - min(t_i_low);
demod_phase_i_low  = repmat(exp(-1i*2*pi*fc*tau_i_low),N,1); % modulation vector after pre-steering

theta_i_high = Theta_vals(2) + pi/120;
% idx_i = B;
% theta_i = Theta(idx_i);
u_i_high = sin(theta_i_high); % normal vector for determining delays
tau_i_high = x_pos*u_i_high/c;
t_i_array_high = n - tau_i_high;
t_i_high = t_i_array_high(:);
T_aperture_i_high = max(t_i_high) - min(t_i_high);
demod_phase_i_high  = repmat(exp(-1i*2*pi*fc*tau_i_high),N,1); % modulation vector after pre-steering

%% select angle to super resolve 
idx_1 = 90;%61;%239;
theta_i2 = Theta_vals(3);
% idx_i = B;
% theta_i = Theta(idx_i);
u_i2 = sin(theta_i2); % normal vector for determining delays
tau_i2 = x_pos*u_i2/c;
t_i_array2 = n - tau_i2;
t_i2 = t_i_array2(:);
T_aperture_i2 = max(t_i2) - min(t_i2);
demod_phase_i2  = repmat(exp(-1i*2*pi*fc*tau_i2),N,1); % modulation vector after pre-steering

theta_i2_low = Theta_vals(3) - pi/120;
% idx_i = B;
% theta_i = Theta(idx_i);
u_i2_low = sin(theta_i2_low); % normal vector for determining delays
tau_i2_low = x_pos*u_i2_low/c;
t_i2_array_low = n - tau_i2_low;
t_i2_low = t_i2_array_low(:);
T_aperture_i2_low = max(t_i2_low) - min(t_i2_low);
demod_phase_i2_low  = repmat(exp(-1i*2*pi*fc*tau_i2_low),N,1); % modulation vector after pre-steering

theta_i2_high = Theta_vals(3) + pi/120;
% idx_i = B;
% theta_i = Theta(idx_i);
u_i2_high = sin(theta_i2_high); % normal vector for determining delays
tau_i2_high = x_pos*u_i2_high/c;
t_i2_array_high = n - tau_i2_high;
t_i2_high = t_i2_array_high(:);
T_aperture_i2_high = max(t_i2_high) - min(t_i2_high);
demod_phase_i2_high  = repmat(exp(-1i*2*pi*fc*tau_i2_high),N,1); % modulation vector after pre-steering

%% correlation matrix
n_sig = [-4*N:4*N]/(2*W);
e_s  = repmat(exp(1i*2*pi*fc*(-tau)),N,1);
Bs = sinc(2*W*(t-n_sig));
Bs_scale = 1/sqrt(trace(Bs'*Bs));
Bs = Bs*Bs_scale;
Rs = (Bs*Bs');

e_i_low  = repmat(exp(1i*2*pi*fc*(-tau_i_low+tau)),N,1);
e_i_high  = repmat(exp(1i*2*pi*fc*(-tau_i_high+tau)),N,1);
e_i2_low  = repmat(exp(1i*2*pi*fc*(-tau_i2_low+tau)),N,1);
e_i2_high  = repmat(exp(1i*2*pi*fc*(-tau_i2_high+tau)),N,1);

e_i  = repmat(exp(1i*2*pi*fc*(-tau_i+tau)),N,1);
e_i2  = repmat(exp(1i*2*pi*fc*(-tau_i2+tau)),N,1);
Bi = bsxfun(@times,e_i,sinc(2*W*(t_i-n_sig)));
Bi_scale = 1/sqrt(trace(Bi'*Bi));
Bi = Bi*Bi_scale;

Bi2 = bsxfun(@times,e_i2,sinc(2*W*(t_i2-n_sig)));
Bi2_scale = 1/sqrt(trace(Bi2'*Bi2));
Bi2 = Bi2*Bi2_scale;

sigma_n = (1/sqrt(MN))*10^(-SNR_dB/20);
sigma_i = 10^(-SIR_dB(1)/20);
sigma_i2 = 10^(-SIR_dB(2)/20);
Rn = (sigma_n^2)*eye(MN);
Ri = (sigma_i^2)*(Bi*Bi') +(sigma_i2^2)*(Bi2*Bi2') ;

%% construct the toeplitz matrix
delta = 1e-5;
freq_samples = exp_fac/(2*T*L)*linspace(-L,L,2*L+1);
F_basis = exp(1i*2*pi*t*freq_samples);
% for ll=1:2*L
%     for kk=1:2*L
%         f_l = (1/(2*T*L))*(ll-1-L);
%         f_k = (1/(2*T*L))*(kk-1-L);
%         A_eplitz(ll,kk) = sum(exp(1i*2*pi*(f_k - f_l)*t));
%     end
% end

% noise covariance
R = Rs + Ri + Rn; % total covariance
Rinv = R^(-1);
Phi_f = (pinv(F_basis'*Rinv*F_basis + 1e-4*eye(2*L+1)))*F_basis'*Rinv;

%% beamspace mvdr (both source and interfer beam)
F_basis_i_low = exp(1i*2*pi*t_i_low*freq_samples).*e_i_low;
F_basis_i2_low = exp(1i*2*pi*t_i2_low*freq_samples).*e_i2_low;
F_basis_i_high = exp(1i*2*pi*t_i_high*freq_samples).*e_i_high;
F_basis_i2_high = exp(1i*2*pi*t_i2_high*freq_samples).*e_i2_high;
beam_space_map = [F_basis'; F_basis_i_low';F_basis_i_high';F_basis_i2_low';F_basis_i2_high'];
beam_space_map = ((beam_space_map*beam_space_map')^(-0.5))*beam_space_map;
Rbeam = beam_space_map*R*beam_space_map';
Rs_inv = inv(Rbeam + 1e-5*eye(size(beam_space_map,1)));
F_square = beam_space_map*F_basis;
Phi_beamspace_f = inv(F_square'*Rs_inv*F_square + 1e-5*eye(2*L+1))*F_square'*Rs_inv*beam_space_map;

%% beamspace mvdr (only source beam)
F_basis_i = exp(1i*2*pi*t_i*freq_samples).*e_i;
beam_space_map = [F_basis'];
beam_space_map = ((beam_space_map*beam_space_map')^(-0.5))*beam_space_map;
Rbeam = beam_space_map*R*beam_space_map';
Rs_inv = inv(Rbeam + 1e-5*eye(2*L+1));
F_square = beam_space_map*F_basis;
Phi_beamspace_source_f = inv(F_square'*Rs_inv*F_square + 1e-4*eye(2*L+1))*F_square'*Rs_inv*beam_space_map;

%% equispaced sampled basis matrix
t_nyq = [0:N-1]/fs; % nyquist samples for testing
psi = exp(1i*2*pi*t_nyq'*freq_samples);

%% compute beampattern
W = 5e9; % bandwidth
f_beam = linspace(-W,1*W,5)*(.5);
theta_beam = linspace(-pi/2,pi/2,2000);
beam_pattern = zeros(length(f_beam),length(theta_beam));
beam_pattern_beamspace = zeros(length(f_beam),length(theta_beam));
beam_pattern_beamspace_2 = zeros(length(f_beam),length(theta_beam));
e_mod = exp(-1i*2*pi*fc*t);

for ii = 1:length(f_beam)
    for jj = 1:length(theta_beam)
        tau_jj = x_pos*sin(theta_beam(jj))/c;
        t_jj = n-tau_jj;
        t_jj =t_jj(:);
        %         e_mod  = repmat(exp(1i*2*pi*(fc+f_beam(ii))*(tau_jj)),2*K,1);
        e_jj = exp(1i*2*pi*(fc + f_beam(ii))*t_jj).*e_mod;
        beam_pattern(ii,jj) = norm(psi*(Phi_f*(e_jj)))/sqrt(length(t_nyq));
        beam_pattern_beamspace(ii,jj) = norm(psi*(Phi_beamspace_f*(e_jj)))/sqrt(length(t_nyq));
        beam_pattern_beamspace_2(ii,jj) = norm(psi*(Phi_beamspace_source_f*(e_jj)))/sqrt(length(t_nyq));
        %         beam_pattern(ii,jj) = norm(V_nu*(Phi*e_jj))/norm(e_jj);%*norm(h_vec));
    end
end


BP_dB = db(beam_pattern.');
BP_dB_beamspace = db(beam_pattern_beamspace.');
BP_dB_beamspace2 = db(beam_pattern_beamspace_2.');

figure(1)
plot(rad2deg(theta_beam),BP_dB)
xlabel('$\theta$','Interpreter','latex')
ylabel('Magnitude Response (dB)','Interpreter','latex')
xlim([-90,90])
sgtitle('MVDR Array space');
% 
figure(2)
plot(rad2deg(theta_beam),BP_dB_beamspace)
xlabel('$\theta$','Interpreter','latex')
ylabel('Magnitude Response (dB)','Interpreter','latex')
xlim([-90,90])
sgtitle('MVDR Beamspace');
exportgraphics(gcf, 'beam_space.pdf', 'ContentType', 'image');

figure(3)
plot(rad2deg(theta_beam),BP_dB_beamspace2)
xlabel('$\theta$','Interpreter','latex')
ylabel('Magnitude Response (dB)','Interpreter','latex')
xlim([-90,90])
sgtitle('MVDR Beamspace (Source)');

            
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