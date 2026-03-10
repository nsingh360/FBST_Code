clear
clc
close all

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

%% array spatial and temporal specificiations
M = 2^6; % array size
N = 2^5; % temporal samples
MN = M*N;
fc = 20e9;  % center frequency
c = physconst('LightSpeed');

W = 5e9; % bandwidth
lambda = c/(fc+W); % wavelngth
fs = 2*W; % sampling frequency
phi = pi/4; % azimuth
theta = pi/3; % elevation
m = [(-M/2+1/2):(M/2-1/2)]';
x_pos = m*lambda/2; % sensor positions
u_s = sin(theta); % normal vector for determining delays


%% Signal specs that do not need to be redifined in loop
n_sinusoids = 100; % number of sinusoids in signal
n = [0:N-1]/fs; % temporal sample vectors
tau = x_pos*u_s/c; % relative delays to phase center
t_array = ones(M,1)*n; % delays across array and time
t = t_array(:); %
T = max(t)-min(t);% temporal aperture of the array

%% Interferer specs that do not need to be redfined in loop
phi_i = -pi/4;
theta_i = pi/6;
u_i = sin(theta_i); % normal vector for determining delays
tau_i = x_pos*u_i/c;
t_i_array = n - tau_i + tau;
t_i = t_i_array(:);
T_i = max(t_i) - min(t_i);
e_i  = repmat(exp(1i*2*pi*fc*(-tau_i+tau)),1,N);

%% MVDR filter design
K =32;
KM = M*K;
N_filter = K + N + 1; % filter length
t_mvdr = ones(M,1)*n(1:K);
t_mvdr = t_mvdr(:);
T_mvdr = max(t_mvdr)-min(t_mvdr);
t_i_mvdr = t_i_array(:,1:K); % temporal points over which interferer is seen
t_i_mvdr =t_i_mvdr(:);
T_i_mvdr = max(t_i_mvdr)-min(t_i_mvdr);
L_s = 1;
L_i = 10;



C = kron(eye(K),ones(M,1)); % constraint matrix
F = eye(K,1); % response vector

delta = 1e-5;

%% simulation specifications
bits = 8;
trials = 50;
SIR_dB = -30;
SNR_dB = 0;


sigma_n = (1/sqrt(MN))*10^(-SNR_dB/20);
sigma_i = 10^(-SIR_dB/20);


t_nyq = [0:N-1]'/fs; % nyquist samples for testing

n_sig = [-2*N:4*N]/(2*W);

B = sinc(2*W*(t-n_sig));
B_scale = 1/sqrt(trace(B'*B));
B = B*B_scale;

Bi = bsxfun(@times,e_i(:),sinc(2*W*(t_i-n_sig)));
Bi_scale = 1/sqrt(trace(Bi'*Bi));
Bi = Bi*Bi_scale;

B_nyq = sinc(2*W*(t_nyq-n_sig))*B_scale;


Rn = (sigma_n^2)*eye(MN);
Rs = (B*B');
Ri = (sigma_i^2)*(Bi*Bi');

R = Rs + Rn + Ri;

R_inv = R^-1;

h = (R_inv*C)*((C'*R_inv*C)^-1)*F;
[h,~,~] = quantizer(h,bits);


t_tap = [0:K-1]/fs;

%% compute beampattern

f_beam = linspace(-W,1*W,50)*.5;
% theta_beam = linspace(-pi/2,pi/2,2000);
theta_beam = linspace(-pi/2,pi/2,2000);
% sin_theta_beam = sin(theta_beam)+sin(theta);
% sin_theta_beam = sin_theta_beam(sin_theta_beam<=1 & sin_theta_beam>=-1);
% theta_beam = asin(sin_theta_beam);
beam_pattern = zeros(length(f_beam),length(theta_beam));
e_mod = exp(-1i*2*pi*(fc)*(ones(M,1)*t_tap));
e_mod =e_mod(:);
% h = exp(1i*2*pi*fc*t_i);

for ii = 1:length(f_beam)
    for jj = 1:length(theta_beam)
        tau_jj = x_pos*(sin(theta_beam(jj))-sin(theta))/c;
        t_jj = t_tap-tau_jj;
        t_jj = t_jj(:);
%         e_mod  = repmat(exp(1i*2*pi*(fc+f_beam(ii))*(tau_jj)),2*K,1);
        e_jj = exp(1i*2*pi*(fc + f_beam(ii))*t_jj).*e_mod;
        beam_pattern(ii,jj) = abs(e_jj'*h);%/(norm(e_jj)*norm(h_vec));
    end
end
save('mpdr_beam.mat','beam_pattern','theta_beam')

figure(1)
plot(rad2deg(theta_beam),zeros(length(theta_beam),1),'r--')
hold on
plot(rad2deg(theta_beam),db(beam_pattern),'k')
grid on
xlabel('$\theta$ (degrees)','Interpreter','latex')
ylabel('Response(dB)','Interpreter','latex')
legend({'Distortionless response'},'location','northwest','Interpreter','latex')
xlim([-90,90])
ylim([-45,max(max(db(beam_pattern(:))),.5)])
saveas(gca,'figs/mvdr_ula_beampattern_8_bit.png')




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