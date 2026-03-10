clear
clc
close all

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

B = 2^5; % no. of beams
M = 2^7; % array size
N = 2^6; % no. of samples
L = ceil(N); % no. of samples in frequency is 2*L
MN = M*N;

Theta = zeros(2*B+1,1);
idx = 40;%randi(B,1);
for i=1:2*B+1
    Theta(i) = asin((i-1-B)/B);
end

fc = 20e9;  % center frequency
Ws = 5e9; % to prevent spatial aliasing
c = physconst('LightSpeed');
lambda = c/(fc+Ws); % wavelngth
W = 5e9; % bandwidth
fs = 2*W; % sampling frequency
T = 1/fs;
theta = 0.8842;%Theta(idx);
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
demod_phase = repmat(exp(-1i*2*pi*fc*tau),N,1);%exp(1i*2*pi*fc*t);

%% Interferer Specs
idx_i = 30;
theta_i = -0.4098;%Theta(idx_i);
u_i = sin(theta_i); % normal vector for determining delays
tau_i = x_pos*u_i/c;
t_i_array = n - tau_i;
t_i = t_i_array(:);
T_aperture_i = max(t_i) - min(t_i);
demod_phase_i  = repmat(exp(-1i*2*pi*fc*(tau_i)),N,1);%repmat(exp(-1i*2*pi*fc*tau_i),N,1); % modulation vector after pre-steering

% %% Interferer Specs
% idx_i2 = 30;
% theta_i2 = pi/6;%Theta(idx_i);
% u_i2 = sin(theta_i2); % normal vector for determining delays
% tau_i2 = x_pos*u_i2/c;
% t_i_array2 = n - tau_i2;
% t_i2 = t_i_array2(:);
% T_aperture_i2 = max(t_i2) - min(t_i2);
% demod_phase_i2  = exp(-1i*2*pi*fc*(t-t_i2));%repmat(exp(-1i*2*pi*fc*tau_i),N,1); % modulation vector after pre-steering

%% construct basis matrix
frequency_samples = (1/(2*T*L))*linspace(-L,L-1,2*L);
F_source = exp(1i*2*pi*t*frequency_samples).*demod_phase;
F_check = F_source;
% F_source = F_source.*demod_phase;
betas_no_null = inv(F_source'*F_source + 1e-5*eye(2*L))*F_source';

F_interference = exp(1i*2*pi*t_i*frequency_samples);
% F_interference2 = exp(1i*2*pi*t_i2*frequency_samples);
F_interference = F_interference.*demod_phase_i;
% F_interference2 = F_interference2.*demod_phase_i2;
beta_null = inv(F_source'*F_source + 1e-8*eye(2*L))*F_source'*(eye(MN,MN) - F_interference*inv(F_interference'*F_interference + 1e-8*eye(2*L))*F_interference');

%% equispaced sampled basis matrix
t_nyq = [0:N-1]/fs; % nyquist samples for testing
psi = exp(1i*2*pi*t_nyq'*frequency_samples);
% for ll=1:2*L
%     f_l = (1/(2*T*L))*(ll-1-L);
%     psi(:,ll) = exp(1i*2*pi*f_l*t_nyq);
% end

%% compute beampattern
W = 5e9; % bandwidth
f_beam = linspace(fc-0.5*W,fc+0.5*W,1);
theta_beam = linspace(-pi/2,pi/2,20000);
beam_pattern_no_null = zeros(length(f_beam),length(theta_beam));
beam_pattern_null = zeros(length(f_beam),length(theta_beam));
e_mod = exp(-1i*2*pi*fc*t);
e_mod_i = exp(1i*2*pi*fc*t_i);

for ii = 1:length(f_beam)
    for jj = 1:length(theta_beam)
        tau_jj = x_pos*sin(theta_beam(jj))/c;
        t_jj = n-tau_jj;
        t_jj =t_jj(:);
        % e_mod  = repmat(exp(-1i*2*pi*(fc)*(tau_jj)),N,1);
        demod_phase_jj = repmat(exp(-1i*2*pi*fc*tau_jj),N,1);%exp(1i*2*pi*fc*t);
        e_jj = exp(1i*2*pi*(f_beam(ii)-fc)*t_jj).*demod_phase_jj;
        beam_pattern_no_null(ii,jj) = norm(psi*(betas_no_null*e_jj))/sqrt(length(t_nyq));
        beam_pattern_null(ii,jj) = norm(psi*(beta_null*e_jj))/sqrt(length(t_nyq));
        %         beam_pattern(ii,jj) = norm(V_nu*(Phi*e_jj))/norm(e_jj);%*norm(h_vec));
    end
end


figure(1)
p0 = plot(rad2deg(theta_beam),zeros(length(theta_beam),1),'r--');
hold on
p1 = plot(rad2deg(theta_beam),db(beam_pattern_no_null),'Color',[0 0 0 0.5]);
hold on
p2 = plot(rad2deg(theta_beam),db(beam_pattern_null),'Color',[0 0 1 0.5]);
grid on
xlabel('$\theta$ (degrees)','Interpreter','latex','FontSize',12)
ylabel('Response(dB)','Interpreter','latex','FontSize',12)

% legendHandles = [repmat(p0,50,1), p1, p2];
% legendLabels = {'Distortionless response','No Nulling', 'Projection Nulling'};
% lgd = legend(legendHandles, legendLabels);
% 
% legendChildren = findobj(lgd, 'Type', 'Line');
% set(legendChildren, 'LineWidth', 2); % Set solid lines for legend
% Dummy plots for solid legend entries with different styles
legend_p0 = plot(-1, -1, 'r--', 'LineWidth', 2);  % Solid red
legend_p1 = plot(-1, -1, 'k-', 'LineWidth', 2); % Dashed blue
legend_p2 = plot(-1, -1, 'b-', 'LineWidth', 2); % Dashed blue

% Legend
legend([legend_p0, legend_p1, legend_p2], {'Distortionless response', 'No interference cancellation', 'Interference nulling'},'location','northwest','Interpreter','latex','FontSize',12);
% legend({'Distortionless response','No Nulling', 'Projection Nulling'},'location','northwest','Interpreter','latex')
xlim([-90,90])
% ylim([-45,.5])
% exportgraphics(gcf, 'beam_pattern_ula.pdf', 'ContentType', 'vector');

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