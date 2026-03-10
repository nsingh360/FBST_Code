clear
clc
close all

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

B = 2^3; % no. of beams
M = 2^7; % array size
N = 2^6; % no. of samples
L = ceil(N); % no. of samples in frequency is 2*L
MN = M*N;

%% Theta grid stores the angles available
Theta = zeros(2*B+1,1);
for i=1:2*B+1
    Theta(i) = asin((i-1-B)/B);
end

fc = 20e9;  % center frequency
Ws = 5e9; % to prevent spatial aliasing
c = physconst('LightSpeed');
lambda = c/(fc+Ws); % wavelngth
W = 5e9; % bandwidth
fs = 2.5*W; % sampling frequency
T = 1/fs;
m = [(-M/2+1/2):(M/2-1/2)]';
x_pos = m*lambda/2; % sensor positions

interference_bias = zeros(2*B+1,2*B+1);

for i=1:2*B+1

    theta = Theta(i); % source angle
    samples = 1000;
    %% source signal specs
    u_s = sin(theta);
    n_sinusoids = 100; % number of sinusoids in signal
    n = [0:N-1]/fs; % temporal sample vectors
    tau = x_pos*u_s/c; % relative delays to phase center
    t_array = n-tau; % delays across array and time
    t = t_array(:); %
    t = t - min(t);
    T_aperture = max(t)-min(t);% temporal aperture of the array
    demod_phase = repmat(exp(-1i*2*pi*fc*tau),N,1);

    frequency_samples = (1/(2*T*L))*linspace(-L,L-1,2*L);
    F = exp(1i*2*pi*t*frequency_samples); % fourier matrix for source
    F = F.*demod_phase;
    P = inv(F'*F + 1e-8*eye(2*L))*F';
    
    t_eq = linspace(min(t),max(t),samples);
    K = ceil(2*W*T_aperture) + 8;
    [V_u,Lambda] = dpss(samples,W*T_aperture,K);
    V_nu = interp1(t_eq,V_u,t,'pchip',0); % non-uniform
    V_nu = V_nu.*demod_phase;

    for j=1:2*B+1

        theta_i = Theta(j); % interferer angle

        fprintf("Source angle: %d, Interferer angle: %d \n",theta*180/pi,theta_i*180/pi);

        %% interferer signal specs
        u_i = sin(theta_i); % normal vector for determining delays
        tau_i = x_pos*u_i/c;
        t_i_array = n - tau_i;
        t_i = t_i_array(:);
        T_aperture_i = max(t_i) - min(t_i);
        t_i = t_i - min(t_i);
        demod_phase_i  = repmat(exp(-1i*2*pi*fc*tau_i),N,1); % modulation vector after pre-steering

        F_i = exp(1i*2*pi*t_i*frequency_samples); % fourier matrix for interferer
        F_i = F_i.*demod_phase_i;
        P_i = inv(F_i'*F_i + 1e-8*eye(2*L))*F_i';
        % [F_i_orth,~] = qr(F_i,0);
        Proj_i = F_i*P_i;%F_i_orth*F_i_orth';

        %% interference bias evaluation
        K_i = ceil(2*W*T_aperture_i) + 8;
        [V_u_i,Lambda_i] = dpss(samples, W*T_aperture_i, K_i);
        V_nu_i = interp1(t_eq,V_u_i, t_i,'pchip',0);
        V_nu_i = V_nu_i.*demod_phase_i;

        F_eq = exp(1i*2*pi*t_eq'*frequency_samples);

        %term1 = F_eq*P*Proj_i*V_nu;
        term1 = F_eq*P*(eye(MN,MN) - Proj_i)*V_nu - V_u;
        term2 = F_eq*P*(eye(MN,MN) - Proj_i)*V_nu_i;

        interference_bias(i,j) = trace(term1'*term1*diag(Lambda/(2*W*T_aperture))) + trace(term2'*term2*diag(Lambda_i/(2*W*T_aperture_i)));
    end
end


figure;
imagesc(10*log10(interference_bias)); % Display matrix as a heat map
c = colorbar; % Add color scale
colormap(jet); % Choose a colormap (e.g., 'jet', 'parula', 'hot')

% Label the colorbar
ylabel(c, 'dB', 'FontSize', 10);

% Define Tick Positions and Labels
theta_values = -90:30:90; % Define tick locations
xticks(linspace(1, size(interference_bias, 2), numel(theta_values))); % Map to matrix index
yticks(linspace(1, size(interference_bias, 1), numel(theta_values)));
xticklabels(theta_values); % Assign labels
yticklabels(theta_values);

% Labels and title
xlabel('$\theta$','Interpreter','latex','FontSize',14,'FontWeight','bold');
ylabel('$\theta_I$','Interpreter','latex','FontSize',14,'FontWeight','bold');