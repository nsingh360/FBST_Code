clear
clc
close all

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

B = 2^7;
M = 2^7; % array size
N = 2^6; % no. of samples
L = ceil(N); % no. of samples in frequency is 2*L
MN = M*N;

Theta = zeros(2*B+1,1);
idx = 240;%120;%randi(B,1);
for i=1:2*B+1
    Theta(i) = asin((i-1-B)/B);
end

fc = 20e9;  % center frequency
Ws = 5e9; % to prevent spatial aliasing
c = physconst('LightSpeed');
lambda = c/(fc+Ws); % wavelngth
eta = lambda/(2*c);
W = 5e9; % bandwidth
fs = 2*W; % sampling frequency
T = 1/fs;
theta = Theta(idx);%Theta(idx); 
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
demod_phase = repmat(exp(-1i*2*pi*fc*tau),N,1);%exp(1i*2*pi*fc*t);
mod_phase = exp(1i*2*pi*(fc)*tau);

%% Interferer Specs
%% select angle to super resolve 
idx_1 = 174;%61;%239;
idx_2 = 175;
alpha = 0.5;
theta_i = asin(alpha*sin(Theta(idx_1)) + (1-alpha)*sin(Theta(idx_2)));
% idx_i = B;
% theta_i = Theta(idx_i);
u_i = sin(theta_i); % normal vector for determining delays
tau_i = x_pos*u_i/c;
t_i_array = n - tau_i;
t_i = t_i_array(:);
T_aperture_i = max(t_i) - min(t_i);
demod_phase_i  = repmat(exp(-1i*2*pi*fc*tau_i),N,1); % modulation vector after pre-steering
mod_phase_i = exp(1i*2*pi*(fc)*tau_i);

%% construct the toeplitz matrix
delta = 1e-5;
freq_samples = 1/(2*T*L)*linspace(-L,L-1,2*L);
F = exp(1i*2*pi*t*freq_samples);
A_eplitz = F'*F;
toeplitz_inv = inv(A_eplitz + (delta)*eye(2*L));
F = F.*demod_phase;

%% construct fourier and toeptliz matrix for interferer
F_i = exp(1i*2*pi*t_i*freq_samples);
A_eplitz_i = F_i'*F_i;
F_i = F_i.*demod_phase_i;

%% equispaced sampled basis matrix
t_nyq = min(t):T/1.1:max(t); % nyquist samples for testing
psi = exp(1i*2*pi*t_nyq'*freq_samples);
t_nyq3 = [0:N-1]/fs;
psi_nyq = exp(1i*2*pi*t_nyq3'*freq_samples);


%% compute beampattern
f_beam = linspace(fc-0.5*W,fc+0.5*W,20);
theta_beam = linspace(-pi/2,pi/2,4000);

beam_pattern_fbst_nulled = zeros(length(f_beam),length(theta_beam));

for ii = 1:length(f_beam)
    for jj = 1:length(theta_beam)
        fprintf('Frequency: %.2f, Angle: %.2f\n',f_beam(ii),theta_beam(jj));
        tau_jj = x_pos*sin(theta_beam(jj))/c;
        t_jj = n-tau_jj;
        t_jj =t_jj(:);
        % e_mod  = repmat(exp(-1i*2*pi*(fc)*(tau_jj)),N,1);
        demod_phase_jj = repmat(exp(-1i*2*pi*fc*tau_jj),N,1);%exp(1i*2*pi*fc*t);
        % mod_phase_jj = exp(1i*2*pi*(fc)*tau_jj);
        e_jj = exp(1i*2*pi*(f_beam(ii)-fc)*t_jj).*demod_phase_jj;

        y = e_jj; % for interpolation using chirp-z
        % y_nulled = y - null_proj*y;
        % y_array_nulled = reshape(y_nulled, [M,N]);
        y_array = reshape(y, [M,N]);

        %% Nulling in beam space
        [M,N] = size(y_array);
        %--- chirp-z transform
        F_transform = zeros(M,2*L);
        A = -1;%exp(-1i*2*pi*T*W);
        W_czt = exp(-1i*pi/L);
        for mm=1:M
            F_transform(mm,:) = czt(y_array(mm,:).',2*L,W_czt,A);
        end
        X_czt = zeros(2*L,2*B+1);
        A= 1;
        B_list = reshape(linspace(0,2*B,2*B+1), [1,2*B+1]);
        for ll=1:2*L
            l_dash = ll-1;
            A = exp(1i*pi*(1/B)*(fc_norm + (1/(2*fc_tilde*T*L))*(l_dash-L))*B);
            W_czt = exp(1i*pi*(1/B)*(fc_norm + (1/(2*fc_tilde*T*L))*(l_dash-L)));
            coeff1 = exp(-1i*pi*(1/B)*(fc_norm + (1/(2*fc_tilde*T*L))*(l_dash-L))*(M/2-1/2)*B_list); % adjusting for array center
            coeff2 = exp(1i*pi*(1/B)*(fc_norm + (1/(2*fc_tilde*T*L))*(l_dash-L))*(M/2-1/2)*B); % adjusting for negative angles
            X_czt(ll,:) = czt(F_transform(:,ll),2*B+1,W_czt,A);
            X_czt(ll,:) = coeff2*X_czt(ll,:).*coeff1;
        end
        w_l = X_czt(:,idx);
        w_l_i = X_czt(:,idx_1);

        beta_s = toeplitz_inv*w_l;
        beta_i = inv(F_i'*F_i + delta*eye(2*L))*w_l_i;
        G = toeplitz_inv*F'*F_i;
        G_ = psi*G*pinv(psi);

        %------ create interpolating vectors in time.
        theta_interps = Theta(idx_1-3:idx_2+3);
        y_nyqs_interps = zeros(length(theta_interps),length(t_nyq));
        count_indx = -3;
        for it=1:length(theta_interps)
            theta_interp = theta_interps(it);
            t_array_interp = n - x_pos*sin(theta_interp)/c;
            t_interp = t_array_interp(:);
            F_interp = exp(1i*2*pi*t_interp*freq_samples);
            y_nyqs_interps(it,:) = psi*inv(F_interp'*F_interp + delta*eye(2*L))*X_czt(:,idx_1+count_indx);
            count_indx = count_indx + 1;
        end
        y_nyq_i_interp = interp1(Theta(idx_1-3:idx_2+3),y_nyqs_interps,theta_i,'spline').';

        y_nyq_s = psi*beta_s;
        y_nyq_i = psi*beta_i;
        y_nyq_time = y_nyq_s - G_*y_nyq_i_interp;
        y_nyq_time = y_nyq_time(t_nyq<=(N-1)*T & t_nyq>=0);
        beam_pattern_fbst_nulled(ii,jj) = norm(y_nyq_time)/sqrt(length(y_nyq_time));

    end
end


figure(1)
p0 = plot(rad2deg(theta_beam),zeros(length(theta_beam),1),'r--');
hold on
grid on
p1 = plot(rad2deg(theta_beam),db(beam_pattern_fbst_nulled),'Color',[0 0 0 1],'LineWidth',1);