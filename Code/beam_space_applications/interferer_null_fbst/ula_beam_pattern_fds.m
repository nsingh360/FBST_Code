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

S = log2(M); % stages in the beamformer
R = 2; % 2-radix beamformer
L_fac = 2; % spatial downsampling factor

Theta = zeros(2*B+1,1);
for i=1:2*B+1
    Theta(i) = asin((i-1-B)/B);
end

%% Angle sampling for final stage fast DS
r_lin = linspace(0,R^S-1,R^S);
Theta_ds = asin(1 - (2*r_lin + 1)/R^S);

fc = 20e9;  % center frequency
Ws = 5e9; % to prevent spatial aliasing
c = physconst('LightSpeed');
lambda = c/(fc+Ws); % wavelngth
eta = lambda/(2*c);
W = 5e9; % bandwidth
fs = 2*W; % sampling frequency
T = 1/fs;
idx_theta = 9;
theta = Theta_ds(idx_theta);%Theta(idx); 
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
mod_phase = exp(1i*2*pi*(fc)*tau);

taps = 8; % for ds
edge_lim = taps + ceil(fs*(max(tau)-min(tau))/2);% clip signal to account for edge distortion
comp_idxs = (edge_lim+1):(N-edge_lim-1); % comparison indicies

taps2 = 8; % for fast ds
edge_lim2 = taps2 + ceil(fs*(max(tau)-min(tau))/2);% clip signal to account for edge distortion
comp_idxs2 = (edge_lim2+1):(N-edge_lim2-1); % comparison indicies

%% Interferer Specs
idx_1_i = 174;%61;%239;
idx_2_i = 175;
alpha = 0.5;
theta_i = asin(alpha*sin(Theta(idx_1_i)) + (1-alpha)*sin(Theta(idx_2_i)));

idx_2 = 42;%239;
idx_1 = 41;
% idx_i = B;
% theta_i = Theta(idx_i);
u_i = sin(theta_i); % normal vector for determining delays
tau_i = x_pos*u_i/c;
t_i_array = n - tau_i;
t_i = t_i_array(:);
T_aperture_i = max(t_i) - min(t_i);
demod_phase_i  = repmat(exp(-1i*2*pi*fc*tau_i),N,1);%exp(-1i*2*pi*fc*(t-t_i)); % modulation vector after pre-steering
mod_phase_i = exp(1i*2*pi*(fc)*tau_i);

%% construct basis matrix
delta=1e-5;
frequency_samples = (1/(2*T*L))*linspace(-L,L-1,2*L);
F = exp(1i*2*pi*t*frequency_samples).*demod_phase;
toeplitz_inv = inv(F'*F + delta*eye(2*L));
F_i = exp(1i*2*pi*t_i*frequency_samples).*demod_phase_i;
A_eplitz_i = F_i'*F_i;
null_proj = F_i*inv(A_eplitz_i + delta*eye(2*L))*F_i';
% F_source = F_source.*demod_phase;

%% equispaced sampled basis matrix
t_nyq = [0:N-1]/fs;
t_nyq2 = [0:N-1]/fs + 1/fs; % nyquist samples for testing
t_nyq3 = [0:N-1]/fs + S/fs;
psi_ds = exp(1i*2*pi*t_nyq2'*frequency_samples);
psi_fds = exp(1i*2*pi*t_nyq3'*frequency_samples);

%% compute beampattern
f_beam = linspace(fc-0.5*W,fc+0.5*W,20);
theta_beam = linspace(-pi/2,pi/2,4000);

beam_pattern_ds_nulled = zeros(length(f_beam),length(theta_beam));
beam_pattern_ds_nulled_beamspace = zeros(length(f_beam),length(theta_beam));
beam_pattern_fds_nulled_beamspace = zeros(length(f_beam),length(theta_beam));

% e_mod = exp(-1i*2*pi*fc*t);
% e_mod_i = exp(-1i*2*pi*fc*t_i);

% e_jj_s = exp(1i*2*pi*(f_beam(1)-fc)*t).*demod_phase;
% e_jj_i = exp(1i*2*pi*(f_beam(1)-fc)*t_i).*demod_phase_i;

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

        G = toeplitz_inv*F'*F_i;
        G_ds = psi_ds*G*pinv(psi_ds);
        G_fds = psi_fds*G*pinv(psi_fds);

        % y_filter_ttd_nulled = TTD_Beamformer((y_array_nulled.*mod_phase).',-tau, fs, @Psi_sinc, 'R',taps).';
        % y_ds_nulled = sum(y_filter_ttd_nulled(:,1:N),1)/M;

        y_filter_ttd = TTD_Beamformer((y_array.*mod_phase).',-tau, fs, @Psi_sinc, 'R',taps).';
        y_filter_ttd_i = TTD_Beamformer((y_array.*mod_phase_i).', -tau_i, fs, @Psi_sinc, 'R', taps).';
        y_ds_s = sum(y_filter_ttd(:,1:N),1)/M;
        y_ds_i = sum(y_filter_ttd_i(:,1:N),1)/M;
        y_ds_nyq = y_ds_s.' - G_ds*(y_ds_i.');

        %% fast delay and sum (for source beam)
        v_now = y_array.*mod_phase;
        for stage=1:S
            v_prev = v_now;
            num_virtual_sensors = M/2^stage;
            num_sectors = R^stage;
            v_now = zeros(num_virtual_sensors,N,num_sectors); % init partial beamforming at stage s

            for mm=0:num_virtual_sensors-1
                for rr = 0:num_sectors-1
                    m = mm+1;
                    r = floor((rr)/R)+1;
                    delays = ((-1)^(rr))*(1/(R^stage))*eta*[2^(stage-1)*(2*mm+0.5)-M/2; 2^(stage-1)*(2*mm+1.5)-M/2];
                    y_interp = reshape(v_prev([2*mm+1,2*mm+2],:,r),[2,N]);
                    v_prev_interped = TTD_Beamformer(y_interp.',-delays, fs, @Psi_sinc, 'R',taps).';
                    v_now(m,:,rr+1) = sum(v_prev_interped,1)./2;
                end
            end
        end
        % v_now = squeeze(v_now);
        y_fds_s = v_now(1,:,idx_theta).';

        %% fast delay and sum (for interferer beam)
        v_now = y_array.*mod_phase_i;
        for stage=1:S
            v_prev = v_now;
            num_virtual_sensors = M/2^stage;
            num_sectors = R^stage;
            v_now = zeros(num_virtual_sensors,N,num_sectors); % init partial beamforming at stage s

            for mm=0:num_virtual_sensors-1
                for rr = 0:num_sectors-1
                    m = mm+1;
                    r = floor((rr)/R)+1;
                    delays = ((-1)^(rr))*(1/(R^stage))*eta*[2^(stage-1)*(2*mm+0.5)-M/2; 2^(stage-1)*(2*mm+1.5)-M/2];
                    y_interp = reshape(v_prev([2*mm+1,2*mm+2],:,r),[2,N]);
                    v_prev_interped = TTD_Beamformer(y_interp.',-delays, fs, @Psi_sinc, 'R',taps).';
                    v_now(m,:,rr+1) = sum(v_prev_interped,1)./2;
                end
            end
        end
        %------ create interpolating vectors in time.
        theta_interps = Theta_ds(idx_1-3:idx_2+3);
        y_nyqs_interps = zeros(length(theta_interps),length(t_nyq3));
        count_indx = -3;
        for it=1:length(theta_interps)
            y_nyqs_interps(it,:) = v_now(1,:,idx_1+count_indx).';
            count_indx = count_indx + 1;
        end
        y_fds_i_interp = interp1(Theta_ds(idx_1-3:idx_2+3),y_nyqs_interps,theta_i,'spline').';
        % v_now = squeeze(v_now);
        y_fds_i = v_now(1,:,idx_1).';
        y_fds_nyq = y_fds_s - G_fds*(y_fds_i_interp); 
        % beam_pattern_ds_nulled(ii,jj) = norm(y_ds_nulled(comp_idxs))/sqrt(length(comp_idxs));
        beam_pattern_ds_nulled_beamspace(ii,jj) = norm(y_ds_nyq(comp_idxs))/sqrt(length(comp_idxs));
        beam_pattern_fds_nulled_beamspace(ii,jj) = norm(y_fds_nyq(comp_idxs))/sqrt(length(comp_idxs));
    end
end

figure(1)
p0 = plot(rad2deg(theta_beam),zeros(length(theta_beam),1),'r--');
hold on
grid on
p1 = plot(rad2deg(theta_beam),db(beam_pattern_ds_nulled_beamspace),'Color',[0 0 0 1],'LineWidth',1);
hold on
grid on
p2 = plot(rad2deg(theta_beam),db(beam_pattern_fds_nulled_beamspace),'Color',[0 0 1 0.5],'LineWidth',1);
xlabel('$\theta$ (degrees)','Interpreter','latex','FontSize',12)
ylabel('Response(dB)','Interpreter','latex','FontSize',12)
legend_p0 = plot(-1, -1, 'r--', 'LineWidth', 2);  % Solid red
legend_p1 = plot(-1, -1, 'k-', 'LineWidth', 2); % Dashed blue
legend_p2 = plot(-1, -1, 'b-', 'LineWidth', 2); % Dashed blue
% Legend
legend([legend_p0, legend_p1, legend_p2], {'Distortionless response', 'TTD Nulled beamspace', 'FDS Nulled Beam space'},'location','south','Interpreter','latex','FontSize',12);
% xlim([min(theta_beam),max(theta_beam)])
% ylim([-45,.5])
% y1 = ylim;
% ylim([-30, 1.5]);
% exportgraphics(gcf, 'super_res_ula_beam_pattern.pdf', 'ContentType', 'vector')
% legend({'Distortionless response','No Nulling', 'Projection Nulling'},'location','northwest','Interpreter','latex')

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