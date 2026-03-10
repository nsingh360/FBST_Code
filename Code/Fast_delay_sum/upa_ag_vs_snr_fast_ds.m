clear
clc
close all

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

S = 4; % no. of stages in the FDS
M = 2^(2*S); % sqrt(M) x sqrt(M) dimensional planar array
N = 2^6;
L = 2; % downsampling factor (in spatial domain)
R = 2; % R^(2*S) total beams
MN = M*N;

%% Angle sampling for final stage
r_lin = linspace(0,R^S-1,R^S);
a1 = 1 - (2*r_lin + 1)/R^S; % sin(theta)cos(phi)
b1 = a1; % sin(theta)sin(phi)

a_idx = 10;
b_idx = 15;
phi = atan(b1(b_idx)/a1(a_idx));
theta = asin(sqrt(a1(a_idx)^2 + b1(b_idx)^2));

fc = 20e9;  % center frequency
c = physconst('LightSpeed');
W = 5e9; % bandwidth
fs = 2.5*W; % sampling frequency
Ws = 5e9; % to prevent spacial aliasing
T = 1/fs;
lambda = c/(fc+Ws); % wavelngth
m = (-sqrt(M)/2+1/2):(sqrt(M)/2-1/2);
[X_mesh, Y_mesh] = meshgrid(m,m);
x_pos = [X_mesh(:),Y_mesh(:)]*lambda/2; % sensor positions
u_s = sin(theta)*[cos(phi);sin(phi)]; % normal vector for determining delays
tau = x_pos*u_s/c; % relative delays to phase center
eta = lambda/(2*c);

%% Signal specs that do not need to be redifined in loop
n_sinusoids = 100; % number of sinusoids in signal
n = [0:N-1]/fs; % temporal sample vectors
t_array = n-tau; % delays across array and time
t = t_array(:); %
T_aperture = max(t)-min(t);% temporal aperture of the array

demod_phase = repmat(exp(-1i*2*pi*fc*tau),N,1);
mod_phase = exp(1i*2*pi*fc*tau);

taps = 8; % for conventional delay and sum
edge_lim = taps + ceil(fs*(max(tau)-min(tau))/2);% clip signal to account for edge distortion
comp_idxs = (edge_lim+1):(N-edge_lim-1); % comparison indicies

taps2 = 16; % for fast delay and sum
edge_lim2 = taps2 + ceil(fs*(max(tau)-min(tau))/2);% clip signal to account for edge distortion
comp_idxs2 = (edge_lim2+1):(N-edge_lim2-1); % comparison indicies

t_nyq = [0:N-1]/fs; % nyquist samples for testing

%% simulation specifications
trials = 2;
SNRs = -30:15:30;
SNR_ds = zeros(trials, length(SNRs));
SNR_fast_ds = zeros(trials, length(SNRs));

for ii=1:trials
    for jj=1:length(SNRs)

        SNR_dB = SNRs(jj); % signal to noise ratio
        a = (randn(n_sinusoids,1) + 1i*randn(n_sinusoids,1))/sqrt(2); % sinusoid amplitudes
        f = W*(2*rand(n_sinusoids,1)-1); % sinusoid frequencies
        s = sig_gen(a,f,n_sinusoids,t);
        sigma = (norm(s)/sqrt(MN))*10^(-SNR_dB/20);
        noise = sigma*(randn(MN,1)+1i*randn(MN,1))/sqrt(2);
        SNR_check = db(norm(s)/norm(noise));
        fprintf('Trial: %d, Set SNR: %.2f, Measured SNR: %.2f\n',ii,SNR_dB,SNR_check)

        y = s + noise;
        y_2d = reshape(y,[M,N]);
        y_3d = reshape(y, [sqrt(M),sqrt(M),N,1,1]); % sqrt(M)x sqrt(M) x N x R x R

        y_multibeamformed_first = zeros(sqrt(M),N,1,R^S);
        %% first FDS along b1 
        for m_a1 = 1:sqrt(M)
            v_now = reshape(y_3d(m_a1,:,:,1,1),[sqrt(M),N,1]);
            for stage=1:S
                v_prev = v_now;
                num_virtual_sensors = sqrt(M)/2^stage;
                num_sectors = R^stage;
                v_now = zeros(num_virtual_sensors,N,num_sectors); % init partial beamforming at stage s
    
                for mm=0:num_virtual_sensors-1
                    for rr = 0:num_sectors-1
                        mi = mm+1;
                        r = floor((rr)/R)+1;
                        delays = ((-1)^(rr))*(1/(R^stage))*eta*[2^(stage-1)*(2*mm+0.5)-sqrt(M)/2; 2^(stage-1)*(2*mm+1.5)-sqrt(M)/2];
                        y_interp = reshape(v_prev([2*mm+1,2*mm+2],:,r),[2,N]);
                        v_prev_interped = TTD_Beamformer(y_interp.',-delays, fs, @Psi_sinc, 'R',taps2).';
                        v_now(mi,:,rr+1) = sum(v_prev_interped,1)./2;
                    end
                end
            end
            y_multibeamformed_first(m_a1,:,1,:) = v_now;
        end

        y_multibeamformed_total = zeros(N,R^S,R^S);
        %% second FDS along a1
        for r_b1=1:R^S
            v_now = reshape(y_multibeamformed_first(:,:,1,r_b1),[sqrt(M),N,1]);
            for stage=1:S
                v_prev = v_now;
                num_virtual_sensors = sqrt(M)/2^stage;
                num_sectors = R^stage;
                v_now = zeros(num_virtual_sensors,N,num_sectors); % init partial beamforming at stage s
    
                for mm=0:num_virtual_sensors-1
                    for rr = 0:num_sectors-1
                        mi = mm+1;
                        r = floor((rr)/R)+1;
                        delays = ((-1)^(rr))*(1/(R^stage))*eta*[2^(stage-1)*(2*mm+0.5)-sqrt(M)/2; 2^(stage-1)*(2*mm+1.5)-sqrt(M)/2];
                        y_interp = reshape(v_prev([2*mm+1,2*mm+2],:,r),[2,N]);
                        v_prev_interped = TTD_Beamformer(y_interp.',-delays, fs, @Psi_sinc, 'R',taps2).';
                        v_now(mi,:,rr+1) = sum(v_prev_interped,1)./2;
                    end
                end
            end
            y_multibeamformed_total(:,:,r_b1) = v_now;
        end

        %--- adjust for positive a1 constraint
        if a_idx<=R^S/2 
            v_pred = y_multibeamformed_total(:,b_idx,a_idx);
        else
            v_pred = y_multibeamformed_total(:,R^S - b_idx+1,R^S - a_idx+1);
        end
        s_nyq_1 = sig_gen(a,f,n_sinusoids,t_nyq+1/fs);
        s_nyq_2 = sig_gen(a,f,n_sinusoids,t_nyq+2*S/fs);

        y_filter_ttd = TTD_Beamformer((y_2d).',-tau, fs, @Psi_sinc, 'R',taps).';
        y_bf_ttd = sum(y_filter_ttd(:,1:N),1)/M;

        %% offset 
        SNR_ds(ii,jj) = norm(s_nyq_1(comp_idxs))/norm(s_nyq_1(comp_idxs) - y_bf_ttd(comp_idxs));
        SNR_fast_ds(ii,jj) = norm(s_nyq_2(comp_idxs2))/norm(s_nyq_2(comp_idxs2) - v_pred(comp_idxs2).');

    end
end

figure(1)
plot(SNRs,SNRs + db(M)/2,'LineWidth',1)
hold on
grid on
plot(SNRs,db(mean(SNR_ds)),'-^','LineWidth',1)
hold on
grid on
plot(SNRs, db(mean(SNR_fast_ds)),'-o','LineWidth',1)
xlabel('Nominal SNR (dB)','Interpreter','latex')
ylabel('Beamformed SNR (dB)','Interpreter','latex')
legend({'Ideal', 'Delay and Sum R=16', 'Fast Delay and Sum (R=32)'},'Location','northwest','Interpreter','latex')
% exportgraphics(gcf, 'upa_snr_fast_ds.pdf', 'ContentType', 'vector')


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