clear
clc
close all

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

S = 7; % stages in the beamformer
M = 2^S; % array size
R = 2; % 2-radix beamformer
L = 2; % spatial downsampling factor
N = 2^6; % no. of samples
MN = M*N;


%% Angle sampling for final stage
r_lin = linspace(0,R^S-1,R^S);
Theta = asin(1 - (2*r_lin + 1)/R^S);
idx = 50;

fc = 20e9;  % center frequency
Ws = 5e9; % to prevent spacial aliasing
c = physconst('LightSpeed');
lambda = c/(fc+Ws); % wavelngth
W = 5e9; % bandwidth
fs = 2.5*W; % sampling frequency
T = 1/fs;
theta = Theta(idx); % select an on-grid sampling angle for testing
m = [(-M/2+1/2):(M/2-1/2)]';
x_pos = m*lambda/2; % sensor positions
u_s = sin(theta); % normal vector for determining delays
eta = lambda/(2*c);

%% Signal specs that do not need to be redifined in loop
n_sinusoids = 100; % number of sinusoids in signal
n = [0:N-1]/fs; % temporal sample vectors
tau = x_pos*u_s/c; % relative delays to phase center
t_array = n-tau; % delays across array and time
t = t_array(:); %
T_aperture = max(t)-min(t);% temporal aperture of the array
demod_phase = repmat(exp(-1i*2*pi*(fc)*tau),N,1);
mod_phase = exp(1i*2*pi*(fc)*tau);
% demod_phase = repmat(demod_phase,N,1);
t_nyq = [0:N-1]/fs;

taps = 8;
edge_lim = taps + ceil(fs*(max(tau)-min(tau))/2);% clip signal to account for edge distortion
comp_idxs = (edge_lim+1):(N-edge_lim-1); % comparison indicies

taps2 = 16;
edge_lim2 = taps2 + ceil(fs*(max(tau)-min(tau))/2);% clip signal to account for edge distortion
comp_idxs2 = (edge_lim2+1):(N-edge_lim2-1); % comparison indicies

%% simulation specifications
trials = 2;
SNRs = -30:6:30;
SNR_ds = zeros(trials, length(SNRs));
SNR_fast_ds = zeros(trials,length(SNRs));

% %% final stage fast ds
% delays2 = zeros(M,1);
% for ll=1:2^S
%     l = ll-1;
%     r = 50-1;
% 
%     delay = 0;
% 
%     for k=0:S-1
%         delay = delay + eta*(2^(k)*(floor(l/(2^(k))) + 0.5) - M/2)*(2^(-k-1))*((-1)^(floor((2^(-S+1+k))*r)));
%     end
%     delays2(ll) = delay;
% end

for ii = 1:trials
    for jj = 1:length(SNRs)

        SNR_dB = SNRs(jj); % signal to noise ratio
        a = (randn(n_sinusoids,1) + 1i*randn(n_sinusoids,1))/sqrt(2); % sinusoid amplitudes
        f = W*(2*rand(n_sinusoids,1)-1); % sinusoid frequencies
        s = sig_gen(a,f,n_sinusoids,t);
        sigma = (norm(s)/sqrt(MN))*10^(-SNR_dB/20);
        noise = sigma*(randn(MN,1)+1i*randn(MN,1))/sqrt(2);
        SNR_check = db(norm(s)/norm(noise));
        fprintf('Trial: %d, Set SNR: %.2f, Measured SNR: %.2f\n',ii,SNR_dB,SNR_check)

        y = s + noise;
        y = reshape(y, [M,N]); % stage 0 init

        v_now = reshape(y, [M,N,1]);
        fs_ds1 = tic;
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
                    v_prev_interped = TTD_Beamformer(y_interp.',-delays, fs, @Psi_sinc, 'R',taps2).';
                    v_now(m,:,rr+1) = sum(v_prev_interped,1)./2;
                end
            end
        end
        fs_time = toc(fs_ds1);

        ds1 = tic;
        y_filter_ttd = TTD_Beamformer(y.',-tau, fs, @Psi_sinc, 'R',taps).';
        y_bf_ttd = sum(y_filter_ttd(:,1:N),1)/M;
        ds_time = toc(ds1);

        s_nyq_1 = sig_gen(a,f,n_sinusoids,t_nyq + 1/fs);
        s_nyq_2 = sig_gen(a,f,n_sinusoids,t_nyq + (S)/fs);
        v_pred = v_now(1,:,idx);

        % %% final stage fast ds
        % y_filter_fast_ttd2 = TTD_Beamformer(y.',-delays2, fs, @Psi_sinc, 'R',taps2).';
        % y_bf_fast_ttd2 = sum(y_filter_fast_ttd2(:,1:N),1)/M;

        %% offset 
        SNR_ds(ii,jj) = norm(s_nyq_1(comp_idxs))/norm(s_nyq_1(comp_idxs) - y_bf_ttd(comp_idxs));
        SNR_fast_ds(ii,jj) = norm(s_nyq_2(comp_idxs2))/norm(s_nyq_2(comp_idxs2) - v_pred(comp_idxs2));
        
    end
end


figure(1)
plot(SNRs,SNRs + db(M)/2,'LineWidth',1)
hold on
grid on
plot(SNRs,db(mean(SNR_ds)),'-^','LineWidth',1)
hold on
grid on
plot(SNRs, db(mean(SNR_fast_ds)),'-*','LineWidth',1)
xlabel('Nominal SNR (dB)','Interpreter','latex','FontSize',12)
ylabel('Beamformed SNR (dB)','Interpreter','latex','FontSize',12)
legend({'Ideal','Delay and Sum R=16', 'Fast Delay and Sum R=32'},'Location','northwest','Interpreter','latex','FontSize',12)
% exportgraphics(gcf, 'ula_snr_fast_ds.pdf', 'ContentType', 'vector');


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