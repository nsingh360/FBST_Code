clear
clc
close all

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

B = 2^7; % no. of beams
M = 2^7; % array size
N_list = [2^4:8:2^7]; % no. of samples

Theta = zeros(2*B+1,1);
idx = B;%randi(B,1);
for i=1:2*B+1
    Theta(i) = asin((i-1-B)/B);
end

%% signal and specs
n_sinusoids = 100; % number of sinusoids in signal
fc = 20e9;  % center frequency
Ws = 5e9; % to prevent spatial aliasing
c = physconst('LightSpeed');
lambda = c/(fc+Ws); % wavelngth
W = 5e9; % bandwidth
fs = 2.5*W; % sampling frequency
T = 1/fs;
theta = Theta(idx);
m = [(-M/2+1/2):(M/2-1/2)]';
x_pos = m*lambda/2; % sensor positions
u_s = sin(theta); % normal vector for determining delays
tau = x_pos*u_s/c;

%% interferer specs
interference_struct.list = [1,3];

totalNumbers = length(interference_struct.list) * 5; % Total unique numbers needed
randomNumbers = randperm(2*B-9, totalNumbers) + 4; % removing edge indices to account for interpolation
interference_struct.idx = reshape(randomNumbers, 1, length(interference_struct.list), 5); % Reshape into matrix
interference_struct.theta = asin(0.5*sin(Theta(interference_struct.idx)) + 0.5*sin(Theta(interference_struct.idx+1)));
interference_struct.u_i = sin(interference_struct.theta);
interference_struct.tau_i  = x_pos.*interference_struct.u_i/c;

trials = 5;
Runtime_beam_space = zeros(trials, length(N_list),length(interference_struct.list));
Runtime_array_space = zeros(trials, length(N_list), length(interference_struct.list));
SNR_dB = 30;
SIR_dB = -10;

for ii=1:trials
    for jj=1:length(N_list)
        N = N_list(jj);
        L = ceil(N);
        MN = M*N;
        sigma = 10^(-SNR_dB/20);
        % fprintf('Trial: %d, Snapshots: %.2f\n',ii,N);

        %% snapshot specs
        n = [0:N-1]/fs; % temporal sample vectors
        t_array = n-tau; % delays across array and time
        t = t_array(:); %
        T_aperture = max(t)-min(t);% temporal aperture of the array
        demod_phase = repmat(exp(-1i*2*pi*fc*tau),N,1);

        frequency_samples = (1/(2*T*L))*linspace(-L,L-1,2*L);
        F_s = exp(1i*2*pi*t*frequency_samples);
        F_s = F_s.*demod_phase;
        
        %% generate random signals
        a = (randn(n_sinusoids,1) + 1i*randn(n_sinusoids,1))/sqrt(2); % sinusoid amplitudes
        f = W*(2*rand(n_sinusoids,1)-1); % sinusoid frequencies
        s = sig_gen(a,f,n_sinusoids,t);
        s_scale = 1/rms(s);
        s_demod = s_scale*s.*demod_phase;

        %% equispaced sampled basis matrix
        psi = zeros(N,2*L);
        t_nyq = [0:N-1]/fs; % nyquist samples for testing
        for ll=1:2*L
            f_l = (1/(2*T*L))*(ll-1-L);
            psi(:,ll) = exp(1i*2*pi*f_l*t_nyq);
        end

        s_nyq = sig_gen(a,f,n_sinusoids,t_nyq')*s_scale;

        %% interference time generation
        for kk=1:length(interference_struct.list)
            num_interferers = interference_struct.list(kk);
            fprintf('Trial: %d, Snapshots: %.2f, interferers: %d\n',ii,N, num_interferers);
            t_interferers = zeros(MN,num_interferers);
            interferer_aperture = zeros(num_interferers,1);
            interferer_demod_phases = zeros(MN,num_interferers);
            interferer_fourier_matrix = zeros(MN,2*L,num_interferers);

            s_interference = zeros(MN,1);
            
            for ll=1:num_interferers
                temp_t = n - interference_struct.tau_i(:,kk,ll);
                t_interferers(:,ll) = temp_t(:);
                interferer_aperture(ll) = max(t_interferers(:,ll)) - min(t_interferers(:,ll));
                interferer_demod_phases(:,ll) = repmat(exp(-1i*2*pi*fc*interference_struct.tau_i(:,kk,ll)),N,1);

                a_i = (randn(n_sinusoids,1) + 1i*randn(n_sinusoids,1))/sqrt(2); % sinusoid amplitudes
                f_i = W*(2*rand(n_sinusoids,1)-1); % sinusoid frequencies
                s_i = sig_gen(a_i,f_i,n_sinusoids,temp_t(:));
                s_scale_i = 1/rms(s_i);
                s_demod_i = s_scale_i*s_i*10^(-SIR_dB/20).*interferer_demod_phases(:,ll);
                s_interference = s_interference + s_demod_i;

                F_i = exp(1i*2*pi*temp_t(:)*frequency_samples);
                interferer_fourier_matrix(:,:,ll) = F_i.*interferer_demod_phases(:,ll);
            end

            noise = sigma*(randn(MN,1)+1i*randn(MN,1))/sqrt(2);
            y = s_demod + s_interference + noise;

            y_array = reshape(y, [M,N]);
            %% chirp-z transform
            F_transform = zeros(M,2*L);
            A = -1;%exp(-1i*2*pi*T*W);
            W_czt = exp(-1i*pi/L);
            for mm=1:M
                F_transform(mm,:) = czt(y_array(mm,:),2*L,W_czt,A);
            end
            X_czt = zeros(2*L,2*B+1);
            A= 1;
            B_list = reshape(linspace(0,2*B,2*B+1), [1,2*B+1]);
            for ll=1:2*L
                l_dash = ll-1;
                A = exp(1i*pi*(1/B)*(fc/(fc+Ws) + (1/(2*(fc+Ws)*T*L))*(l_dash-L))*B);
                W_czt = exp(1i*pi*(1/B)*(fc/(fc+Ws) + (1/(2*(fc+Ws)*T*L))*(l_dash-L)));
                coeff1 = exp(-1i*pi*(1/B)*(fc/(fc+Ws) + (1/(2*(fc+Ws)*T*L))*(l_dash-L))*(M/2-1/2)*B_list); % adjusting for array center
                coeff2 = exp(1i*pi*(1/B)*(fc/(fc+Ws) + (1/(2*(fc+Ws)*T*L))*(l_dash-L))*(M/2-1/2)*B); % adjusting for negative angles
                X_czt(ll,:) = czt(F_transform(:,ll),2*B+1,W_czt,A);
                X_czt(ll,:) = coeff2*X_czt(ll,:).*coeff1;
            end

            %% measure array space runtime
            ts_array_space = tic;
            y_tilde = y;
            for ll=1:num_interferers
                y_tilde = y_tilde - interferer_fourier_matrix(:,:,ll)*inv(interferer_fourier_matrix(:,:,ll)'*interferer_fourier_matrix(:,:,ll) + 1e-8*eye(2*L))*interferer_fourier_matrix(:,:,ll)'*y;
            end
            w_array_space = F_s'*y_tilde;
            runtime_array_space = toc(ts_array_space);

            %% measure beam space runtime
            ts_beam_space = tic;
            w = X_czt(:,idx);
            for ll=1:num_interferers
                interp_idx = interference_struct.idx(1,kk,ll);
                w_l_interp = interp1(Theta(interp_idx-3:interp_idx+3),X_czt(:,interp_idx-3:interp_idx+3).',interference_struct.theta(1,kk,ll),'spline').';
                w = w - F_s'*interferer_fourier_matrix(:,:,ll)*inv(interferer_fourier_matrix(:,:,ll)'*interferer_fourier_matrix(:,:,ll) + 1e-8*eye(2*L))*w_l_interp;
            end
            runtime_beam_space = toc(ts_beam_space);

            Runtime_array_space(ii,jj,kk) = runtime_array_space;
            Runtime_beam_space(ii,jj,kk) = runtime_beam_space;
        end

    end
end

%% plot results

figure(1)
legendi = {};
symb = {'--' , '-'};
count = 1;
for i=1:length(interference_struct.list)
    plot(N_list, mean(Runtime_array_space(:,:,i)),symb{i},'LineWidth',1,'Color','#0072BD');
    hold on
    plot(N_list,mean(Runtime_beam_space(:,:,i)),symb{i},'LineWidth',1,'Color','#D95319');
    hold on 
    legendi{count} = sprintf('Array space K = %d',interference_struct.list(i));
    legendi{count+1} = sprintf('Beam space K = %d',interference_struct.list(i));
    count = count+2;
end

xlim([16 128])
set(gca, 'XScale', 'log'); % Set x-axis to logarithmic scale
xticks([16 32 64 128]); % Set specific log2 tick marks
xticklabels({'2^4', '2^5', '2^6', '2^7'}); % Label as powers of 2

xlabel('Snapshots (N)','Interpreter','latex')
ylabel('Runtime (seconds)','Interpreter','latex')
legend(legendi,'Location','northwest','Interpreter','latex','FontSize',12)
grid on
% exportgraphics(gcf, 'interference_runtime_vs_snapshots.pdf', 'ContentType', 'vector');


%% supporting functions
function [X] = sig_gen(a,f,n_sinusoids,t)
X = zeros(size(t));
for ii = 1:n_sinusoids % For each antenna
    sig = (a(ii)*exp(1i*2*pi*(f(ii))*(t)));
    X = X + sig;
end
end