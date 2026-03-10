clear
clc
close all

%% set random number generator seed
seed = 1;%randi([1,10000]);
rng(seed)

B_lst = 2^3:6:2^6;

fc = 20e9;  % center frequency
c = physconst('LightSpeed');
lambda = c/fc; % wavelngth
W = 5e9; % bandwidth
fs = 2.5*W; % sampling frequency
T = 1/fs;
trials = 20;

Runtime_fbst = zeros(trials, length(B_lst));
Runtime_fbst_superfast_toeplitz = zeros(trials, length(B_lst));

Runtime_ttd16 = zeros(trials, length(B_lst));
Runtime_ttd32 = zeros(trials, length(B_lst));

for b=1:length(B_lst)

    B = B_lst(b); % no. of beams
    M = 2^8; % array size
    N = 2^6; % no. of samples
    L = ceil(N); % no. of samples in frequency is 2*L
    MN = M*N;
    
    J = linspace(0,B-1,B);
    K = linspace(0,B-1,B);
    j_idx = randsample(1:B,1);
    k_idx = randsample(1:B,1);
    theta = asin(sqrt(J(j_idx)^2 + K(k_idx)^2)/(sqrt(2)*B));
    phi = asin(K(k_idx)/(sqrt(J(j_idx)^2 + K(k_idx)^2)));

    m = (-sqrt(M)/2+1/2):(sqrt(M)/2-1/2);
    [X_mesh, Y_mesh] = meshgrid(m,m);
    x_pos = [X_mesh(:),Y_mesh(:)]*lambda/2; % sensor positions
    u_s = sin(theta)*[cos(phi);sin(phi)]; % normal vector for determining delays
    
    %% Signal specs that do not need to be redifined in loop
    n_sinusoids = 100; % number of sinusoids in signal
    n = [0:N-1]/fs; % temporal sample vectors
    tau = x_pos*u_s/c; % relative delays to phase center
    t_array = n-tau; % delays across array and time
    t = t_array(:); %
    T_aperture = max(t)-min(t);% temporal aperture of the array
    demod_phase = repmat(exp(-1i*2*pi*fc*tau),N,1);
    mod_phase = exp(1i*2*pi*fc*tau);
    % demod_phase = repmat(demod_phase,N,1);
    R = 8;
    edge_lim = R + ceil(fs*(max(tau)-min(tau))/2);% clip signal to account for edge distortion
    comp_idxs = (edge_lim+1):(N-edge_lim-1); % comparison indicies

    R2 = 16;
    edge_lim2 = R2 + ceil(fs*(max(tau)-min(tau))/2);% clip signal to account for edge distortion
    comp_idxs2 = (edge_lim2+1):(N-edge_lim2-1); % comparison indicies
    
    %% construct the toeplitz matrix
    delta = 1e-5;
    A_eplitz = zeros(2*L,2*L);
    for ll=1:2*L
        for kk=1:2*L
            f_l = (1/(2*T*L))*(ll-1-L);
            f_k = (1/(2*T*L))*(kk-1-L);
            A_eplitz(ll,kk) = sum(exp(1i*2*pi*(f_k - f_l)*t));
        end
    end
    toeplitz_inv = (A_eplitz + delta*eye(2*L))^-1;

    %% get canonical vectors of Toeplitz inverse via brute force
    [X1,X2,X3,X4] = eval_canonical_vecs_brute_force(A_eplitz + delta*eye(2*L),2*L);

    %% equispaced sampled basis matrix
    psi = zeros(N,2*L);
    t_nyq = [0:N-1]'/fs; % nyquist samples for testing
    for nn=1:N
        for ll=1:2*L
            f_l = (1/(2*T*L))*(ll-1-L);
            psi(nn,ll) = exp(1i*2*pi*f_l*t_nyq(nn));
        end
    end

    %% interpolator for TTD beamformer
    [out] = Kernel_Filter_Gen(fs, @Psi_sinc,'R',R);
    [h, h_idx] = out.Filter_Gen(-tau);

    [out2] = Kernel_Filter_Gen(fs, @Psi_sinc,'R',R2);
    [h2, h_idx2] = out2.Filter_Gen(-tau);
    SNR = 10;
    for ii=1:trials
        fprintf('Trial: %d, Beam no.: %d\n',ii,b)
        SNR_dB = SNR; % signal to noise ratio
        a = (randn(n_sinusoids,1) + 1i*randn(n_sinusoids,1))/sqrt(2); % sinusoid amplitudes
        f = W*(2*rand(n_sinusoids,1)-1); % sinusoid frequencies
        s = sig_gen(a,f,n_sinusoids,t);
        sigma = (norm(s)/sqrt(MN))*10^(-SNR_dB/20);
        noise = sigma*(randn(MN,1)+1i*randn(MN,1))/sqrt(2);
        SNR_check = db(norm(s)/norm(noise));

        s_demod = s.*demod_phase;

        y = s_demod + noise;
        y_2d = reshape(y, [M,N]);
        y_3d = reshape(y, [sqrt(M),sqrt(M),N]);

        %% runtime for fbst (precompute inverse)
        ts_fbst = tic;
        F_transform = zeros(sqrt(M),sqrt(M),2*L); % transform w.r.t. temporal frequency
        A = -1;%exp(-1i*2*pi*T*W);
        W_czt = exp(-1i*pi/L);
        for mm1=1:sqrt(M)
            for mm2=1:sqrt(M)
                F_transform(mm1,mm2,:) = czt(y_3d(mm1,mm2,:),2*L,W_czt,A);
            end
        end

        X2_transform = zeros(sqrt(M),B,2*L); % transform w.r.t. second spatial frequency
        A= 1;
        B_list = reshape(linspace(0,B-1,B), [1,B]);
        for ll=1:2*L
            for mm1=1:sqrt(M)
                l_dash = ll-1;
                % A = exp(-1i*pi*(1/(sqrt(2)*B))*(1 + (1/(2*fc*T*L))*(l_dash-L)));%exp(-1i*pi*(1 + 2*(W/fc)*(l_dash-L/2)/L));
                W_czt = exp(1i*pi*(1/(sqrt(2)*B))*(1 + (1/(2*fc*T*L))*(l_dash-L)));
                coeff = exp(-1i*pi*(1/(sqrt(2)*B))*(1 + (1/(2*fc*T*L))*(l_dash-L))*(sqrt(M)/2-1/2)*B_list);
                X2_transform(mm1,:,ll) = czt(F_transform(mm1,:,ll),B,W_czt,A);
                X2_transform(mm1,:,ll) = X2_transform(mm1,:,ll).*coeff;
            end
        end

        X1_transform = zeros(B,B,2*L); % final transform matrix (after taking czt along the first spatial frequency)
        for ll=1:2*L
            for b2=1:B
                l_dash = ll-1;
                A = 1;%exp(-1i*pi*(1/(sqrt(2)*B))*(1 + (1/(2*fc*T*L))*(l_dash-L))); % check
                W_czt = exp(1i*pi*(1/(sqrt(2)*B))*(1 + (1/(2*fc*T*L))*(l_dash-L)));
                coeff = exp(-1i*pi*(1/(sqrt(2)*B))*(1 + (1/(2*fc*T*L))*(l_dash-L))*(sqrt(M)/2-1/2)*B_list);
                X1_transform(:,b2,ll) = czt(X2_transform(:,b2,ll),B,W_czt,A);
                X1_transform(:,b2,ll) = X1_transform(:,b2,ll).*coeff.';
            end
        end

        w_l = reshape(X1_transform(k_idx,j_idx,:),[2*L,1]);
        t_inv = tic;
        for bb=1:B^2
            alpha = toeplitz_inv*w_l;
            y_nyq = psi*alpha;
        end
        inv_time = toc(t_inv);
        runtime_fbst = toc(ts_fbst);

        %% superfast toeplitz 
        t_inv_superfast = tic;
        for bb=1:B^2
            b1 = toep_mult(X2,w_l);
            b2 = toep_mult(X4,w_l);
            b3 = toep_mult(X1,b1);
            b4 = toep_mult(X3,b2);
            alpha_superfast = b3 - b4;
            y_nyq_superfast = psi*alpha_superfast;
        end
        inv_time_superfast = toc(t_inv_superfast);

        %% runtime for ttd with 16 taps
        ts2 = tic;
        y_ttd_in = (y_2d.*mod_phase).';
        N_filter = N+2*R+1;
        
        % [out] = Kernel_Filter_Gen(fs, Psi,'R',R);
        % [h, h_idx] = out.Filter_Gen(delays);
        h = squeeze(h);
        h_idx = squeeze(h_idx);
        
        wrapN = @(x, n) (1 + mod(x-1, n));
        
        H = zeros(N_filter,M);
        
        for i=1:M
            H(wrapN(h_idx(i,:),N_filter),i) = h(i,:);
        end
        
        Hf = fft(H,[],1);
        for bb=1:B^2
            % y_filter_ttd = TTD_Beamformer((y.*mod_phase).',-tau, fs, @Psi_sinc,h, h_idx, 'R',R).';
            Xf = fft(y_ttd_in,N_filter,1);
            Y = ifft(Hf.*Xf,[],1);
            y_filter_ttd = Y(1:N,:).';
            y_bf_ttd = reshape(sum(y_filter_ttd(:,1:N),1)/M, [N,1]);
        end
        runtime_ttd = toc(ts2);

        %% runtime for ttd with 32 taps
        ts3 = tic;
        y_ttd_in2 = (y_2d.*mod_phase).';
        N_filter2 = N+2*R2+1;
        
        % [out] = Kernel_Filter_Gen(fs, Psi,'R',R);
        % [h, h_idx] = out.Filter_Gen(delays);
        h2 = squeeze(h2);
        h_idx2 = squeeze(h_idx2);
        
        wrapN2 = @(x, n) (1 + mod(x-1, n));
        
        H2 = zeros(N_filter2,M);
        
        for i=1:M
            H2(wrapN2(h_idx2(i,:),N_filter2),i) = h2(i,:);
        end
        
        Hf2 = fft(H2,[],1);
        for bb=1:B^2
            % y_filter_ttd = TTD_Beamformer((y.*mod_phase).',-tau, fs, @Psi_sinc,h, h_idx, 'R',R).';
            Xf2 = fft(y_ttd_in2,N_filter2,1);
            Y2 = ifft(Hf2.*Xf2,[],1);
            y_filter_ttd2 = Y2(1:N,:).';
            y_bf_ttd2 = reshape(sum(y_filter_ttd2(:,1:N),1)/M, [N,1]);
        end
        runtime_ttd16 = toc(ts3);

        Runtime_fbst(ii,b) = runtime_fbst;
        Runtime_fbst_superfast_toeplitz(ii,b) =runtime_fbst - inv_time + inv_time_superfast;
        Runtime_ttd16(ii,b) = runtime_ttd;
        Runtime_ttd32(ii,b) = runtime_ttd16;

    end

end

figure(1)
plot(B_lst.^2,mean(Runtime_fbst),'-^','LineWidth',1)
hold on
grid on
plot(B_lst.^2, mean(Runtime_fbst_superfast_toeplitz),'-^','LineWidth',1);
hold on
grid on
plot(B_lst.^2, mean(Runtime_ttd16),'-*','LineWidth',1)
hold on
grid on
plot(B_lst.^2, mean(Runtime_ttd32),'-*','LineWidth',1)
xlabel('no. of beams (B)','Interpreter','latex')
ylabel('runtime (s)','Interpreter','latex')
legend({'FBST(Inverse pre-computed)','FBST (superfast toeplitz)','Delay and Sum (R=16)','Delay and Sum (R=32)'},'Location','northwest','Interpreter','latex')


%% supporting functions
function [X] = sig_gen(a,f,n_sinusoids,t)
X = zeros(size(t));
for ii = 1:n_sinusoids % For each antenna
    sig = (a(ii)*exp(1i*2*pi*(f(ii))*(t)));
    X = X + sig;
end
end

function y = Psi_sinc(z)
y = sinc(z);
end