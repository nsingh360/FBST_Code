clear
clc
close all

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

M_lst = [2^1 2^2 2^3 2^4 2^5 2^6 2^7 2^8 2^9 2^10 2^11 2^12];
trials = 20;

Runtime_fbst_superfast = zeros(trials, length(M_lst));
Runtime_fbst_precompute = zeros(trials, length(M_lst));
Runtime_ds_16_tap = zeros(trials, length(M_lst));
Runtime_ds_32_tap = zeros(trials, length(M_lst));
Runtime_fast_ds_16_tap = zeros(trials, length(M_lst));
Runtime_fast_ds_32_tap = zeros(trials, length(M_lst));

Runtime_chirp_z = zeros(trials, length(M_lst));

for mm1 = 1:length(M_lst)
    M = M_lst(mm1);
    N = 2^6;
    B = 2^(log2(M)-1); % create a total of M+1 fbst beams
    R = 2^(log2(M)); % beams for ds and fast-ds
    L = ceil(N); % no. of samples in frequency is 2*L
    MN = M*N;

    Theta = zeros(2*B+1,1);
    idx = 2*B;%randi(B,1);
    for i=1:2*B+1
        Theta(i) = asin((i-1-B)/B);
    end

    fc = 20e9;  % center frequency
    Ws = 5e9; % to prevent spacial aliasing
    fc_tilde = (fc+Ws); % adjustment for spatial aliasing
    fc_norm = fc/fc_tilde;
    c = physconst('LightSpeed');
    lambda = c/(fc+Ws); % wavelngth
    W = 5e9; % bandwidth
    fs = 2.5*W; % sampling frequency
    T = 1/fs;
    theta = Theta(idx);
    m = [(-M/2+1/2):(M/2-1/2)]';
    x_pos = m*lambda/2; % sensor positions
    u_s = sin(theta); % normal vector for determining delays

    %% parameters for fast delay and sum
    S = log2(M);
    R_fac = 2; % 2-radix fast delay and sum
    L_fac = 2; % downsampling factor in spatial dimension
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
    taps = 8;
    edge_lim = taps + ceil(fs*(max(tau)-min(tau))/2);% clip signal to account for edge distortion
    comp_idxs = (edge_lim+1):(N-edge_lim-1); % comparison indicies
    
    taps2 = 16;
    edge_lim2 = taps2 + ceil(fs*(max(tau)-min(tau))/2);% clip signal to account for edge distortion
    comp_idxs2 = (edge_lim2+1):(N-edge_lim2-1); % comparison indicies

    taps3 = 128;

    %% construct the toeplitz matrix
    delta = 1e-5;
    freq_samples = 1/(2*T*L)*linspace(-L,L-1,2*L);
    F_basis = exp(1i*2*pi*t*freq_samples).*demod_phase;
    A_eplitz = F_basis'*F_basis;
    toeplitz_inv = inv(A_eplitz + (delta)*eye(2*L));

    %% get canonical vectors of Toeplitz inverse via brute force
    [X1,X2,X3,X4] = eval_canonical_vecs_brute_force(A_eplitz + delta*eye(2*L),2*L);

    %% pre-processing for superfast toeplitz
    c_x1 = X1(:, 1);
    r_x1 = X1(1, :);
    p_x1 = [c_x1; 0; flipud(conj(r_x1(2:end)'))];
    pf_x1 = fft(p_x1);

    c_x2 = X2(:, 1);
    r_x2 = X2(1, :);
    p_x2 = [c_x2; 0; flipud(conj(r_x2(2:end)'))];
    pf_x2 = fft(p_x2);

    c_x3 = X3(:, 1);
    r_x3 = X3(1, :);
    p_x3 = [c_x3; 0; flipud(conj(r_x3(2:end)'))];
    pf_x3 = fft(p_x3);

    c_x4 = X4(:, 1);
    r_x4 = X4(1, :);
    p_x4 = [c_x4; 0; flipud(conj(r_x4(2:end)'))];
    pf_x4 = fft(p_x4);

    %% equispaced sampled basis matrix
    t_nyq = [0:N-1]/fs; % nyquist samples for testing
    psi = exp(1i*2*pi*t_nyq'*freq_samples);

    %% interpolator for TTD beamformer
    [out] = Kernel_Filter_Gen(fs, @Psi_sinc,'R',taps);
    [h, h_idx] = out.Filter_Gen(-tau);

    [out2] = Kernel_Filter_Gen(fs, @Psi_sinc,'R',taps2);
    [h2, h_idx2] = out2.Filter_Gen(-tau);

    [out3] = Kernel_Filter_Gen(fs, @Psi_sinc,'R',taps3);

    %% Filters for chirp-z
    %----- stage-1 filters
    A_stage_1 = -1;%exp(-1i*2*pi*T*W);
    W_czt_stage_1 = exp(-1i*pi/L);
    %------ filter design for chirp-z stage-1 (taken for each row of y)
    m_len1 = N; % input signal length
    k1_1 = 2*L; % chirp-z length
    n_len = 1;
    nfft_1 = 2^nextpow2(m_len1+k1_1-1);
    %------- Premultiply data.
    kk = ((-m_len1+1):max(k1_1-1,m_len1-1)).';
    kk2 = (kk .^ 2) ./ 2;
    ww_stage_1 = W_czt_stage_1 .^ (kk2);   % <----- Chirp filter is 1./ww
    nn = (0:(m_len1-1))';
    aa = A_stage_1 .^ ( -nn );
    aa_stage_1 = aa.*ww_stage_1(m_len1+nn);
    fv_stage_1 = fft( 1 ./ ww_stage_1(1:(k1_1-1+m_len1)), nfft_1 );   % <----- Chirp filter.

    %------ stage-2 filters
    m_len2 = M; % input signal length
    k1_2 = 2*B+1; % chirp-z length
    n_len = 1;
    nfft_2 = 2^nextpow2(m_len2+k1_2-1);
    %------- Premultiply data.
    kk = ((-m_len2+1):max(k1_2-1,m_len2-1)).';
    kk2 = (kk .^ 2) ./ 2;
    nn = (0:(m_len2-1))';

    A_coefs_stage_2 = zeros(2*L,1);
    W_czt_coefs_stage_2 = zeros(2*L,1);
    ww_stage2 = zeros(length(kk2),2*L);
    aa_stage2 = zeros(m_len2,2*L);
    FV_filters_stage_2 = zeros(nfft_2,2*L);
    
    for ll=1:2*L
        l_dash = ll-1;
        A_coefs_stage_2(ll,1) = exp(1i*pi*(1/B)*(fc_norm + (1/(2*fc_tilde*T*L))*(l_dash-L))*B);
        aa = A_coefs_stage_2(ll,1) .^ ( -nn );
        W_czt_coefs_stage_2(ll,1) = exp(1i*pi*(1/B)*(fc_norm + (1/(2*fc_tilde*T*L))*(l_dash-L)));
        ww_stage2(:,ll) = W_czt_coefs_stage_2(ll,1) .^ (kk2);   % <----- Chirp filter is 1./ww
        aa_stage2(:,ll) = aa.*ww_stage2(m_len2+nn,ll);
        FV_filters_stage_2(:,ll) = fft( 1 ./ ww_stage2(1:(k1_2-1+m_len2)), nfft_2 );   % <----- Chirp filter.
    end

    SNR = 10;
    for jj = 1:trials

        SNR_dB = SNR; % signal to noise ratio
        a = (randn(n_sinusoids,1) + 1i*randn(n_sinusoids,1))/sqrt(2); % sinusoid amplitudes
        f = W*(2*rand(n_sinusoids,1)-1); % sinusoid frequencies
        s = sig_gen(a,f,n_sinusoids,t);
        sigma = (norm(s)/sqrt(MN))*10^(-SNR_dB/20);
        noise = sigma*(randn(MN,1)+1i*randn(MN,1))/sqrt(2);
        SNR_check = db(norm(s)/norm(noise));
        fprintf('Trial: %d, Set SNR: %.2f, Array size: %.2f\n',jj,SNR_dB,M)

        s_demod = s.*demod_phase;
        y = s_demod + noise;
        y = reshape(y, [M,N]);

        %% chirp-z transform
        t_init_chirpz = tic;
        F_transform = zeros(M,2*L);
        % A = -1;%exp(-1i*2*pi*T*W);
        % W_czt = exp(-1i*pi/L);
        for mm=1:M
            y_in = y(mm,:).';
            y_in = y_in .* aa_stage_1(:,ones(1,n_len));
            fy = fft(  y_in, nfft_1 );
            fy = fy .* fv_stage_1(:,ones(1, n_len));
            g  = ifft( fy );
            g = g( m_len1:(m_len1+k1_1-1), :, :) .* ww_stage_1( m_len1:(m_len1+k1_1-1),ones(1, n_len) );
            F_transform(mm,:) = g;
        end

        X_czt = zeros(2*L,2*B+1);
        A= 1;
        B_list = reshape(linspace(0,2*B,2*B+1), [1,2*B+1]);
        for ll=1:2*L
            l_dash = ll-1;

            y_in = F_transform(:,ll);
            ww = W_czt_coefs_stage_2(ll,1) .^ (kk2);
            y_in = y_in .* aa_stage2(:,ll);
            fy = fft(  y_in, nfft_2 );
            fv = FV_filters_stage_2(:,ll);
            fy = fy .* fv(:,ones(1, n_len));
            g  = ifft( fy );
            g = g( m_len2:(m_len2+k1_2-1), :, :) .* ww_stage2( m_len2:(m_len2+k1_2-1),ll );
            % A = exp(1i*pi*(1/B)*(fc_norm + (1/(2*fc_tilde*T*L))*(l_dash-L))*B);
            % W_czt = exp(1i*pi*(1/B)*(fc_norm + (1/(2*fc_tilde*T*L))*(l_dash-L)));
            coeff1 = exp(-1i*pi*(1/B)*(fc_norm + (1/(2*fc_tilde*T*L))*(l_dash-L))*(M/2-1/2)*B_list); % adjusting for array center
            coeff2 = exp(1i*pi*(1/B)*(fc_norm + (1/(2*fc_tilde*T*L))*(l_dash-L))*(M/2-1/2)*B); % adjusting for negative angles
            X_czt(ll,:) = g;
            X_czt(ll,:) = coeff2*X_czt(ll,:).*coeff1;
        end
        w_l = X_czt(:,idx);
        runtime_chirpz = toc(t_init_chirpz);
        
        %% beamforming inverse problem (precompute inverse)
        t_init_precompute_fbst = tic;
        alpha = toeplitz_inv*w_l;
        y_nyq = psi*alpha;
        runtime_fbst_precompute_beam = toc(t_init_precompute_fbst);

        %% beamforming inverse problem (superfast toeplitz)
        t_init_superfast_fbst = tic;
        b1 = ifft(pf_x2 .* fft([w_l; zeros(length(w_l), 1)]));
        b1 = b1(1:length(w_l));
        % b1 = toep_mult(X2,w_l);
        b2 = ifft(pf_x4 .* fft([w_l; zeros(length(w_l), 1)]));
        b2 = b2(1:length(w_l));
        % b2 = toep_mult(X4,w_l);
        b3 = ifft(pf_x1 .* fft([b1; zeros(length(b1), 1)]));
        b3 = b3(1:length(b1));
        % b3 = toep_mult(X1,b1);
        b4 = ifft(pf_x3 .* fft([b2; zeros(length(b2), 1)]));
        b4 = b4(1:length(b2));
        % b4 = toep_mult(X3,b2);
        alpha_superfast = b3 - b4;
        y_nyq_superfast = psi*alpha_superfast;
        runtime_fbst_superfast_beam = toc(t_init_superfast_fbst);

        Runtime_fbst_precompute(jj,mm1) = runtime_chirpz + (2*B+1)*runtime_fbst_precompute_beam;
        Runtime_fbst_superfast(jj,mm1) = runtime_chirpz + (2*B+1)*runtime_fbst_superfast_beam;
        Runtime_chirp_z(jj,mm1) = runtime_chirpz;

        %% fast delay and sum beamforming 16 taps
        v_now = reshape(y.*mod_phase, [M,N,1]);
        count_filter_time = 0; % counter for subtracting filter initialization

        init_fast_ds = tic;
        for stage=1:S
            v_prev = v_now;
            num_virtual_sensors = M/2^stage;
            num_sectors = R_fac^stage;
            v_now = zeros(num_virtual_sensors,N,num_sectors); % init partial beamforming at stage s

            for mm2=0:num_virtual_sensors-1
                for rr = 0:num_sectors-1
                    m2 = mm2+1;
                    r = floor((rr)/R_fac)+1;
                    delays = ((-1)^(rr))*(1/(R_fac^stage))*eta*[2^(stage-1)*(2*mm2+0.5)-M/2; 2^(stage-1)*(2*mm2+1.5)-M/2];
                    y_interp = reshape(v_prev([2*mm2+1,2*mm2+2],:,r),[2,N]);
                    % v_prev_interped = TTD_Beamformer(y_interp.',-delays, fs, @Psi_sinc, 'R',taps2).';
                    % each stage filter design 
                    init_filter_fast_ds = tic;
                    [h_stage, h_stage_idx] = out.Filter_Gen(-delays);
                    y_stage = (y_interp).';
                    N_filter = N+2*taps+1;
                    h_stage = squeeze(h_stage);
                    h_stage_idx = squeeze(h_stage_idx);
                    
                    wrapN = @(x, n) (1 + mod(x-1, n));
                    
                    H_stage = zeros(N_filter,L_fac);
                    
                    for i=1:L_fac
                        H_stage(wrapN(h_stage_idx(i,:),N_filter),i) = h_stage(i,:);
                    end
                    Hf_stage = fft(H_stage,[],1);
                    count_filter_time = count_filter_time + toc(init_filter_fast_ds);

                    Xf_stage = fft(y_stage,N_filter,1);
                    Y_stage = ifft(Hf_stage.*Xf_stage,[],1);
                    y_filter_stage = Y_stage(1:N,:).';

                    v_now(m2,:,rr+1) = sum(y_filter_stage,1)./2;
                end
            end
        end
        runtime_fast_ds_16_tap = toc(init_fast_ds);
        Runtime_fast_ds_16_tap(jj,mm1) = runtime_fast_ds_16_tap - count_filter_time;

        %% fast delay and sum beamforming 32 taps
        v_now = reshape(y.*mod_phase, [M,N,1]);
        count_filter_time2 = 0; % counter for subtracting filter initialization

        init_fast_ds2 = tic;
        for stage=1:S
            v_prev = v_now;
            num_virtual_sensors = M/2^stage;
            num_sectors = R_fac^stage;
            v_now = zeros(num_virtual_sensors,N,num_sectors); % init partial beamforming at stage s

            for mm2=0:num_virtual_sensors-1
                for rr = 0:num_sectors-1
                    m2 = mm2+1;
                    r = floor((rr)/R_fac)+1;
                    delays = ((-1)^(rr))*(1/(R_fac^stage))*eta*[2^(stage-1)*(2*mm2+0.5)-M/2; 2^(stage-1)*(2*mm2+1.5)-M/2];
                    y_interp = reshape(v_prev([2*mm2+1,2*mm2+2],:,r),[2,N]);
                    % v_prev_interped = TTD_Beamformer(y_interp.',-delays, fs, @Psi_sinc, 'R',taps2).';
                    % each stage filter design 
                    init_filter_fast_ds2 = tic;
                    [h_stage, h_stage_idx] = out2.Filter_Gen(-delays);
                    y_stage = (y_interp).';
                    N_filter = N+2*taps2+1;
                    h_stage = squeeze(h_stage);
                    h_stage_idx = squeeze(h_stage_idx);
                    
                    wrapN = @(x, n) (1 + mod(x-1, n));
                    
                    H_stage = zeros(N_filter,L_fac);
                    
                    for i=1:L_fac
                        H_stage(wrapN(h_stage_idx(i,:),N_filter),i) = h_stage(i,:);
                    end
                    Hf_stage = fft(H_stage,[],1);
                    count_filter_time2 = count_filter_time2 + toc(init_filter_fast_ds2);

                    Xf_stage = fft(y_stage,N_filter,1);
                    Y_stage = ifft(Hf_stage.*Xf_stage,[],1);
                    y_filter_stage = Y_stage(1:N,:).';

                    v_now(m2,:,rr+1) = sum(y_filter_stage,1)./2;
                end
            end
        end
        runtime_fast_ds_32_tap = toc(init_fast_ds2);
        Runtime_fast_ds_32_tap(jj,mm1) = runtime_fast_ds_32_tap - count_filter_time2;


        %% traditional delay and sum beamforming 16 taps
        t_init_ds_filter = tic;
        N_filter = N+2*taps+1;
        
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
        runtime_ds_filter = toc(t_init_ds_filter);
        
        t_init_beam_ds = tic;
        y_ttd_in = (y.*mod_phase).';
        % y_filter_ttd = TTD_Beamformer((y.*mod_phase).',-tau, fs, @Psi_sinc,h, h_idx, 'R',R).';
        Xf = fft(y_ttd_in,N_filter,1);
        Y = ifft(Hf.*Xf,[],1);
        y_filter_ttd = Y(1:N,:).';
        y_bf_ttd = reshape(sum(y_filter_ttd(:,1:N),1)/M, [N,1]);
        runtime_ds_beam = toc(t_init_beam_ds);

        Runtime_ds_16_tap(jj,mm1) = R*runtime_ds_beam;

        %% traditional delay and sum beamforming 32 taps
        t_init_ds_filter2 = tic;

        N_filter2 = N+2*taps2+1;
        
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
        runtime_ds_filter2 = toc(t_init_ds_filter2);

        t_init_beam_ds2 = tic;
        y_ttd_in2 = (y.*mod_phase).';
        % y_filter_ttd = TTD_Beamformer((y.*mod_phase).',-tau, fs, @Psi_sinc,h, h_idx, 'R',R).';
        Xf2 = fft(y_ttd_in2,N_filter2,1);
        Y2 = ifft(Hf2.*Xf2,[],1);
        y_filter_ttd2 = Y2(1:N,:).';
        y_bf_ttd2 = reshape(sum(y_filter_ttd2(:,1:N),1)/M, [N,1]);
        runtime_ds_beam2 = toc(t_init_beam_ds2);

        Runtime_ds_32_tap(jj,mm1) = R*runtime_ds_beam2;
    end
end

figure(1)
semilogy(M_lst,mean(Runtime_fbst_superfast),'-*','LineWidth',1)
hold on
grid on
% plot(B_lst, mean(Runtime_fbst2), 'LineWidth',1)
% hold on 
% grid on
% plot(B_lst, mean(Runtime_fbst_deconv),'-*','LineWidth',1)
% hold on
% grid on
semilogy(M_lst, mean(Runtime_fbst_precompute),'-*','LineWidth',1);
hold on
grid on
% plot(B_lst, mean(Runtime_fbst_slepian_inv),'LineWidth',1)
% hold on
% grid on
semilogy(M_lst,mean(Runtime_ds_16_tap),'-^','LineWidth',1)
hold on
grid on
semilogy(M_lst, mean(Runtime_ds_32_tap),'-^','LineWidth',1)
hold on
grid on
semilogy(M_lst, mean(Runtime_fast_ds_16_tap),'-^','LineWidth',1)
hold on 
grid on
semilogy(M_lst, mean(Runtime_fast_ds_32_tap),'-^','LineWidth',1)
xlim([M_lst(1) M_lst(end)])
set(gca, 'XScale', 'log'); % Set x-axis to logarithmic scale
xticks(M_lst); % Set specific log2 tick marks
xticklabels({'2^1','2^2','2^3','2^4','2^5', '2^6','2^7', '2^8','2^9','2^{10}','2^{11}','2^{12}'}); % Label as powers of 2
xlabel('Array size (M)','Interpreter','latex')
ylabel('Runtime (s)','Interpreter','latex')
legend({'FBST Superfast','FBST Precompute','Delay and Sum (R=16)', 'Delay and Sum (R=32)', 'Fast Delay and Sum (R=16)', 'Fast Delay and Sum (R=32)'},'Location','northwest','Interpreter','latex','FontSize',12)
% exportgraphics(gcf, 'ula_runtime_vs_M.pdf', 'ContentType', 'vector');

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