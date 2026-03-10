clear
clc
close all

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

M_lst = [2^1 2^2 2^3 2^4 2^5 2^6 2^7 2^8 2^9 2^10 2^11 2^12];
exp_fac = 1.05;
trials = 25;

Runtime_fbst_superfast = zeros(trials, length(M_lst));
Runtime_fbst_precompute = zeros(trials, length(M_lst));
Runtime_ds_16_tap = zeros(trials, length(M_lst));
Runtime_ds_32_tap = zeros(trials, length(M_lst));
Runtime_fast_ds_16_tap = zeros(trials, length(M_lst));
Runtime_fast_ds_32_tap = zeros(trials, length(M_lst));

Runtime_DFT_beamformer = zeros(trials, length(M_lst));

Runtime_chirp_z = zeros(trials, length(M_lst));

for mm1 = 1:length(M_lst)
    M = M_lst(mm1);
    % N = 2^6;
    B = 2^(log2(M)-1); % create a total of M+1 fbst beams
    R = 2^(log2(M)); % beams for ds and fast-ds
    N=2^5;
    L = ceil(N); % no. of samples in frequency is 2*L
    MN =M*N;
    

    Theta = zeros(2*B+1,1);
    idx = 2*B;%randi(B,1);
    for i=1:2*B+1
        Theta(i) = asin((i-1-B)/B);
    end

    fc = 20e9;  % center frequency
    Ws = 5e9; % to prevent spacial aliasing
    fc_tilde = (fc+Ws); % adjustment for spatial aliasing
    N = max(32,ceil(1/(2/exp_fac-1)*((M-1)*Ws/(fc+Ws)-1)));
    MN = M*N;
    L = ceil(N); % no. of samples in frequency is 2*L
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

    freq_resolution = fs/N; % frequency grid for sub-band beamformign

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
    freq_samples = exp_fac/(2*T*L)*linspace(-L,L,2*L+1);
    % F_basis = exp(1i*2*pi*t*freq_samples).*demod_phase;
    % A_eplitz = zeros(2*L+1,2*L+1);
    % MN_blk = 60000;
    % if MN<=MN_blk
    %     F_basis = exp(1i*2*pi*t*freq_samples).*demod_phase;
    %     A_eplitz = F_basis'*F_basis;
    % else
    %     blks = floor(MN/MN_blk);
    %     for bb=1:blks+1
    %         if bb<=blks
    %             tb = t((bb-1)*MN_blk+1:bb*MN_blk);
    %             dp_b = demod_phase((bb-1)*MN_blk+1:bb*MN_blk);
    %         else
    %             tb = t((bb-1)*MN_blk+1:end);
    %             dp_b = demod_phase((bb-1)*MN_blk+1:end);
    %         end
    %         F_basis_batch = exp(1i*2*pi*tb*freq_samples).*dp_b;
    %         A_eplitz = A_eplitz + F_basis_batch'*F_basis_batch;
    %     end
    % end
    % toeplitz_inv = inv(A_eplitz + (delta)*eye(2*L+1));
    % % % save data
    % % filename = sprintf('A_eplitz_and_inv_%d.mat', MN);
    % save(['runtime_ULA_data\' filename],'A_eplitz','toeplitz_inv');
    % load data
    load(sprintf('runtime_ULA_data/A_eplitz_and_inv_%d.mat',MN));

    %% get canonical vectors of Toeplitz inverse via brute force
    [X1,X2,X3,X4] = eval_canonical_vecs_brute_force(A_eplitz + delta*eye(2*L+1),2*L+1);
    x0 = X1(1,1);

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

    combined_filter = pf_x1.*pf_x2 - pf_x3.*pf_x4;

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
    A = exp(-1i*pi*exp_fac);%exp(-1i*2*pi*T*W);
    W_czt = exp(-1i*pi*exp_fac/L);
    %------ filter design for chirp-z stage-1 (taken for each row of y)
    m = N; % input signal length
    k1 = 2*L+1; % chirp-z length
    n = 1;
    nfft = 2^nextpow2(m+k1-1);
    %------- Premultiply data.
    kk = ((-m+1):max(k1-1,m-1)).';
    kk2 = (kk .^ 2) ./ 2;
    ww = W_czt .^ (kk2);   % <----- Chirp filter is 1./ww
    nn = (0:(m-1))';
    aa_stage1 = A .^ ( -nn );
    aa_stage1 = aa_stage1.*ww(m+nn);

    fv_stage1 = fft( 1 ./ ww(1:(k1-1+m)), nfft );   % <----- Chirp filter.

    %---------- Filter design for stage-2 (taken for each column of F_transform)
    m_len2 = M; % input signal length
    k1_2 = 2*B+1; % chirp-z length
    n = 1;
    nfft_2 = 2^nextpow2(m_len2+k1_2-1);
    %------- Premultiply data.
    kk = ((-m_len2+1):max(k1_2-1,m_len2-1)).';
    kk2 = (kk .^ 2) ./ 2;
    nn = (0:(m_len2-1))';

    A_coefs_stage2 = zeros(2*L+1,1);
    W_czt_coefs_stage2 = zeros(2*L+1,1);
    FV_filters = zeros(nfft_2,2*L+1);
    ww_stage2 = zeros(length(kk2),2*L+1);
    aa_stage2 = zeros(m_len2,2*L+1);

    for ll=1:2*L+1
        l_dash = ll-1;
        A_coefs_stage2(ll,1) = exp(1i*pi*(1/B)*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L))*B);
        aa = A_coefs_stage2(ll,1) .^ ( -nn );
        W_czt_coefs_stage2(ll,1) = exp(1i*pi*(1/B)*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L)));
        ww_stage2(:,ll) = W_czt_coefs_stage2(ll,1) .^ (kk2);   % <----- Chirp filter is 1./ww
        aa_stage2(:,ll) = aa.*ww_stage2(m_len2+nn,ll);
        FV_filters(:,ll) = fft( 1 ./ ww_stage2(1:(k1_2-1+m_len2),ll), nfft_2 );   % <----- Chirp filter.
    end

    SNR = 30;
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

        s_nyq = sig_gen(a,f,n_sinusoids,t_nyq'); % for Array gain test

        %% chirp-z transform
        t_init_chirpz = tic;
        F_transform = zeros(M,2*L+1);
        for mm=1:M
            y_in = y(mm,:).' .* aa_stage1(:,ones(1,n));
            fy = fft(  y_in, nfft );
            fy = fy .* fv_stage1(:,ones(1, n));
            g  = ifft( fy );
            g = g( m:(m+k1-1), :, :) .* ww( m:(m+k1-1),ones(1, n) );
            F_transform(mm,:) = g;
        end

        X_czt = zeros(2*L+1,2*B+1);
        B_list = reshape(linspace(0,2*B,2*B+1), [1,2*B+1]);
        for ll=1:2*L+1
            l_dash = ll-1;
            y_in = F_transform(:,ll) .* aa_stage2(:,ll);
            fy = fft(  y_in, nfft_2 );
            fv = FV_filters(:,ll);
            fy = fy .* fv(:,ones(1, n));
            g  = ifft( fy );
            g = g( m_len2:(m_len2+k1_2-1), :, :) .* ww_stage2( m_len2:(m_len2+k1_2-1),ll );
            coeff1 = exp(-1i*pi*(1/B)*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L))*(M/2-1/2)*B_list); % adjusting for array center
            coeff2 = exp(1i*pi*(1/B)*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L))*(M/2-1/2)*B); % adjusting for negative angles
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
        fft_wl = fft([w_l; zeros(length(w_l), 1)]);
        b1 = ifft(pf_x2 .* fft_wl);
        b1 = b1(1:length(w_l));
        % b1 = toep_mult(X2,w_l);
        b2 = ifft(pf_x4 .* fft_wl);
        b2 = b2(1:length(w_l));
        % b2 = toep_mult(X4,w_l);
        b3 = ifft(pf_x1 .* fft([b1; zeros(length(b1), 1)]));
        b3 = b3(1:length(b1));
        % b3 = toep_mult(X1,b1);
        b4 = ifft(pf_x3 .* fft([b2; zeros(length(b2), 1)]));
        b4 = b4(1:length(b2));
        % b4 = toep_mult(X3,b2);
        % b_combined = ifft(combined_filter.*fft_wl);
        % b_combined = b_combined(1:length(w_l));
        alpha_superfast = (b3 - b4)/x0;
        y_nyq_superfast = psi*alpha_superfast;
        runtime_fbst_superfast_beam = toc(t_init_superfast_fbst);

        Runtime_fbst_precompute(jj,mm1) = (runtime_chirpz + (2*B+1)*runtime_fbst_precompute_beam)./N;
        Runtime_fbst_superfast(jj,mm1) = (runtime_chirpz + (2*B+1)*runtime_fbst_superfast_beam)./N;
        Runtime_chirp_z(jj,mm1) = runtime_chirpz;

        AG_fbst_precompute = db(norm(s_nyq)/norm(s_nyq-y_nyq)) - SNR_dB;
        AG_fbst_superfast = db(norm(s_nyq)/norm(s_nyq - y_nyq_superfast)) - SNR_dB;
        fprintf("AG Ideal: %.2f AG FBST Precompute: %.2f AG Superfast: %.2f\n",db(M)/2,AG_fbst_precompute,AG_fbst_superfast);

        %% fast delay and sum beamforming 16 taps
        count_filter_time = 0; % counter for subtracting filter initialization
        init_fast_ds = tic;
        v_now = reshape(y.*mod_phase, [M,N,1]);
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
        Runtime_fast_ds_16_tap(jj,mm1) = (runtime_fast_ds_16_tap - count_filter_time)./N;

        %% fast delay and sum beamforming 32 taps
        count_filter_time2 = 0; % counter for subtracting filter initialization

        init_fast_ds2 = tic;
        v_now = reshape(y.*mod_phase, [M,N,1]);
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
        Runtime_fast_ds_32_tap(jj,mm1) = (runtime_fast_ds_32_tap - count_filter_time2)./N;


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

        Runtime_ds_16_tap(jj,mm1) = R*runtime_ds_beam./N;

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

        Runtime_ds_32_tap(jj,mm1) = R*runtime_ds_beam2./N;


        %% sub-band processing beamformer
        t_init_subband = tic;
        Y_f = zeros(N,M);
        for mmm=1:M
            y_temp = y(mmm,:);
            Y_f(:,mmm) = fft(y_temp);
        end

        % Y_f = fft(y.'); % mth column for mth array sample

        Y_nb = zeros(N,1);
        f_set = zeros(N,1);


        Y_freq_beamspace = zeros(N,2*B+1);
        beam_list = linspace(0,2*B,2*B+1);

        for nn=0:N-1
            ym = Y_f(nn+1,:).';
            if nn<N/2
                f_curr = fc + nn*freq_resolution;
            elseif nn>=N/2
                f_curr = fc + (nn)*freq_resolution - N*freq_resolution;
            end
            W_czt = exp(1i*pi*f_curr/((fc+W)*B));
            A_czt = exp(1i*pi*f_curr/(fc+W));
            coeff2 = exp(1i*pi*f_curr*(M/2-1/2)/(fc+W));
            coeff1 = exp(-1i*pi*f_curr*((M/2-1/2)/(fc+W))*(beam_list/B)).';
            temp_var = czt(ym,2*B+1,W_czt,A_czt);
            Y_freq_beamspace(nn+1,:) = coeff2*temp_var.*coeff1;
        end

        Y_time_beamspace = zeros(N,2*B+1);

        for bbb=1:2*B+1
            Y_time_beamspace(:,bbb) = ifft(Y_freq_beamspace(:,bbb))/M;
        end
        % 
        % Y_time_beamspace = ifft(Y_freq_beamspace);
        % Y_bf = Y_time_beamspace(:,idx)/M; % 
        runtime_subband = toc(t_init_subband);
        Runtime_DFT_beamformer(jj,mm1) = runtime_subband./N;
    end
end

figure(1)
semilogy(M_lst(1:end-1),mean(Runtime_fbst_superfast(:,1:end-1)),'-*','LineWidth',1)
hold on
grid on
% plot(B_lst, mean(Runtime_fbst2), 'LineWidth',1)
% hold on 
% grid on
% plot(B_lst, mean(Runtime_fbst_deconv),'-*','LineWidth',1)
% hold on
% grid on
semilogy(M_lst(1:end-1), mean(Runtime_fbst_precompute(:,1:end-1)),'-*','LineWidth',1);
hold on
grid on
% plot(B_lst, mean(Runtime_fbst_slepian_inv),'LineWidth',1)
% hold on
% grid on
semilogy(M_lst(1:end-1),mean(Runtime_ds_16_tap(:,1:end-1)),'-^','LineWidth',1)
hold on
grid on
semilogy(M_lst(1:end-1), mean(Runtime_ds_32_tap(:,1:end-1)),'-^','LineWidth',1)
hold on
grid on
semilogy(M_lst(1:end-1), mean(Runtime_fast_ds_16_tap(:,1:end-1)),'-^','LineWidth',1)
hold on 
grid on
semilogy(M_lst(1:end-1), mean(Runtime_fast_ds_32_tap(:,1:end-1)),'-^','LineWidth',1)
hold on
grid on
semilogy(M_lst(1:end-1),mean(Runtime_DFT_beamformer(:,1:end-1)),'-^','LineWidth',1)
xlim([M_lst(1) M_lst(end-1)])
set(gca, 'XScale', 'log'); % Set x-axis to logarithmic scale
xticks(M_lst); % Set specific log2 tick marks
xticklabels({'2^1','2^2','2^3','2^4','2^5', '2^6','2^7', '2^8','2^9','2^{10}','2^{11}'}); % Label as powers of 2
xlabel('Array size (M)','Interpreter','latex')
ylabel('Runtime (s)','Interpreter','latex')
legend({'FBST Superfast','FBST Precompute','Delay and Sum (R=16)', 'Delay and Sum (R=32)', 'Fast Delay and Sum (R=16)', 'Fast Delay and Sum (R=32)','Sub-Band Processing'},'Location','northwest','Interpreter','latex','FontSize',12)
% exportgraphics(gcf, 'ula_runtime_vs_M_subband.pdf', 'ContentType', 'vector');

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