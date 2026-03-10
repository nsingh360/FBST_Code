clear
clc
close all

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

exp_fac = 1.05;
M_lst = [2^2 2^4 2^6 2^8 2^10 2^12];
trials = 25;


fc = 20e9;  % center frequency
Ws = 5e9;
fc_tilde = fc+Ws;
fc_norm = fc/fc_tilde;
c = physconst('LightSpeed');
lambda = c/(fc+Ws); % wavelngth
W = 5e9; % bandwidth
fs = 2*W; % sampling frequency
T = 1/fs;

Runtime_fbst_superfast = zeros(trials, length(M_lst));
Runtime_fbst_precompute = zeros(trials, length(M_lst));
Runtime_ds_16_tap = zeros(trials, length(M_lst));
Runtime_ds_32_tap = zeros(trials, length(M_lst));
Runtime_fast_ds_16_tap = zeros(trials, length(M_lst));
Runtime_fast_ds_32_tap = zeros(trials, length(M_lst));
Runtime_chirp_z = zeros(trials, length(M_lst));

for m_num = 1:length(M_lst)

    M = M_lst(m_num); % sqrt(M) x sqrt(M) planar array
    B = ceil((-3 + sqrt(9 + 8*(M-1)))/4); % approximately M beams are formed for FBST
    R = 2^(log2(sqrt(M))); % beams for ds and fast-ds are R^2
    N = max(32,ceil(2/(2/exp_fac-1)*(sqrt(2)*(sqrt(M)-1)*W/(fc+W)-1)));
    L = ceil(N); % no. of samples in frequency is 2*L+1
    MN = M*N;

    J = linspace(0,B,B+1);
    K = linspace(-B,B,2*B+1);
    j_idx = randsample(1:B+1,1);
    k_idx = randsample(1:2*B+1,1);
    theta = asin(sqrt(J(j_idx)^2 + K(k_idx)^2)/(B));
    phi = asin(K(k_idx)/(sqrt(J(j_idx)^2 + K(k_idx)^2)));

    while ((~isreal(theta) || ~isreal(phi)) || (isnan(theta) || isnan(phi)))
        j_idx = randsample(1:B+1,1);
        k_idx = randsample(1:2*B+1,1);
        theta = asin(sqrt(J(j_idx)^2 + K(k_idx)^2)/(B));
        phi = asin(K(k_idx)/(sqrt(J(j_idx)^2 + K(k_idx)^2)));
    end

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
    taps = 8;
    edge_lim = taps + ceil(fs*(max(tau)-min(tau))/2);% clip signal to account for edge distortion
    comp_idxs = (edge_lim+1):(N-edge_lim-1); % comparison indicies

    taps2 = 16;
    edge_lim2 = taps2 + ceil(fs*(max(tau)-min(tau))/2);% clip signal to account for edge distortion
    comp_idxs2 = (edge_lim2+1):(N-edge_lim2-1); % comparison indicies

    %% parameters for fast delay and sum
    S = log2(sqrt(M)); % 2 S stage 1D FDS
    R_fac = 2; % 2-radix fast delay and sum
    L_fac = 2; % downsampling factor in spatial dimension
    eta = lambda/(2*c);

    %% construct the toeplitz matrix
    delta = 1e-5;
    freq_samples = exp_fac/(2*T*L)*linspace(-L,L,2*L+1);
    F_basis = exp(1i*2*pi*t*freq_samples).*demod_phase;
    A_eplitz = F_basis'*F_basis;
    toeplitz_inv = inv(A_eplitz + (delta)*eye(2*L+1));

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

    %% equispaced sampled basis matrix
    t_nyq = [0:N-1]/fs; % nyquist samples for testing
    psi = exp(1i*2*pi*t_nyq'*freq_samples);

    %% interpolator for TTD beamformer
    [out] = Kernel_Filter_Gen(fs, @Psi_sinc,'R',taps);
    [h, h_idx] = out.Filter_Gen(-tau);

    [out2] = Kernel_Filter_Gen(fs, @Psi_sinc,'R',taps2);
    [h2, h_idx2] = out2.Filter_Gen(-tau);

    %% chirp filter initialization
    %------- first stage chirp-z
    A_stage1 = exp(-1i*pi*exp_fac);%exp(-1i*2*pi*T*W);
    W_czt_stage1 = exp(-1i*pi*exp_fac/L);
    m_stage1 = N; % input signal length
    k1_1 = 2*L+1; % chirp-z length
    n_len = sqrt(M);
    nfft = 2^nextpow2(m_stage1+k1_1-1);
    %------- Premultiply data.
    kk = ((-m_stage1+1):max(k1_1-1,m_stage1-1)).';
    kk2 = (kk .^ 2) ./ 2;
    ww = W_czt_stage1 .^ (kk2);   % <----- Chirp filter is 1./ww
    nn = (0:(m_stage1-1))';
    aa = A_stage1 .^ ( -nn );
    aa_stage1 = aa.*ww(m_stage1+nn);
    
    fv_stage1 = fft( 1 ./ ww(1:(k1_1-1+m_stage1)), nfft );   % <----- Chirp filter.
    
    %--------- second stage chirp-z 
    m_stage2 = sqrt(M); % input signal length
    k1_2 = B+1; % chirp-z length
    n_len2 = sqrt(M);
    nfft_2 = 2^nextpow2(m_stage2+k1_2-1);
    
    A_stage2 = 1;
    kk = ((-m_stage2+1):max(k1_2-1,m_stage2-1)).';
    kk2 = (kk .^ 2) ./ 2;
    nn = (0:(m_stage2-1))';
    
    W_czt_stage2_coeffs = zeros(2*L+1,1);
    aa_stage2 = zeros(m_stage2,2*L+1);
    FV_stage2 = zeros(nfft_2,2*L+1);
    ww_stage2 = zeros(length(kk2),2*L+1);
    for ll=1:2*L+1
        l_dash = ll - 1;
        W_czt_stage2_coeffs(ll,1) = exp(1i*pi*(1/(B))*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L)));
        ww_stage2(:,ll) = W_czt_stage2_coeffs(ll,1) .^ (kk2);   % <----- Chirp filter is 1./ww
        aa = A_stage2 .^ ( -nn );
        aa_stage2(:,ll) = aa.*ww_stage2(m_stage2+nn,ll);
        FV_stage2(:,ll) = fft( 1 ./ ww_stage2(1:(k1_2-1+m_stage2),ll), nfft_2 );   % <----- Chirp filter.
    end
    
    %------------ third stage chirp-z 
    m_stage3 = sqrt(M); % input signal length
    k1_3 = 2*B+1;
    n_len3 = B+1;
    nfft_3 = 2^nextpow2(m_stage3+k1_3-1);
    kk = ((-m_stage3+1):max(k1_3-1,m_stage3-1)).';
    kk2 = (kk .^ 2) ./ 2;
    nn = (0:(m_stage3-1))';
    
    A_coeff_stage3 = zeros(2*L+1,1);
    W_czt_stage3_coeffs = zeros(2*L+1,1);
    aa_stage3 = zeros(m_stage3,2*L+1);
    FV_stage3 = zeros(nfft_3,2*L+1);
    ww_stage3 = zeros(length(kk2),2*L+1);
    for ll=1:2*L+1
        l_dash = ll - 1;
        A_coeff_stage3(ll,1) = exp(1i*pi*(1/(B))*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L))*B);
        W_czt_stage3_coeffs(ll,1) = exp(1i*pi*(1/(B))*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L)));
        ww_stage3(:,ll) = W_czt_stage3_coeffs(ll,1) .^ (kk2);   % <----- Chirp filter is 1./ww
        aa = A_coeff_stage3(ll,1) .^ ( -nn );
        aa_stage3(:,ll) = aa.*ww_stage3(m_stage3+nn,ll);
        FV_stage3(:,ll) = fft( 1 ./ ww_stage3(1:(k1_3-1+m_stage3),ll), nfft_3 );   % <----- Chirp filter.
    end

    SNR = 30;
    for jj=1:trials

        fprintf('Trial: %d, Array size.: %d\n',jj,M_lst(m_num))
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

        s_nyq = sig_gen(a,f,n_sinusoids,t_nyq'); % for Array gain test

        %% runtime for FBST 
        %% chirp-z transform 
        t_init_chirp_z = tic;
        F_transform = zeros(sqrt(M),sqrt(M),2*L+1); % transform w.r.t. temporal frequency
        % A = -1;%exp(-1i*2*pi*T*W);
        % W_czt = exp(-1i*pi/L);
        for mm1=1:sqrt(M)
            y_in = reshape(y_3d(mm1,:,:),n_len,N).';
            y_in = y_in .* aa_stage1(:,ones(1,n_len));
            fy = fft(  y_in, nfft );
            fy = fy .* fv_stage1(:,ones(1, n_len));
            g  = ifft( fy );
            g = g( m_stage1:(m_stage1+k1_1-1), :, :) .* ww( m_stage1:(m_stage1+k1_1-1),ones(1, n_len) );
            F_transform(mm1,:,:) = g.';
        end

        X2_transform = zeros(sqrt(M),B+1,2*L+1); % transform w.r.t. second spatial frequency
        A= 1;
        B_list1 = reshape(linspace(0,2*B,2*B+1), [1,2*B+1]);
        B_list2 = reshape(linspace(0,B,B+1), [1,B+1]);
        for ll=1:2*L+1
                l_dash = ll-1;

                y_in = reshape(F_transform(:,:,ll),sqrt(M),sqrt(M)).';
                y_in = y_in .* reshape(aa_stage2(:,ll,ones(1,n_len2)),sqrt(M),n_len2);
                fy = fft(  y_in, nfft_2 );
                fy = fy .* reshape(FV_stage2(:,ll,ones(1, n_len2)),nfft_2,n_len2);
                g  = ifft( fy );
                g = g( m_stage2:(m_stage2+k1_2-1), :, :) .* reshape(ww_stage2( m_stage2:(m_stage2+k1_2-1),ll,ones(1, n_len2) ),k1_2,n_len2);

                % A = exp(1i*pi*(1/(B))*(fc_norm + (1/(2*fc_tilde*T*L))*(l_dash-L))*B);%exp(-1i*pi*(1 + 2*(W/fc)*(l_dash-L/2)/L));
                % W_czt = exp(1i*pi*(1/(B))*(fc_norm + (1/(2*fc_tilde*T*L))*(l_dash-L)));
                coeff = exp(-1i*pi*(1/(B))*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L))*(sqrt(M)/2-1/2)*B_list2); % adjusting for array center
                % coeff2 = exp(1i*pi*(1/B)*(fc_norm + (1/(2*fc_tilde*T*L))*(l_dash-L))*(sqrt(M)/2-1/2)*B); % adjusting for negative azimuthal
                X2_transform(:,:,ll) = g.';
                X2_transform(:,:,ll) = X2_transform(:,:,ll).*coeff;
        end

        X1_transform = zeros(2*B+1,B+1,2*L+1); % final transform matrix (after taking czt along the first spatial frequency)
        for ll=1:2*L+1
                l_dash = ll-1;
                
                y_in = reshape(X2_transform(:,:,ll),sqrt(M),B+1);
                y_in = y_in .* reshape(aa_stage3(:,ll,ones(1,n_len3)),sqrt(M),n_len3);
                fy = fft(  y_in, nfft_3 );
                fy = fy .* reshape(FV_stage3(:,ll,ones(1, n_len3)),nfft_3,n_len3);
                g  = ifft( fy );
                g = g( m_stage3:(m_stage3+k1_3-1), :, :) .* reshape(ww_stage3( m_stage3:(m_stage3+k1_3-1),ll,ones(1, n_len3) ),k1_3,n_len3);

                % A = exp(1i*pi*(1/(B))*(fc_norm + (1/(2*fc_tilde*T*L))*(l_dash-L))*B); % check
                % W_czt = exp(1i*pi*(1/(B))*(fc_norm + (1/(2*fc_tilde*T*L))*(l_dash-L)));
                coeff = exp(-1i*pi*(1/(B))*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L))*(sqrt(M)/2-1/2)*B_list1); % adjusting for array center
                coeff2 = exp(1i*pi*(1/B)*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L))*(sqrt(M)/2-1/2)*B); % adjusting for negative azimuthal
                X1_transform(:,:,ll) = g;
                X1_transform(:,:,ll) = coeff2*X1_transform(:,:,ll).*coeff.';
        end

        w_l = reshape(X1_transform(k_idx,j_idx,:),[2*L+1,1]);
        runtime_chirp_z = toc(t_init_chirp_z);

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
        alpha_superfast = (b3 - b4)/x0;
        y_nyq_superfast = psi*alpha_superfast;
        runtime_fbst_superfast_beam = toc(t_init_superfast_fbst);

        Runtime_fbst_precompute(jj,m_num) = (runtime_chirp_z + (B+1)*(2*B+1)*runtime_fbst_precompute_beam)/N;
        Runtime_fbst_superfast(jj,m_num) = (runtime_chirp_z + (B+1)*(2*B+1)*runtime_fbst_superfast_beam)/N;
        Runtime_chirp_z(jj,m_num) = runtime_chirp_z;

        AG_fbst_precompute = db(norm(s_nyq)/norm(s_nyq-y_nyq)) - SNR_dB;
        AG_fbst_superfast = db(norm(s_nyq)/norm(s_nyq - y_nyq_superfast)) - SNR_dB;
        fprintf("AG Ideal: %.2f AG FBST Precompute: %.2f AG Superfast: %.2f\n",db(M)/2,AG_fbst_precompute,AG_fbst_superfast);


        %% runtime for fast delay and sum 16 tap
        init_fast_ds = tic;
        y_multibeamformed_first = zeros(sqrt(M),N,1,R_fac^S);
        counter_loop1 = 0;
        %% first FDS along b1 
        for m_a1 = 1:sqrt(M)
            v_now = reshape(y_3d(m_a1,:,:,1,1),[sqrt(M),N,1]);
            for stage=1:S
                v_prev = v_now;
                num_virtual_sensors = sqrt(M)/2^stage;
                num_sectors = R_fac^stage;
                v_now = zeros(num_virtual_sensors,N,num_sectors); % init partial beamforming at stage s
    
                for mm=0:num_virtual_sensors-1
                    for rr = 0:num_sectors-1
                        m2 = mm+1;
                        r = floor((rr)/R_fac)+1;
                        delays = ((-1)^(rr))*(1/(R_fac^stage))*eta*[2^(stage-1)*(2*mm+0.5)-sqrt(M)/2; 2^(stage-1)*(2*mm+1.5)-sqrt(M)/2];
                        y_interp = reshape(v_prev([2*mm+1,2*mm+2],:,r),[2,N]);
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
                        counter_loop1 = counter_loop1 + toc(init_filter_fast_ds);
    
                        Xf_stage = fft(y_stage,N_filter,1);
                        Y_stage = ifft(Hf_stage.*Xf_stage,[],1);
                        y_filter_stage = Y_stage(1:N,:).';
                        v_now(m2,:,rr+1) = sum(y_filter_stage,1)./2;
                    end
                end
            end
            y_multibeamformed_first(m_a1,:,1,:) = v_now;
        end

        y_multibeamformed_total = zeros(N,R_fac^S,R_fac^S);
        counter_loop2 = 0;
        %% second FDS along a1
        for r_b1=1:R_fac^S
            v_now = reshape(y_multibeamformed_first(:,:,1,r_b1),[sqrt(M),N,1]);
            for stage=1:S
                v_prev = v_now;
                num_virtual_sensors = sqrt(M)/2^stage;
                num_sectors = R_fac^stage;
                v_now = zeros(num_virtual_sensors,N,num_sectors); % init partial beamforming at stage s
    
                for mm=0:num_virtual_sensors-1
                    for rr = 0:num_sectors-1
                        m2 = mm+1;
                        r = floor((rr)/R_fac)+1;
                        delays = ((-1)^(rr))*(1/(R_fac^stage))*eta*[2^(stage-1)*(2*mm+0.5)-sqrt(M)/2; 2^(stage-1)*(2*mm+1.5)-sqrt(M)/2];
                        y_interp = reshape(v_prev([2*mm+1,2*mm+2],:,r),[2,N]);
                        % v_prev_interped = TTD_Beamformer(y_interp.',-delays, fs, @Psi_sinc, 'R',taps2).';
                        % each stage filter design
                        init_filter_fast_ds_b1 = tic;
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
                        counter_loop2 = counter_loop2 + toc(init_filter_fast_ds_b1);
    
                        Xf_stage = fft(y_stage,N_filter,1);
                        Y_stage = ifft(Hf_stage.*Xf_stage,[],1);
                        y_filter_stage = Y_stage(1:N,:).';
                        v_now(m2,:,rr+1) = sum(y_filter_stage,1)./2;
                    end
                end
            end
            y_multibeamformed_total(:,:,r_b1) = v_now;
        end
        runtime_fast_ds = toc(init_fast_ds);
        Runtime_fast_ds_16_tap(jj,m_num) = (runtime_fast_ds - counter_loop1 - counter_loop2)/N; % subtract filter init time


        %% runtime for fast delay and sum 32 tap
        init_fast_ds_2 = tic;
        y_multibeamformed_first = zeros(sqrt(M),N,1,R_fac^S);
        counter_loop1 = 0;
        %% first FDS along b1 
        for m_a1 = 1:sqrt(M)
            v_now = reshape(y_3d(m_a1,:,:,1,1),[sqrt(M),N,1]);
            for stage=1:S
                v_prev = v_now;
                num_virtual_sensors = sqrt(M)/2^stage;
                num_sectors = R_fac^stage;
                v_now = zeros(num_virtual_sensors,N,num_sectors); % init partial beamforming at stage s
    
                for mm=0:num_virtual_sensors-1
                    for rr = 0:num_sectors-1
                        m2 = mm+1;
                        r = floor((rr)/R_fac)+1;
                        delays = ((-1)^(rr))*(1/(R_fac^stage))*eta*[2^(stage-1)*(2*mm+0.5)-sqrt(M)/2; 2^(stage-1)*(2*mm+1.5)-sqrt(M)/2];
                        y_interp = reshape(v_prev([2*mm+1,2*mm+2],:,r),[2,N]);
                        % v_prev_interped = TTD_Beamformer(y_interp.',-delays, fs, @Psi_sinc, 'R',taps2).';
                        % each stage filter design
                        init_filter_fast_ds = tic;
                        [h_stage, h_stage_idx] = out2.Filter_Gen(-delays);
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
                        counter_loop1 = counter_loop1 + toc(init_filter_fast_ds);
    
                        Xf_stage = fft(y_stage,N_filter,1);
                        Y_stage = ifft(Hf_stage.*Xf_stage,[],1);
                        y_filter_stage = Y_stage(1:N,:).';
                        v_now(m2,:,rr+1) = sum(y_filter_stage,1)./2;
                    end
                end
            end
            y_multibeamformed_first(m_a1,:,1,:) = v_now;
        end

        y_multibeamformed_total = zeros(N,R_fac^S,R_fac^S);
        counter_loop2 = 0;
        %% second FDS along a1
        for r_b1=1:R_fac^S
            v_now = reshape(y_multibeamformed_first(:,:,1,r_b1),[sqrt(M),N,1]);
            for stage=1:S
                v_prev = v_now;
                num_virtual_sensors = sqrt(M)/2^stage;
                num_sectors = R_fac^stage;
                v_now = zeros(num_virtual_sensors,N,num_sectors); % init partial beamforming at stage s
    
                for mm=0:num_virtual_sensors-1
                    for rr = 0:num_sectors-1
                        m2 = mm+1;
                        r = floor((rr)/R_fac)+1;
                        delays = ((-1)^(rr))*(1/(R_fac^stage))*eta*[2^(stage-1)*(2*mm+0.5)-sqrt(M)/2; 2^(stage-1)*(2*mm+1.5)-sqrt(M)/2];
                        y_interp = reshape(v_prev([2*mm+1,2*mm+2],:,r),[2,N]);
                        % v_prev_interped = TTD_Beamformer(y_interp.',-delays, fs, @Psi_sinc, 'R',taps2).';
                        % each stage filter design
                        init_filter_fast_ds_b1 = tic;
                        [h_stage, h_stage_idx] = out2.Filter_Gen(-delays);
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
                        counter_loop2 = counter_loop2 + toc(init_filter_fast_ds_b1);
    
                        Xf_stage = fft(y_stage,N_filter,1);
                        Y_stage = ifft(Hf_stage.*Xf_stage,[],1);
                        y_filter_stage = Y_stage(1:N,:).';
                        v_now(m2,:,rr+1) = sum(y_filter_stage,1)./2;
                    end
                end
            end
            y_multibeamformed_total(:,:,r_b1) = v_now;
        end
        runtime_fast_ds_2 = toc(init_fast_ds_2);
        Runtime_fast_ds_32_tap(jj,m_num) = (runtime_fast_ds_2 - counter_loop1 - counter_loop2)/N; % subtract filter init time


        %% runtime for delay and sum 16 tap
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
        y_ttd_in = (y_2d.*mod_phase).';
        % y_filter_ttd = TTD_Beamformer((y.*mod_phase).',-tau, fs, @Psi_sinc,h, h_idx, 'R',R).';
        Xf = fft(y_ttd_in,N_filter,1);
        Y = ifft(Hf.*Xf,[],1);
        y_filter_ttd = Y(1:N,:).';
        y_bf_ttd = reshape(sum(y_filter_ttd(:,1:N),1)/M, [N,1]);
        runtime_ds_beam = toc(t_init_beam_ds);
        Runtime_ds_16_tap(jj,m_num) = runtime_ds_beam*(R^2)/N;

        %% runtime for delay and sum 32 taps
        init_filter_ds2 = tic;
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
        runtime_filter_ds2 = toc(init_filter_ds2);

        init_beam_ds2 = tic;
        y_ttd_in2 = (y_2d.*mod_phase).';
        % y_filter_ttd = TTD_Beamformer((y.*mod_phase).',-tau, fs, @Psi_sinc,h, h_idx, 'R',R).';
        Xf2 = fft(y_ttd_in2,N_filter2,1);
        Y2 = ifft(Hf2.*Xf2,[],1);
        y_filter_ttd2 = Y2(1:N,:).';
        y_bf_ttd2 = reshape(sum(y_filter_ttd2(:,1:N),1)/M, [N,1]);
        runtime_beam_ds2 = toc(init_beam_ds2);
        Runtime_ds_32_tap(jj,m_num) = runtime_beam_ds2*(R^2)/N;
        
    end
end

figure(1)
semilogy(M_lst,mean(Runtime_fbst_superfast),'-*','LineWidth',1)
hold on
grid on
semilogy(M_lst, mean(Runtime_fbst_precompute),'-*','LineWidth',1);
hold on
grid on
semilogy(M_lst,mean(Runtime_ds_16_tap),'-^','LineWidth',1)
hold on
grid on
semilogy(M_lst, mean(Runtime_ds_32_tap),'-^','LineWidth',1)
hold on
grid on
semilogy(M_lst,mean(Runtime_fast_ds_16_tap),'-^','LineWidth',1)
hold on
grid on
semilogy(M_lst, mean(Runtime_fast_ds_32_tap),'-^','LineWidth',1)
% semilogy(M_lst, mean(Runtime_fast_ds_16_tap),'-^','LineWidth',1)
% hold on 
% grid on
% semilogy(M_lst, mean(Runtime_fast_ds_32_tap),'-^','LineWidth',1)
xlim([M_lst(1) M_lst(end)])
set(gca, 'XScale', 'log'); % Set x-axis to logarithmic scale
xticks(M_lst); % Set specific log2 tick marks
xticklabels({'2^2','2^4','2^6', '2^8', '2^{10}', '2^{12}', '2^{14}'}); % Label as powers of 2
xlabel('Array size (M)','Interpreter','latex')
ylabel('Runtime (s)','Interpreter','latex')
legend({'FBST Superfast','FBST Precompute','Delay and Sum (R=16)', 'Delay and Sum (R=32)', 'Fast Delay and Sum (R=16)', 'Fast Delay and Sum (R=32)'},'Location','northwest','Interpreter','latex','FontSize',12)
% exportgraphics(gcf, 'upa_runtime_vs_M.pdf', 'ContentType', 'vector');

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