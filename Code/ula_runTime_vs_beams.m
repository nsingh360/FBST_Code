clear
clc
close all

%% set random number generator seed
seed = 1;%randi([1,10000]);
rng(seed)

B_lst = 2^4:8:2^6;
trials = 25;
Runtime_fbst = zeros(trials, length(B_lst));
Runtime_fbst2 = zeros(trials, length(B_lst));
Runtime_fbst_cgd = zeros(trials, length(B_lst));
Runtime_fbst_slepian = zeros(trials, length(B_lst));
Runtime_fbst_slepian_inv = zeros(trials, length(B_lst));
% Runtime_fbst_deconv = zeros(trials, length(B_lst));
Runtime_fbst_superfast_toeplitz = zeros(trials, length(B_lst));
Runtime_ttd8 = zeros(trials, length(B_lst));
Runtime_ttd16 = zeros(trials, length(B_lst));

Runtime_chirpZ = zeros(trials, length(B_lst));
Runtime_sinusoid_inv = zeros(trials, length(B_lst));
Runtime_slepian_inv = zeros(trials, length(B_lst));
Runtime_sinusoid_cgd = zeros(trials, length(B_lst));
Runtime_sinusoid_toeptliz_superfast = zeros(trials, length(B_lst));

for b=1:length(B_lst)

    B = B_lst(b); % no. of beams
    M = 2^5; % array size
    N = 2^7; % no. of samples
    L = ceil(N); % no. of samples in frequency is 2*L
    MN = M*N;
    
    Theta = zeros(B,1);
    idx = B;%randi(B,1);
    for i=1:B
        Theta(i) = asin((i-1)/B);
    end
    
    fc = 20e9;  % center frequency
    c = physconst('LightSpeed');
    lambda = c/fc; % wavelngth
    W = 5e9; % bandwidth
    fs = 2.5*W; % sampling frequency
    T = 1/fs;
    theta = Theta(idx);
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
    % ls = reshape(linspace(-L,L-1,2*L),[2*L,1]);
    % F_L = (ls)*(1/2)*(1/(L*T));
    % A_col = sum(exp(1i*2*pi*(-1/(2*T) - F_L)*t'),2);
    % A_row = A_col';
    % A_eplitz = toeplitz(A_col, A_row);
    A_col = A_eplitz(:,1);
    A_row = A_eplitz(1,:);
    A_matMul = [A_col; 0; flipud(conj(A_row(2:end)'))];
    toeplitz_inv = (A_eplitz + delta*eye(2*L))^-1;

    %% get canonical vectors of Toeplitz inverse via brute force
    [X1,X2,X3,X4] = eval_canonical_vecs_brute_force(A_eplitz + delta*eye(2*L),2*L);

    %% pre-processing for superfast toeplitz
    c_x1 = X1(:, 1);
    r_x1 = X1(1, :);
    p_x1 = [c_x1; 0; flipud(conj(r_x1(2:end)'))];

    c_x2 = X2(:, 1);
    r_x2 = X2(1, :);
    p_x2 = [c_x2; 0; flipud(conj(r_x2(2:end)'))];

    c_x3 = X3(:, 1);
    r_x3 = X3(1, :);
    p_x3 = [c_x3; 0; flipud(conj(r_x3(2:end)'))];

    c_x4 = X4(:, 1);
    r_x4 = X4(1, :);
    p_x4 = [c_x4; 0; flipud(conj(r_x4(2:end)'))];

    %% construct slepian matrix
    [A_basis, V_nyq, K] = construct_fast_slepian_basis(W, t, T, L, M, N, T_aperture);
    Phi = ((A_basis'*A_basis + delta*eye(K))^-1)*(A_basis');
    
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

    %% convolution filter for fbst_deconv
    ls = reshape(linspace(-L,L-1,2*L),[2*L,1]);
    F_L = (ls)*(1/2)*(1/(L*T));
    h_conv = sum(exp(-1i*2*pi*F_L*t'),2);

    SNR = 10;
    for ii=1:trials
        SNR_dB = SNR; % signal to noise ratio
        a = (randn(n_sinusoids,1) + 1i*randn(n_sinusoids,1))/sqrt(2); % sinusoid amplitudes
        f = W*(2*rand(n_sinusoids,1)-1); % sinusoid frequencies
        s = sig_gen(a,f,n_sinusoids,t);
        sigma = (norm(s)/sqrt(MN))*10^(-SNR_dB/20);
        noise = sigma*(randn(MN,1)+1i*randn(MN,1))/sqrt(2);
        SNR_check = db(norm(s)/norm(noise));
        fprintf('Trial: %d, Set SNR: %.2f, Measured SNR: %.2f\n',ii,SNR_dB,SNR_check)

        s_demod = s.*demod_phase;

        y = s_demod + noise;
        y = reshape(y, [M,N]);
        
        %% runtime for fbst
        ts = tic;
        F_transform = zeros(M,2*L);
        A = -1;%exp(-1i*2*pi*T*W);
        W_czt = exp(-1i*pi/L);
        for mm=1:M
            F_transform(mm,:) = czt(y(mm,:),2*L,W_czt,A);
        end

        X_czt = zeros(2*L,B);
        A= 1;
        B_list = reshape(linspace(0,B-1,B), [1,B]);
        for ll=1:2*L
            l_dash = ll-1;
            % A = exp(-1i*pi*(1 + 2*(W/fc)*(l_dash-L/2)/L));
            W_czt = exp(1i*pi*(1/B)*(1 + (1/(2*fc*T*L))*(l_dash-L)));
            coeff = exp(-1i*pi*(1/B)*(1 + (1/(2*fc*T*L))*(l_dash-L))*(M/2-1/2)*B_list);
            X_czt(ll,:) = czt(F_transform(:,ll),B,W_czt,A);
            X_czt(ll,:) = X_czt(ll,:).*coeff;
        end
        
        w_l = X_czt(:,idx);
        t_inv = tic;
        for bb=1:B
            alpha = toeplitz_inv*w_l;
            y_nyq = psi*alpha;
        end
        inv_time = toc(t_inv);
        runtime_fbst = toc(ts);
        
        % t_inv2 = tic;
        % for bb=1:B
        %     alpha = inv(A_eplitz + delta*eye(2*L))*w_l;
        %     y_nyq = psi*alpha;
        % end
        % inv_time2 = toc(t_inv2);

        % t_inv3 = tic;
        % for bb=1:B
        %     [alpha_deconv,r] = deconv(w_l,h_conv,"same",Method="least-squares",RegularizationFactor=1e-5);
        %     y_nyq_deconv = psi*alpha_deconv;
        % end
        % inv_time3 = toc(t_inv3);

        t_inv4 = tic;
        for bb=1:B
            % Phi = ((A_basis'*A_basis + delta*eye(K))^-1)*(A_basis');
            alpha_slepian = Phi*w_l;
            y_nyq_slepian = V_nyq*alpha_slepian;
        end
        inv_time4 = toc(t_inv4);

        % t_inv5 = tic;
        % for bb=1:B
        %     Phi = ((A_basis'*A_basis + delta*eye(K))^-1)*(A_basis');
        %     alpha_slepian_inv = Phi*w_l;
        %     y_nyq_slepian_inv = V_nyq*alpha_slepian_inv;
        % end
        % inv_time5 = toc(t_inv5);

        %% CGD for fourier extension based fbst
        t_inv6 = tic;
        for bb=1:B
            alpha_cgd = zeros(2*L, 1);    % Initial guess
          
    
            % Use FFT for efficient Toeplitz matrix-vector multiplication
            matmul = A_eplitz*alpha_cgd;%ifft(fft(A_matMul) .* fft([alpha_cgd; zeros(2*L, 1)]));
            matmul = matmul(1:2*L);
            res = w_l - matmul;  % Initial residual
            p = res;              % Initial search direction
            rsold = res' * res;     % Initial residual norm squared
            
            for iter = 1:5
                matmul = ifft(fft(A_matMul) .* fft([p; zeros(2*L, 1)]));
                matmul = matmul(1:2*L);
                Ap = matmul;  % Efficient multiplication with Toeplitz matrix
                step = rsold / (p' * Ap);
                alpha_cgd = alpha_cgd + step * p;
                res = res - step * Ap;
                rsnew = res' * res;
                
                if sqrt(rsnew) < 1e-6
                    break;
                end
                
                p = res + (rsnew / rsold) * p;
                rsold = rsnew;
            end
            y_nyq_cgd = psi*alpha_cgd;
        end
        inv_time6 = toc(t_inv6);

        %% superfast toeplitz 
        t_inv7 = tic;
        for bb=1:B
            b1 = ifft(fft(p_x2) .* fft([w_l; zeros(length(w_l), 1)]));
            b1 = b1(1:length(w_l));
            % b1 = toep_mult(X2,w_l);
            b2 = ifft(fft(p_x4) .* fft([w_l; zeros(length(w_l), 1)]));
            b2 = b2(1:length(w_l));
            % b2 = toep_mult(X4,w_l);
            b3 = ifft(fft(p_x1) .* fft([b1; zeros(length(b1), 1)]));
            b3 = b3(1:length(b1));
            % b3 = toep_mult(X1,b1);
            b4 = ifft(fft(p_x3) .* fft([b2; zeros(length(b2), 1)]));
            b4 = b4(1:length(b2));
            % b4 = toep_mult(X3,b2);
            alpha_superfast = b3 - b4;
            y_nyq_superfast = psi*alpha_superfast;
        end
        inv_time7 = toc(t_inv7);

        % runtime_fbst2 = runtime_fbst - inv_time + inv_time2;
        % runtime_fbst3 = runtime_fbst - inv_time + inv_time3;
        runtime_fbst_cgd = runtime_fbst - inv_time + inv_time6;
        runtime_fbst_superfast = runtime_fbst - inv_time + inv_time7;
        runtime_slepian = runtime_fbst - inv_time + inv_time4;
        % runtime_slepian_inv = runtime_fbst - inv_time + inv_time5;

        runtime_chirpz = runtime_fbst - inv_time;
        runtime_sinusoid_inv = inv_time;
        runtime_slepian_inv = inv_time4;
        runtime_sinusoid_cgd = inv_time6;
        runtime_sinusoid_superfast = inv_time7;

        s_nyq = sig_gen(a,f,n_sinusoids,t_nyq);
        s_nyq_2 = sig_gen(a,f,n_sinusoids,t_nyq+1/fs);

        %% runtime for ttd with 16 taps
        ts2 = tic;
        y_ttd_in = (y.*mod_phase).';
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
        for bb=1:B
            % y_filter_ttd = TTD_Beamformer((y.*mod_phase).',-tau, fs, @Psi_sinc,h, h_idx, 'R',R).';
            Xf = fft(y_ttd_in,N_filter,1);
            Y = ifft(Hf.*Xf,[],1);
            y_filter_ttd = Y(1:N,:).';
            y_bf_ttd = reshape(sum(y_filter_ttd(:,1:N),1)/M, [N,1]);
        end
        runtime_ttd = toc(ts2);

        %% runtime for ttd with 32 taps
        ts3 = tic;
        y_ttd_in2 = (y.*mod_phase).';
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
        for bb=1:B
            % y_filter_ttd = TTD_Beamformer((y.*mod_phase).',-tau, fs, @Psi_sinc,h, h_idx, 'R',R).';
            Xf2 = fft(y_ttd_in2,N_filter2,1);
            Y2 = ifft(Hf2.*Xf2,[],1);
            y_filter_ttd2 = Y2(1:N,:).';
            y_bf_ttd2 = reshape(sum(y_filter_ttd2(:,1:N),1)/M, [N,1]);
        end
        runtime_ttd16 = toc(ts3);

        Runtime_fbst(ii,b) = runtime_fbst;
        % Runtime_fbst2(ii,b) = runtime_fbst2;
        % Runtime_fbst_deconv(ii,b) = runtime_fbst3;
        Runtime_fbst_cgd(ii,b) = runtime_fbst_cgd;
        Runtime_fbst_slepian(ii,b) = runtime_slepian;
        % Runtime_fbst_slepian_inv(ii,b) = runtime_slepian_inv;
        Runtime_fbst_superfast_toeplitz(ii,b) = runtime_fbst_superfast;
        Runtime_ttd8(ii,b) = runtime_ttd;
        Runtime_ttd16(ii,b) = runtime_ttd16;

        Runtime_chirpZ(ii,b) = runtime_chirpz;
        Runtime_sinusoid_inv(ii,b) = runtime_sinusoid_inv;
        Runtime_slepian_inv(ii,b) = runtime_slepian_inv;
        Runtime_sinusoid_cgd(ii,b) = runtime_sinusoid_cgd;
        Runtime_sinusoid_toeptliz_superfast(ii,b) = runtime_sinusoid_superfast;
    end

end

figure(1)
plot(B_lst,mean(Runtime_fbst),'-*','LineWidth',1)
hold on
grid on
% plot(B_lst, mean(Runtime_fbst2), 'LineWidth',1)
% hold on 
% grid on
% plot(B_lst, mean(Runtime_fbst_deconv),'-*','LineWidth',1)
% hold on
% grid on
plot(B_lst, mean(Runtime_fbst_cgd),'-*','LineWidth',1)
hold on
grid on
plot(B_lst, mean(Runtime_fbst_superfast_toeplitz),'-*','LineWidth',1);
hold on
grid on
plot(B_lst, mean(Runtime_fbst_slepian),'-*','LineWidth',1)
hold on
grid on
% plot(B_lst, mean(Runtime_fbst_slepian_inv),'LineWidth',1)
% hold on
% grid on
plot(B_lst,mean(Runtime_ttd8),'-^','LineWidth',1)
hold on
grid on
plot(B_lst, mean(Runtime_ttd16),'-^','LineWidth',1)
xlabel('no. of beams (B)','Interpreter','latex')
ylabel('runtime (s)','Interpreter','latex')
legend({'FBST (Sinusoid basis)','FBST (Sinusoid basis) CGD', 'FBST (Sinusoid basis) Superfast toeplitz','FBST (Slepian basis)','Delay and Sum (R=16)', 'Delay and Sum (R=32)'},'Location','northwest','Interpreter','latex')


figure(2);
plot(B_lst, mean(Runtime_chirpZ),'LineWidth',1);
hold on
grid on
plot(B_lst, mean(Runtime_sinusoid_inv),'LineWidth',1)
hold on 
grid on
plot(B_lst, mean(Runtime_sinusoid_cgd),'LineWidth',1)
hold on
grid on
plot(B_lst, mean(Runtime_sinusoid_toeptliz_superfast),'LineWidth',1)
hold on 
grid on 
plot(B_lst, mean(Runtime_slepian_inv),'LineWidth',1)
xlabel('no. of beams (B)','Interpreter','latex')
ylabel('runtime (s)','Interpreter','latex')
legend({'Chirp Z Step','FBST (Sinusoid basis) INV Step', 'FBST (Sinusoid basis) CGD Step','FBST (Sinusoid basis) Superfast Toepltiz Step', 'FBST (Slepian basis) INV Step'},'Location','northwest','Interpreter','latex')



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

function y = toep_mult(T, x)
    % Efficient multiplication of Toeplitz matrix T with vector x
    %
    % T: Toeplitz matrix
    % x: vector to be multiplied
    
    % Get the first column and first row of the Toeplitz matrix
    c = T(:, 1);
    r = T(1, :);
    
    % Use FFT for efficient Toeplitz matrix-vector multiplication
    n = length(x);
    p = [c; 0; flipud(conj(r(2:end)'))];
    y = ifft(fft(p) .* fft([x; zeros(n, 1)]));
    y = y(1:n);
end