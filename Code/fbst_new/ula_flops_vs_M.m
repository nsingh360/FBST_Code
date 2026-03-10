clear
clc
close all

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

M_lst = 2.^(1:14);
exp_fac = 1.005;
trials = 1;

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

    fc = 20e9;  % center frequency
    Ws = 5e9; % to prevent spacial aliasing
    fc_tilde = (fc+Ws); % adjustment for spatial aliasing
 
    fc_norm = fc/fc_tilde;
    c = physconst('LightSpeed');
    lambda = c/(fc+Ws); % wavelngth
    W = 5e9; % bandwidth
    fs = 2.5*W; % sampling frequency
    T = 1/fs;

    T_aperture_max = (M-1)*lambda/(2*c);
    N = max(32,ceil(1.1*T_aperture_max*fs));
    MN = M*N;
    L = ceil(N); % no. of samples in frequency is 2*L

    %% parameters for fast delay and sum
    S = log2(M);
    R_fac = 2; % 2-radix fast delay and sum
    L_fac = 2; % downsampling factor in spatial dimension
    eta = lambda/(2*c);

    freq_resolution = fs/N; % frequency grid for sub-band beamformig

    edge_lim = ceil(fs*T_aperture_max/2);
    N_eff = length(edge_lim+1:(N-(edge_lim+1)));
    fprintf("Array size: %d, Beams: %d, DFT effective snapshots: %d\n",M,2*B+1,N_eff);

    t_beams = 2*B+1;
    t_basis = 2*L+1;

    %% FLOPS for FBST (Chirp-z) step
    chirpz_time_domain = M*chirpz_flops(N,t_basis);
    chirpz_beam_domain = t_basis*chirpz_flops(M,t_beams);

    %% FLOPS for FBST (TOEPLITZ step)
    precompute_flops = t_beams*matvec_flops(t_basis,t_basis); % flops when applying the inverse directly
    canonical_vec_flops = t_beams*(2*fft_flops(t_basis,2*t_basis));

    %% FLOPS for FBST (reconstruction)
    reconstruc_flops = (t_beams)*chirpz_flops(t_basis,N);

    Runtime_fbst_precompute(mm1) = (chirpz_time_domain + chirpz_beam_domain + precompute_flops + reconstruc_flops)/N;
    Runtime_fbst_superfast(mm1) = (chirpz_time_domain + chirpz_beam_domain + canonical_vec_flops + reconstruc_flops)/N;

    %% FLOPS for sub-band processing (no zero padding)
    freq_fft_flops = M*fft_flops(N,N);
    chirpz_freq_beams_flops = N*chirpz_flops(M,t_beams);
    time_fft_flops = t_beams*fft_flops(N,N);

    Runtime_DFT_beamformer(mm1) = (freq_fft_flops + chirpz_freq_beams_flops + time_fft_flops)/N;


    %% FLOPS fast delay and sum
    Nf = N + 2*16+1;
    Runtime_fast_ds_16_tap(mm1) = M*log2(M)*(30*Nf*log2(Nf)+12*Nf+8*N)/N;

    Nf2 = N + 2*32+1;
    Runtime_fast_ds_32_tap(mm1) = M*log2(M)*(30*Nf2*log2(Nf2)+12*Nf2+8*N)/N;

    %% FLOPS for brute-force delay and sum
    Runtime_ds_16_tap(mm1) = t_beams*M*(30*Nf*log2(Nf)+12*Nf+8*N)/N;
    Runtime_ds_32_tap(mm1) = t_beams*M*(30*Nf2*log2(Nf2)+12*Nf2+8*N)/N;

    
end

figure(1)
semilogy(M_lst,Runtime_fbst_superfast,'-*','LineWidth',1)
hold on
grid on
% plot(B_lst, mean(Runtime_fbst2), 'LineWidth',1)
% hold on 
% grid on
% plot(B_lst, mean(Runtime_fbst_deconv),'-*','LineWidth',1)
% hold on
% grid on
semilogy(M_lst, Runtime_fbst_precompute,'-*','LineWidth',1);
hold on
grid on
% plot(B_lst, mean(Runtime_fbst_slepian_inv),'LineWidth',1)
% hold on
% grid on
semilogy(M_lst,Runtime_ds_16_tap,'-^','LineWidth',1)
hold on
grid on
semilogy(M_lst, Runtime_ds_32_tap,'-^','LineWidth',1)
hold on
grid on
semilogy(M_lst, Runtime_fast_ds_16_tap,'-^','LineWidth',1)
hold on 
grid on
semilogy(M_lst, Runtime_fast_ds_32_tap,'-^','LineWidth',1)
hold on
grid on
semilogy(M_lst,Runtime_DFT_beamformer,'-^','LineWidth',1)
xlim([M_lst(1) M_lst(end)])
set(gca, 'XScale', 'log'); % Set x-axis to logarithmic scale
xticks(M_lst); % Set specific log2 tick marks
xticklabels({'2^1','2^2','2^3','2^4','2^5', '2^6','2^7', '2^8','2^9','2^{10}','2^{11}','2^{12}','2^{13}','2^{14}'}); % Label as powers of 2
xlabel('Array size (M)','Interpreter','latex')
ylabel('Floating Point Operations (FLOPS)','Interpreter','latex')
legend({'FBST Superfast','FBST Precompute','Delay and Sum (R=16)', 'Delay and Sum (R=32)', 'Fast Delay and Sum (R=16)', 'Fast Delay and Sum (R=32)','Sub-Band Processing'},'Location','northwest','Interpreter','latex','FontSize',11)
exportgraphics(gcf, 'ula_flops_vs_M_subband.pdf', 'ContentType', 'vector');

%% supporting functions
function y = chirpz_flops(M,N)
% N : output length 
% M : input length 
L = N + M - 1;
y = 10*L*log2(L) + 6*(L+M+N);
end

function y = fft_flops(M,N)
% N : output length
% M : input length
y = 5*N*log2(N);
end



function y = matvec_flops(M,N)
% flops for multiplying MxN complex matrix with Nx1 complex vector
y = 8*M*N;
end