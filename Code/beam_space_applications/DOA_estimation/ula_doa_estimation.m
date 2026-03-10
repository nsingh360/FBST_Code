clear
clc
close all

%% set: Stress testing with large number of sources seems to work fine

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

B = 2^7; % no. of beams
M = 2^8; % array size
N = 2^6; % no. of samples
L = ceil(N); % no. of samples in frequency is 2*L
exp_fac = 1.005; % factor to include freqs slightly out of band;
MN = M*N;

I = 40; % no. of sources 

Theta = zeros(2*B+1,1);
idx = [10,50,100,200,250];%randi(B,1);
for i=1:2*B+1
    Theta(i) = asin((i-1-B)/B);
end
theta_off_grid = linspace(-pi/3,pi/3,I);

%% parameters for fast delay and sum
S = log2(M);
R_fac = 2; % 2-radix fast delay and sum
L_fac = 2; % downsampling factor in spatial dimension

fc = 20e9;  % center frequency
Ws = 5e9; % to prevent spacial aliasing
c = physconst('LightSpeed');
lambda = c/(fc+Ws); % wavelngth
eta = lambda/(2*c);
W = 5e9; % bandwidth
fs = 2*W; % sampling frequency
T = 1/fs;
theta = 56*pi/180;%Theta(idx);
m = [(-M/2+1/2):(M/2-1/2)]';
x_pos = m*lambda/2; % sensor positions
u_s = sin(theta); % normal vector for determining delays
fc_tilde = (fc+Ws); % adjustment for spatial aliasing
fc_norm = fc/fc_tilde;

hpbw = 2/(M*0.5);

t = zeros(MN,I);
demod_phase = zeros(MN,I);
T_aperture = zeros(I);
n_sinusoids = 100; % number of sinusoids in signal
n = [0:N-1]/fs; % temporal sample vectors
% V_NU = zeros(MN,K,I);
K=0;

freq_samples = exp_fac/(2*T*L)*linspace(-L,L,2*L+1);


Theta_rand = [-pi/3,-pi/6,0,pi/6, pi/3];%unifrnd(-pi/2+pi/10,pi/2-pi/10,I,1); % first index is source
for i=1:I

    theta = theta_off_grid(i);%iTheta_rand(i);
    u_s = sin(theta);
    tau = x_pos*u_s/c; % relative delays to phase center
    t_array = n-tau; % delays across array and time
    t(:,i) = t_array(:); %
    T_aperture(i) = max(t(:,i))-min(t(:,i));% temporal aperture of the array
    demod_phase(:,i) = repmat(exp(-1i*2*pi*(fc)*tau),N,1);

    if i==1
        K = ceil(2*W*T_aperture(i))+10;% subspace dimension
    end

    %% slepian function for signal generation
    [V_u,Lambda] = dpss(MN+1,W*T_aperture(i),K); % firt K uniformly sampled PSWF's
    t_u = [0:(MN)]'; % uniform sampling points over the interval
    t_nu = [M*N*(t(:,i)-min(t(:,i)))/T_aperture(i)]; % non-uniform sample points, scaled to match interval
    V_NU(:,:,i) = interp1(t_u,V_u,t_nu,'pchip',0); % non-uniform
    LAMBDA(:,i) = Lambda;
end


% toepltiz_invs = zeros(2*L+1,2*L+1,2*B+1);
% for bb=1:2*B+1
%     theta_b = Theta(bb);
%     u_sb = sin(theta_b);
%     taub = x_pos*u_sb/c; % relative delays to phase center
%     t_arrayb = n-taub; % delays across array and time
%     tb = t_arrayb(:); %
%     demod_phaseb = repmat(exp(-1i*2*pi*(fc)*taub),N,1);
%     Fb = exp(1i*2*pi*tb*freq_samples).*demod_phaseb;
%     toepltiz_invs(:,:,bb) = inv(Fb'*Fb +1e-5*eye(2*L+1));
% end



%% simulation specifications
trials = 10;
SNRs = -10;
SIR_dB = 0*ones(I-1,1);

doa_estimate = zeros(trials,2*B+1);
doa_estimate_ortho = zeros(trials, 2*B+1);

for ii = 1:trials
    for jj = 1:length(SNRs)

        SNR_dB = SNRs(jj); % signal to noise ratio
        % f = linspace(-W,W,n_sinusoids); % sinusoid frequencies
        fprintf('Trial: %d\n',ii)

        s_demod = zeros(MN,1);
        for i=1:I
            a = (randn(n_sinusoids,K) + 1i*randn(n_sinusoids,K))/sqrt(2).*(LAMBDA(:,i).'); % sinusoid amplitudes
            s = slepian_sig_gen(squeeze(V_NU(:,:,i)),a,n_sinusoids,t(:,i));
            if i==1 % source track
                s = s/rms(s);
            else 
                s = s/rms(s)*10^(-SIR_dB(i-1)/20);
            end
            s_demod = s_demod + s.*demod_phase(:,i);
        end

        sigma = 10^(-SNR_dB/20);


        noise = sigma*(randn(MN,1)+1i*randn(MN,1))/sqrt(2);

        y = s_demod + noise;
        y_vec = y;
        y = reshape(y, [M,N]);
        % F_transform = zeros(M,2*L);
        % A = -1;%exp(-1i*2*pi*T*W);
        % W_czt = exp(-1i*pi/L);
        % tic;
        % for mm=1:M
        %     F_transform(mm,:) = czt(y(mm,:),2*L,W_czt,A);
        % end
        % 
        % X_czt = zeros(2*L,B);
        % A= 1;
        % for ll=1:2*L
        %     l_dash = ll-1;
        %     % A = exp(-1i*pi*(1 + 2*(W/fc)*(l_dash-L/2)/L));
        %     W_czt = exp(1i*pi*(1/B)*(1 + (1/(2*fc*T*L))*(l_dash-L)));
        %     X_czt(ll,:) = czt(F_transform(:,ll),B,W_czt,A);
        % end
        % 
        % w_l = X_czt(:,idx);

        %% sanity check
        % y_check = zeros(M,N);
        % ls = reshape(linspace(-L,L-1,2*L),[2*L,1]);
        % w1 = -pi*(1 + (1/(2*fc*T*L))*ls)*sin(Theta(idx));
        % w2 = (pi/L)*ls;
        % F_L = ls*(1/2)*(1/(L*T));
        % for mm=1:M
        %     for nn=1:N
        %         y_check(mm,nn) = sum(w_l.*exp(1i*w1*mm).*exp(1i*w2*nn));
        %     end
        % end
        % w_l_check = zeros(2*L,1);
        % for ll=1:2*L
        %     sum_l = 0;
        %     for mm=1:M
        %         for nn=1:N
        %             sum_l = sum_l + y(mm,nn)*exp(-1i*w1(ll)*(m(mm)))*exp(-1i*w2(ll)*(nn-1));
        %         end
        %     end
        %     w_l_check(ll) = sum_l;
        % end
        % 
        % inner_prod_basis = zeros(2*L,MN);
        % for ll=1:2*L
        %     for mn=1:MN
        %         inner_prod_basis(ll,mn) = exp(-1i*2*pi*F_L(ll)*t(mn));
        %     end
        % end
        % w_l_check2 = inner_prod_basis*s;

        %% chirp-z transform
        F_transform = zeros(M,2*L+1);
        A = exp(-1i*pi*exp_fac);%exp(-1i*2*pi*T*W);
        W_czt = exp(-1i*pi*exp_fac/L);
        for mm=1:M
            F_transform(mm,:) = czt(y(mm,:).',2*L+1,W_czt,A);
        end
        X_czt = zeros(2*L+1,2*B+1);
        A= 1;
        B_list = reshape(linspace(0,2*B,2*B+1), [1,2*B+1]);
        for ll=1:2*L+1
            l_dash = ll-1;
            A = exp(1i*pi*(1/B)*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L))*B);
            W_czt = exp(1i*pi*(1/B)*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L)));
            coeff1 = exp(-1i*pi*(1/B)*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L))*(M/2-1/2)*B_list); % adjusting for array center
            coeff2 = exp(1i*pi*(1/B)*(fc_norm + (exp_fac/(2*fc_tilde*T*L))*(l_dash-L))*(M/2-1/2)*B); % adjusting for negative angles
            X_czt(ll,:) = czt(F_transform(:,ll),2*B+1,W_czt,A);
            X_czt(ll,:) = coeff2*X_czt(ll,:).*coeff1;
        end
        % 
        % for bb=1:2*B+1
        %     theta_b = Theta(bb);
        %     u_sb = sin(theta_b);
        %     taub = x_pos*u_sb/c; % relative delays to phase center
        %     t_arrayb = n-taub; % delays across array and time
        %     tb = t_arrayb(:); %
        %     demod_phaseb = repmat(exp(-1i*2*pi*(fc)*taub),N,1);
        %     Fb = exp(1i*2*pi*tb*freq_samples).*demod_phaseb;
        %     tinv = squeeze(toepltiz_invs(:,:,bb));
        %     doa_estimate_ortho(ii,bb) = norm(Fb*tinv*X_czt(:,bb));
        % end

        doa_var = vecnorm(X_czt,1);
        doa_estimate(ii,:) = doa_var./max(doa_var);

        % alpha_fbst = inv(A_eplitz + delta*eye(2*L))*F_basis'*y(:);
        % y_nyq = psi*alpha_fbst;

        % X_czt_ortho = zeros(2*L+1,2*B+1);
        % for bb=1:2*B+1
        %     theta = Theta(bb);
        %     u_s = sin(theta);
        %     tau = x_pos*u_s/c; % relative delays to phase center
        %     t_array = n-tau; % delays across array and time
        %     tt = t_array(:); %
        %     T_aperture = max(tt)-min(tt);% temporal aperture of the array
        %     demod_phase_f = repmat(exp(-1i*2*pi*(fc)*tau),N,1);
        %     F_basis = exp(1i*2*pi*tt*freq_samples).*demod_phase_f;
        %     X_czt_ortho(:,bb) = ((F_basis'*F_basis)^(-0.5))*X_czt(:,bb);
        % end
        % doa_var = vecnorm(X_czt_ortho,1);
        % doa_estimate_ortho(ii,:) = doa_var./max(doa_var);

    end
end

[pks,locs] = findpeaks(db(mean(doa_estimate)),'MinPeakProminence',2.5);


figure(1)
plot(Theta*180/pi,db(mean(doa_estimate)));
hold on
plot(Theta*180/pi,db(mean(doa_estimate_ortho)));
hold on
for i=1:length(locs)
    xline(Theta(locs(i))*180/pi,'r--', 'LineWidth', 1);
    hold on;
end
% exportgraphics(gcf, 'ula_snr_new.pdf', 'ContentType', 'vector');


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

function [X] = slepian_sig_gen(V,a,n_slepians,t)
X = sum(a*V.').';
end