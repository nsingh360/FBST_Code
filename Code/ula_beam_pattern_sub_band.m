clear
clc
close all

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

B = 2^7; % no. of beams
M = 2^7; % array size
N = 2^6; % no. of samples
L = ceil(N); % no. of samples in frequency is 2*L
MN = M*N;

S = log2(M); % stages in the beamformer
R = 2; % 2-radix beamformer
L_fac = 2; % spatial downsampling factor

Theta_fbst = zeros(2*B+1,1);
idx = 240;%randi(B,1);
for i=1:2*B+1
    Theta_fbst(i) = asin((i-1-B)/B);
end

%% Angle sampling for final stage fast DS
r_lin = linspace(0,R^S-1,R^S);
Theta_ds = asin(1 - (2*r_lin + 1)/R^S);

fc = 20e9;  % center frequency
Ws = 5e9; % to prevent spatial aliasing
c = physconst('LightSpeed');
lambda = c/(fc+Ws); % wavelngth
eta = lambda/(2*c);
W = 5e9; % bandwidth
fs = 2*W; % sampling frequency
T = 1/fs;
idx_theta = 9;
theta = Theta_fbst(idx);%Theta(idx); 
m = [(-M/2+1/2):(M/2-1/2)]';
x_pos = m*lambda/2; % sensor positions
u_s = sin(theta); % normal vector for determining delays

freq_resolution = fs/N;

%% Signal specs that do not need to be redifined in loop
n_sinusoids = 100; % number of sinusoids in signal
n = [0:N-1]/fs; % temporal sample vectors
tau = x_pos*u_s/c; % relative delays to phase center
t_array = n-tau; % delays across array and time
t = t_array(:); %
T_aperture = max(t)-min(t);% temporal aperture of the array
demod_phase = repmat(exp(-1i*2*pi*fc*tau),N,1);%exp(1i*2*pi*fc*t);

taps = 8; % for ds
edge_lim = taps + ceil(fs*(max(tau)-min(tau))/2);% clip signal to account for edge distortion
comp_idxs = (edge_lim+1):(N-edge_lim-1); % comparison indicies

taps2 = 8; % for fast ds
edge_lim2 = taps2 + ceil(fs*(max(tau)-min(tau))/2);% clip signal to account for edge distortion
comp_idxs2 = (edge_lim2+1):(N-edge_lim2-1); % comparison indicies

%% construct basis matrix
frequency_samples = (1/(2*T*L))*linspace(-L,L-1,2*L);
F_source = exp(1i*2*pi*t*frequency_samples);
F_check = F_source;
% F_source = F_source.*demod_phase;
betas_no_null = inv(F_source'*F_source + 1e-5*eye(2*L))*F_source';

%% equispaced sampled basis matrix
t_nyq = [0:N-1]/fs; % nyquist samples for testing
psi = exp(1i*2*pi*t_nyq'*frequency_samples);
% for ll=1:2*L
%     f_l = (1/(2*T*L))*(ll-1-L);
%     psi(:,ll) = exp(1i*2*pi*f_l*t_nyq);
% end

%% final stage fast ds
delays2 = zeros(M,1);
for ll=1:2^S
    l = ll-1;
    r = 50;

    delay = 0;

    for k=0:S-1
        delay = delay + eta*(2^(k)*(floor(l/(2^(-k))) + 0.5) - M/2)*(2^(-k-1))*((-1)^(floor((2^(-S+1+k))*r)));
    end
    delays2(ll) = delay;
end

%% compute beampattern
W = 5e9; % bandwidth
f_beam = linspace(fc-0.5*W,fc+0.5*W,20);
theta_beam = linspace(-pi/2,pi/2,300);
theta_beam_zoomed = linspace(theta - 0.05,theta + 0.05, 50);
beam_pattern_subBand = zeros(length(f_beam),length(theta_beam));
e_mod = exp(-1i*2*pi*fc*t);

beam_pattern_subBand_zoomed = zeros(length(f_beam),length(theta_beam_zoomed));

for ii = 1:length(f_beam)
    for jj = 1:length(theta_beam)
        fprintf('Frequency: %.2f, Angle: %.2f\n',f_beam(ii),theta_beam(jj));
        tau_jj = x_pos*sin(theta_beam(jj))/c;
        t_jj = n-tau_jj;
        t_jj =t_jj(:);
        e_mod  = repmat(exp(-1i*2*pi*(fc)*(tau_jj)),N,1);
        e_jj = exp(1i*2*pi*(f_beam(ii)-fc)*t_jj).*e_mod;

        y = reshape(e_jj, [M,N]);

        % subband FFT across array elements 
        Y_f = fft(y.'); % mth column for mth array sample

        Y_nb = zeros(N,1);
        f_set = zeros(N,1);


        Y_freq_beamspace = zeros(N,2*B+1);
        beam_list = linspace(0,2*B,2*B+1);

        Y_beam = zeros(N,1);

        for nn=0:N-1
            ym = Y_f(nn+1,:).';
            if nn<N/2
                f_curr = fc + nn*freq_resolution;
            elseif nn>=N/2
                f_curr = fc + (nn)*freq_resolution - N*freq_resolution;
            end
            a_theta = exp(-1i*2*pi*f_curr*tau);
            Y_beam(nn+1) = a_theta'*ym;
            % W_czt = exp(1i*pi*f_curr/((fc+W)*B));
            % A_czt = exp(1i*pi*f_curr/(fc+W));
            % coeff2 = exp(1i*pi*f_curr*(M/2-1/2)/(fc+W));
            % coeff1 = exp(-1i*pi*f_curr*((M/2-1/2)/(fc+W))*(beam_list/B)).';
            % temp_var = czt(ym,2*B+1,W_czt,A_czt);
            % Y_freq_beamspace(nn+1,:) = coeff2*temp_var.*coeff1;
        end

        Y_time_beamspace = ifft(Y_beam);
        Y_bf = Y_time_beamspace/M; %

        beam_pattern_subBand(ii,jj) = norm(Y_bf)/sqrt(length(t_nyq));


        % %------- zoomed angles 
        % for jj = 1:length(theta_beam_zoomed)
        %     fprintf('Frequency: %.2f, Zoomed-Angle: %.2f\n',f_beam(ii),theta_beam_zoomed(jj));
        %     tau_jj = x_pos*sin(theta_beam_zoomed(jj))/c;
        %     t_jj = n-tau_jj;
        %     t_jj =t_jj(:);
        %     e_mod  = repmat(exp(-1i*2*pi*(fc)*(tau_jj)),N,1);
        %     e_jj = exp(1i*2*pi*(f_beam(ii)-fc)*t_jj).*e_mod;
        % 
        %     y = reshape(e_jj, [M,N]);
        % 
        %     % subband FFT across array elements 
        %     Y_f = fft(y.'); % mth column for mth array sample
        % 
        %     Y_nb = zeros(N,1);
        %     f_set = zeros(N,1);
        % 
        % 
        %     Y_freq_beamspace = zeros(N,2*B+1);
        %     beam_list = linspace(0,2*B,2*B+1);
        % 
        %     for nn=0:N-1
        %         ym = Y_f(nn+1,:).';
        %         if nn<N/2
        %             f_curr = fc + nn*freq_resolution;
        %         elseif nn>=N/2
        %             f_curr = fc + (nn)*freq_resolution - N*freq_resolution;
        %         end
        %         W_czt = exp(1i*pi*f_curr/((fc+W)*B));
        %         A_czt = exp(1i*pi*f_curr/(fc+W));
        %         coeff2 = exp(1i*pi*f_curr*(M/2-1/2)/(fc+W));
        %         coeff1 = exp(-1i*pi*f_curr*((M/2-1/2)/(fc+W))*(beam_list/B)).';
        %         temp_var = czt(ym,2*B+1,W_czt,A_czt);
        %         Y_freq_beamspace(nn+1,:) = coeff2*temp_var.*coeff1;
        %     end
        % 
        %     Y_time_beamspace = ifft(Y_freq_beamspace);
        %     Y_bf = Y_time_beamspace(:,idx)/M; %
        % 
        %     beam_pattern_subBand_zoomed(ii,jj) = norm(Y_bf)/sqrt(length(t_nyq));
        % end
    
    end
end



figure(1)
p0 = plot(rad2deg(theta_beam),zeros(length(theta_beam),1),'r--');
hold on
grid on
p2 = plot(rad2deg(theta_beam),db(beam_pattern_subBand),'Color',[0 0 1 0.6]);
xlabel('$\theta$ (degrees)','Interpreter','latex','FontSize',12)
ylabel('Response(dB)','Interpreter','latex','FontSize',12)

% legendHandles = [repmat(p0,50,1), p1, p2];
% legendLabels = {'Distortionless response','No Nulling', 'Projection Nulling'};
% lgd = legend(legendHandles, legendLabels);
% 
% legendChildren = findobj(lgd, 'Type', 'Line');
% set(legendChildren, 'LineWidth', 2); % Set solid lines for legend
% Dummy plots for solid legend entries with different styles
legend_p0 = plot(-1, -1, 'r--', 'LineWidth', 2);  % Solid red
legend_p1 = plot(-1, -1, 'k-', 'LineWidth', 2); % Dashed blue

% Legend
legend([legend_p0, legend_p1], {'Distortionless response', 'Sub band'},'location','northwest','Interpreter','latex','FontSize',12);
% legend({'Distortionless response','No Nulling', 'Projection Nulling'},'location','northwest','Interpreter','latex')
xlim([-90,90])
ylim([-45,.5])

% First Inset
ax1 = axes('Position', [0.55, 0.7, 0.2, 0.2]); % [x, y, width, height]
box on;
plot(rad2deg(theta_beam_zoomed(15:end-15)),zeros(length(theta_beam_zoomed(15:end-15)),1),'r--');
hold on;
plot(rad2deg(theta_beam_zoomed(15:end-15)), db(beam_pattern_subBand_zoomed(:,15:end-15)), 'Color',[0 0 0 1], 'LineWidth', 1);
y1 = ylim;
ylim([y1(1), 1.5]);
set(gca, 'YTick', []);
% exportgraphics(gcf, 'ula_beam_pattern_combined.pdf', 'ContentType', 'vector');

%% supporting functions
function [X] = sig_gen(a,f,n_sinusoids,t)
X = zeros(size(t));
for ii = 1:n_sinusoids % For each antenna
    sig = (a(ii)*exp(1i*2*pi*f(ii)*(t)));
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