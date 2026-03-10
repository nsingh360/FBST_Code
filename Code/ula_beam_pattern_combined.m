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
theta = Theta_ds(idx_theta);%Theta(idx); 
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
beam_pattern_fbst = zeros(length(f_beam),length(theta_beam));
beam_pattern_ds = zeros(length(f_beam),length(theta_beam));
beam_pattern_fast_ds = zeros(length(f_beam),length(theta_beam));
e_mod = exp(-1i*2*pi*fc*t);

beam_pattern_fbst_zoomed = zeros(length(f_beam),length(theta_beam_zoomed));
beam_pattern_ds_zoomed = zeros(length(f_beam),length(theta_beam_zoomed));
beam_pattern_fast_ds_zoomed = zeros(length(f_beam),length(theta_beam_zoomed));

for ii = 1:length(f_beam)
    for jj = 1:length(theta_beam)
        fprintf('Frequency: %.2f, Angle: %.2f\n',f_beam(ii),theta_beam(jj));
        tau_jj = x_pos*sin(theta_beam(jj))/c;
        t_jj = n-tau_jj;
        t_jj =t_jj(:);
        % e_mod  = repmat(exp(-1i*2*pi*(fc)*(tau_jj)),N,1);
        e_jj = exp(1i*2*pi*(f_beam(ii))*t_jj).*e_mod;

        %% fast ds beam pattern
        y = reshape(e_jj,[M,N,1]);
        v_now = y;
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
        v_pred = v_now(1,:,idx_theta);

        %% final stage fast ds
        y_filter_fast_ttd2 = TTD_Beamformer(y.',-delays2, fs, @Psi_sinc, 'R',taps).';
        y_bf_fast_ttd2 = sum(y_filter_fast_ttd2(:,1:N),1)/M;

        y_filter_ttd = TTD_Beamformer(y.',-tau, fs, @Psi_sinc, 'R',taps).';
        y_bf_ttd = sum(y_filter_ttd(:,1:N),1)/M;

        beam_pattern_fast_ds(ii,jj) = norm(v_pred(comp_idxs2))/sqrt(length(comp_idxs2));
        beam_pattern_ds(ii,jj) = norm(y_bf_ttd(comp_idxs))/sqrt(length(comp_idxs));
        beam_pattern_fbst(ii,jj) = norm(psi*(betas_no_null*e_jj))/sqrt(length(t_nyq));
                % beam_pattern(ii,jj) = norm(V_nu*(Phi*e_jj))/norm(e_jj);%*norm(h_vec));
    end

    %------- zoomed angles 
    for jj = 1:length(theta_beam_zoomed)
        fprintf('Frequency: %.2f, Zoomed-Angle: %.2f\n',f_beam(ii),theta_beam_zoomed(jj));
        tau_jj = x_pos*sin(theta_beam_zoomed(jj))/c;
        t_jj = n-tau_jj;
        t_jj =t_jj(:);
        % e_mod  = repmat(exp(-1i*2*pi*(fc)*(tau_jj)),N,1);
        e_jj = exp(1i*2*pi*(f_beam(ii))*t_jj).*e_mod;

        %% fast ds beam pattern
        y = reshape(e_jj,[M,N,1]);
        v_now = y;
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
        v_pred = v_now(1,:,idx_theta);
        %% final stage fast ds
        y_filter_fast_ttd2 = TTD_Beamformer(y.',-delays2, fs, @Psi_sinc, 'R',taps).';
        y_bf_fast_ttd2 = sum(y_filter_fast_ttd2(:,1:N),1)/M;

        y_filter_ttd = TTD_Beamformer(y.',-tau, fs, @Psi_sinc, 'R',taps).';
        y_bf_ttd = sum(y_filter_ttd(:,1:N),1)/M;

        beam_pattern_fast_ds_zoomed(ii,jj) = norm(v_pred(comp_idxs2))/sqrt(length(comp_idxs2));
        beam_pattern_ds_zoomed(ii,jj) = norm(y_bf_ttd(comp_idxs))/sqrt(length(comp_idxs));
        beam_pattern_fbst_zoomed(ii,jj) = norm(psi*(betas_no_null*e_jj))/sqrt(length(t_nyq));
                % beam_pattern(ii,jj) = norm(V_nu*(Phi*e_jj))/norm(e_jj);%*norm(h_vec));
    end
end



figure(1)
p0 = plot(rad2deg(theta_beam),zeros(length(theta_beam),1),'r--');
hold on
grid on
p2 = plot(rad2deg(theta_beam),db(beam_pattern_fast_ds),'Color',[0 0 1 0.6]);
hold on
grid on
p3 = plot(rad2deg(theta_beam),db(beam_pattern_ds),'Color',[0 1 0 0.2]);
hold on
grid on
p1 = plot(rad2deg(theta_beam),db(beam_pattern_fbst),'Color',[0 0 0 0.5]);
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
legend_p2 = plot(-1, -1, 'b-', 'LineWidth', 2); % Dashed blue
legend_p3 = plot(-1, -1, 'g-', 'LineWidth', 2); % Dashed blue

% Legend
legend([legend_p0, legend_p1, legend_p2, legend_p3], {'Distortionless response', 'FBST', 'FDS (R=16)', 'DS (R=16)'},'location','northwest','Interpreter','latex','FontSize',12);
% legend({'Distortionless response','No Nulling', 'Projection Nulling'},'location','northwest','Interpreter','latex')
xlim([-90,90])
ylim([-45,.5])

% First Inset
ax1 = axes('Position', [0.55, 0.7, 0.2, 0.2]); % [x, y, width, height]
box on;
plot(rad2deg(theta_beam_zoomed(15:end-15)),zeros(length(theta_beam_zoomed(15:end-15)),1),'r--');
hold on;
plot(rad2deg(theta_beam_zoomed(15:end-15)), db(beam_pattern_fbst_zoomed(:,15:end-15)), 'Color',[0 0 0 1], 'LineWidth', 1);
hold on
plot(rad2deg(theta_beam_zoomed(15:end-15)), db(beam_pattern_fast_ds_zoomed(:,15:end-15)), 'Color',[0 0 1 0.2],'LineWidth',1);
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