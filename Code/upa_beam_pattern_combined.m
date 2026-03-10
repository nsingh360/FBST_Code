%%%% TODO: ASK IF THIS IS NECESSARY %%%
clear
clc
close all

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

S = 4; % no. of stages in the FDS
M = 2^(2*S); % sqrt(M) x sqrt(M) dimensional planar array
N = 2^6;
R_fac = 2; % R^(2*S) total beams
MN = M*N;

fc = 20e9;  % center frequency
Ws = 5e9;
fc_tilde = fc+Ws;
fc_norm = fc/fc_tilde;
c = physconst('LightSpeed');
lambda = c/(fc+Ws); % wavelngth
W = 5e9; % bandwidth
fs = 2.5*W; % sampling frequency
T = 1/fs;

%% Angle sampling for final stage
r_lin = linspace(0,R_fac^S-1,R_fac^S);
a1 = 1 - (2*r_lin + 1)/R_fac^S; % sin(theta)cos(phi)
b1 = a1; % sin(theta)sin(phi)

a_idx = 4;
b_idx = 7;
phi = atan(b1(b_idx)/a1(a_idx));
theta = asin(sqrt(a1(a_idx)^2 + b1(b_idx)^2));

m = (-sqrt(M)/2+1/2):(sqrt(M)/2-1/2);
[X_mesh, Y_mesh] = meshgrid(m,m);
x_pos = [X_mesh(:),Y_mesh(:)]*lambda/2; % sensor positions
u_s = sin(theta)*[cos(phi);sin(phi)]; % normal vector for determining delays

taps = 8; % for ds
edge_lim = taps + ceil(fs*(max(tau)-min(tau))/2);% clip signal to account for edge distortion
comp_idxs = (edge_lim+1):(N-edge_lim-1); % comparison indicies

taps2 = 16; % for fast ds
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

%% compute beampattern
W = 5e9; % bandwidth
f_beam = linspace(fc-0.5*W,fc+0.5*W,1);
theta_beam = linspace(-pi/2,pi/2,20);
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
        y_3d = reshape(e_jj, [sqrt(M),sqrt(M),N,1,1]); % sqrt(M)x sqrt(M) x N x R x R
        y_multibeamformed_first = zeros(sqrt(M),N,1,R_fac^S);
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
                        N_filter = N+2*taps2+1;
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
                        N_filter = N+2*taps2+1;
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

        v_pred = y_multibeamformed_total(:,b_idx,a_idx);
        
        % %% final stage fast ds
        % y_filter_fast_ttd2 = TTD_Beamformer(y.',-delays2, fs, @Psi_sinc, 'R',taps).';
        % y_bf_fast_ttd2 = sum(y_filter_fast_ttd2(:,1:N),1)/M;

        y_filter_ttd = TTD_Beamformer(y.',-tau, fs, @Psi_sinc, 'R',taps).';
        y_bf_ttd = sum(y_filter_ttd(:,1:N),1)/M;

        beam_pattern_fast_ds(ii,jj) = norm(v_pred)/sqrt(length(t_nyq));
        beam_pattern_ds(ii,jj) = norm(y_bf_ttd)/sqrt(length(t_nyq));
        beam_pattern_fbst(ii,jj) = norm(psi*(betas_no_null*e_jj))/sqrt(length(t_nyq));
        %         beam_pattern(ii,jj) = norm(V_nu*(Phi*e_jj))/norm(e_jj);%*norm(h_vec));
    end

    % %------- zoomed angles 
    % for jj = 1:length(theta_beam_zoomed)
    %     fprintf('Frequency: %.2f, Zoomed-Angle: %.2f\n',f_beam(ii),theta_beam_zoomed(jj));
    %     tau_jj = x_pos*sin(theta_beam_zoomed(jj))/c;
    %     t_jj = n-tau_jj;
    %     t_jj =t_jj(:);
    %     % e_mod  = repmat(exp(-1i*2*pi*(fc)*(tau_jj)),N,1);
    %     e_jj = exp(1i*2*pi*(f_beam(ii))*t_jj).*e_mod;
    % 
    %     %% fast ds beam pattern
    %     y = reshape(e_jj,[M,N,1]);
    %     v_now = y;
    %     fs_ds1 = tic;
    %     for stage=1:S
    %         v_prev = v_now;
    %         num_virtual_sensors = M/2^stage;
    %         num_sectors = R^stage;
    %         v_now = zeros(num_virtual_sensors,N,num_sectors); % init partial beamforming at stage s
    %         for mm=0:num_virtual_sensors-1
    %             for rr = 0:num_sectors-1
    %                 m = mm+1;
    %                 r = floor((rr)/R)+1;
    %                 delays = ((-1)^(rr))*(1/(R^stage))*eta*[2^(stage-1)*(2*mm+0.5)-M/2; 2^(stage-1)*(2*mm+1.5)-M/2];
    %                 y_interp = reshape(v_prev([2*mm+1,2*mm+2],:,r),[2,N]);
    %                 v_prev_interped = TTD_Beamformer(y_interp.',-delays, fs, @Psi_sinc, 'R',taps).';
    %                 v_now(m,:,rr+1) = sum(v_prev_interped,1)./2;
    %             end
    %         end
    %     end
        
    %     % %% final stage fast ds
    %     % y_filter_fast_ttd2 = TTD_Beamformer(y.',-delays2, fs, @Psi_sinc, 'R',taps).';
    %     % y_bf_fast_ttd2 = sum(y_filter_fast_ttd2(:,1:N),1)/M;
    % 
    %     y_filter_ttd = TTD_Beamformer(y.',-tau, fs, @Psi_sinc, 'R',taps).';
    %     y_bf_ttd = sum(y_filter_ttd(:,1:N),1)/M;
    % 
    %     beam_pattern_fast_ds_zoomed(ii,jj) = norm(v_now(1,:,idx_theta))/sqrt(length(t_nyq));
    %     beam_pattern_ds_zoomed(ii,jj) = norm(y_bf_ttd)/sqrt(length(t_nyq));
    %     beam_pattern_fbst_zoomed(ii,jj) = norm(psi*(betas_no_null*e_jj))/sqrt(length(t_nyq));
    %     %         beam_pattern(ii,jj) = norm(V_nu*(Phi*e_jj))/norm(e_jj);%*norm(h_vec));
    % end
end



figure(1)
p0 = plot(rad2deg(theta_beam),zeros(length(theta_beam),1),'r--');
hold on
p1 = plot(rad2deg(theta_beam),db(beam_pattern_fbst),'Color',[0 0 0 0.5]);
hold on
p2 = plot(rad2deg(theta_beam),db(beam_pattern_fast_ds),'Color',[0 0 1 0.5]);
hold on
p3 = plot(rad2deg(theta_beam),db(beam_pattern_ds),'Color',[0 1 0 0.5]);
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
legend([legend_p0, legend_p1, legend_p2, legend_p3], {'Distortionless response', 'FBST', 'FDS', 'DS'},'location','northwest','Interpreter','latex','FontSize',12);
% legend({'Distortionless response','No Nulling', 'Projection Nulling'},'location','northwest','Interpreter','latex')
xlim([-90,90])
ylim([-45,.5])

% First Inset
ax1 = axes('Position', [0.6, 0.6, 0.25, 0.25]); % [x, y, width, height]
box on;
plot(rad2deg(theta_beam_zoomed),zeros(length(theta_beam_zoomed),1),'r--');
hold on;
plot(rad2deg(theta_beam_zoomed), db(beam_pattern_fbst_zoomed), 'Color',[0 0 0 0.5], 'LineWidth', 1);
hold on
plot(rad2deg(theta_beam_zoomed), db(beam_pattern_fast_ds_zoomed), 'Color',[0 0 1 0.5],'LineWidth',1);
y1 = ylim;
ylim([y1(1), 0.5]);
% exportgraphics(gcf, 'beam_pattern_ula.pdf', 'ContentType', 'vector');

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