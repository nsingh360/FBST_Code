clear
clc
close all

%% set random number generator seed
seed = randi([1,10000]);
rng(seed)

%% array spatial and temporal specificiations
M = 2^6; % array size
N = 2^5; % temporal samples
MN = M*N;
fc = 20e9;  % center frequency
c = physconst('LightSpeed');
W = 5e9; % bandwidth
lambda = c/(fc+W); % wavelngth
% lambda = c/(fc); % wavelngth
fs = 2*W; % sampling frequency
phi = pi/4; % azimuth
theta = pi/3; % elevation
m = [(-M/2+1/2):(M/2-1/2)]';
x_pos = m*lambda/2; % sensor positions
u_s = sin(theta); % normal vector for determining delays


%% Signal specs that do not need to be redifined in loop
n_sinusoids = 100; % number of sinusoids in signal
n = [0:N-1]/fs; % temporal sample vectors
tau = x_pos*u_s/c; % relative delays to phase center
t_array = n-tau; % delays across array and time
t = t_array(:); %
T = max(t)-min(t);% temporal aperture of the array
L = 2;
K = ceil(2*W*T)+L;% subspace dimension

%% Interferer specs that do not need to be redfined in loop
phi_i = -pi/4;
theta_i = pi/6;
u_i = sin(theta_i); % normal vector for determining delays
tau_i = x_pos*u_i/c;
t_i_array = n - tau_i;
t_i = t_i_array(:);
T_i = max(t_i) - min(t_i);
e_i  = repmat(exp(1i*2*pi*fc*(-tau_i+tau)),N,1); % modulation vector after pre-steering


%% PSWF interpolation
OS_fac = 10;
t_nyq = [0:N-1]'/fs; % nyquist samples for testing

V_u = dpss(OS_fac*MN+1,(W*T),K); % firt K uniformly sampled PSWF's
t_u = [0:(OS_fac*MN)]'; % uniform sampling points over the interval
t_nu = [OS_fac*MN*(t-min(t(:)))/T]; % non-uniform sample points, scaled to match interval
V_nu = interp1(t_u,V_u,t_nu,'pchip',0); % non-uniform
V_nu = V_nu;

V_nyq = interp1(t_u,V_u,OS_fac*MN*(t_nyq-min(t))/T,'pchip',0); % interpolate to nyquist


OS_fac = 1;
V2_u = dpss(OS_fac*MN+1,W*T,K); % firt K uniformly sampled PSWF's
t2_u = [0:(OS_fac*MN)]'; % uniform sampling points over the interval
t2_nu = [OS_fac*MN*(t-min(t(:)))/T]; % non-uniform sample points, scaled to match interval
V2_nu = interp1(t2_u,V2_u,t2_nu,'spline',0); % non-uniform



V2_nyq = interp1(t2_u,V2_u,OS_fac*MN*(t_nyq-min(t))/T,'pchip',0); % interpolate to nyquist



%% simulation specifications
trials = 50;
SIR_dB = -60;%-30:6:30;
SNR_dB = 0;%:6:30;
samples = 10*100000;

sigma_n = (1/sqrt(MN))*10^(-SNR_dB/20);
sigma_i = 10^(-SIR_dB/20);




cutoff = 1e-8;


V_orth = orth(V_nu);

n_sig = [-2*N:4*N]/(2*W);

B = sinc(2*W*(t-n_sig));
B_scale = 1/sqrt(trace(B'*B));
B = B*B_scale;

Bi = bsxfun(@times,e_i,sinc(2*W*(t_i-n_sig)));
Bi_scale = 1/sqrt(trace(Bi'*Bi));
Bi = Bi*Bi_scale;

B_nyq = sinc(2*W*(t_nyq-n_sig))*B_scale;


Rn = (sigma_n^2)*eye(MN);
Rs = (B*B');
Ri = (sigma_i^2)*(Bi*Bi');

R = Rs + Rn + Ri;

R_inv = R^-1;
UR_invU = V_nu'*(R_inv*V_nu);
Phi_full = (UR_invU^-1)*(V_nu'*R_inv);
V_pinv = ((V_nu'*V_nu)^-1)*V_nu';
% Psi = Phi_full*V_nu;
% [U,~,~] = svd(Psi);
% Psi_null = U(:,K+1:end);


% Psi_inv = ((Psi'*Psi)^-1)*Psi';
% 
% R_null = Psi_null*((Psi_null'*R_phi*Psi_null)^-1)*Psi_null';
% R_null = (R_null + R_null)/2;
% C = (R_null*R_phi)*Psi_inv';
% Phi_tilde = Psi_inv - C';


var_nom = real(trace((UR_invU^-1))-trace(V_pinv*Rs*V_pinv'));
Rs_nom = V_pinv*Rs*V_pinv';


V = orth(R_inv*V_nu);

D = 100;
bits = 1;1:12;
num_meas = zeros(1,length(bits));
thresh = 1e-9;



for jj = 1:length(bits)
    phi = zeros(MN,10*K);
phi_orth = zeros(MN,10*K);
P = 1;
sin_prev = inf;
    for ii = 1:100*K
%         if ii == 1
%             phi_orth = zeros(M*N,1);
%         else
%             phi_orth = orth(phi(:,1:ii-1));
%         end
        
        U = (randn(K,D) + 1i*randn(K,D))/sqrt(2);
        U = U;%./sqrt(sum(abs(U).^2,1));
        Q = V*U - phi_orth(:,1:(P))*(phi_orth(:,1:(P))'*V)*U;%-phi_orth*(phi_orth'*V*U);
%         Q = Q./max(abs(Q));
        kappa = max(sqrt(sum(abs(V).^2,2)))/sqrt(2);
%         Q = (Q)./abs(Q);
        Q = (quantizer_fixed(Q,bits(jj),1,0));
        
        Q_scale = Q./sqrt(sum(abs(Q).^2,1));
        
        [H,S,W] = svd(phi_orth(:,1:(P))'*V,'econ');
        sig_max = max(max(diag(S)));
        sig_min = min(min(diag(S)));
        Ht = H(:,diag(S)>thresh*sig_max);
        Wt = W(:,diag(S)>thresh*sig_max);
        St = S(diag(S)>thresh*sig_max,diag(S)>thresh*sig_max);
        Vi = V*Wt;
%         G = ((Vi'*Q_scale));
%         phi_hat = 
%         G = V'*Q_scale - V'*phi_orth(:,1:(P-1))*(phi_orth(:,1:(P-1))'*Q_scale);
        G = (V'*Q_scale);
%         G = V'*Q_scale - Wt*St*(Vi'*Q_scale);
        G_norm = (sqrt(sum(abs(G).^2,1)));
        
%         G_norm = G_norm.^2 + (1-min(min(diag(S))))*(sum(abs(V'*Q_scale).^2,1));
        [~,idx] = max(G_norm);
        
        if P>K
            phi_ii = [phi(:,1:(P-1)), Q(:, idx(1))];
            phi_ii_orth = orth(phi_ii);
            [~,S,~] = svd(phi_ii_orth'*V,'econ');
            Re = (((phi_ii'*V_nu)'*(((phi_ii'*R)*phi_ii)^-1)*(phi_ii'*V_nu)))^-1;
            var_ii = trace(Re)-trace(Rs_nom);
            sin_ii = real(norm(acos(diag(S)))^2);
            
            if db(var_ii/var_nom)/2<.1%.8*K*pi/2
                fprintf('Number of %d-bit measurements = %d\n',bits(jj),P-1)
                phi = phi_ii(:,1:(P));
                num_meas(jj) = P;
                break
            end
            
            if (sin_ii<.98*sin_prev)      
                phi(:,P) = Q(:, idx(1));
                phi_orth(:,P) = phi(:,P) - phi_orth(:,1:P-1)*(phi_orth(:,1:P-1)'*phi(:,P));
                phi_orth(:,P) = phi_orth(:,P)/norm(phi_orth(:,P));
                P = P+1;
                fprintf('%d-bit, iteration = %d, Number measurements = %d, sin = %.2f, variance = %.3f\n',bits(jj),ii, P, sin_ii,db(var_ii/var_nom))
                sin_prev = sin_ii;
            end
            
            
            
        else
            phi(:,P) = Q(:, idx(1));
            if ii==1
                phi_orth(:,P) = phi(:,P)/norm(phi(:,P));
            else
                phi_orth(:,P) = phi(:,P) - phi_orth(:,1:P-1)*(phi_orth(:,1:P-1)'*phi(:,P));
                phi_orth(:,P) = phi_orth(:,P)/norm(phi_orth(:,P));
            end
            P = P+1;
        end
        
    end
end




Phi = phi';

% Phi = ([R_inv*V_nu, randn(MN,10)]*randn(K+10))';
% Phi = ([R_inv*V_nu])';


Psi = Phi*V_nu;
[U,~,~] = svd(Psi);
Psi_null = U(:,K+1:end);
R_phi = (Phi*(R*Phi'));
R_phi = (R_phi + R_phi')/2;
R_phi_inv = R_phi^-1;

Psi_inv = ((Psi'*(R_phi_inv)*Psi + 1e-5*eye(K))^-1)*Psi'*R_phi_inv;
Phi_tilde = Psi_inv;


% R_phi_inv = R_phi^-1;
% Phi_tilde = ((Psi'*(R_phi_inv*Psi))^-1)*(Psi'*R_phi_inv);



% Phi = quantizer(V_nu',bits);
% Psi = Phi*V_nu;
% Phi_tilde = ((Psi'*Psi)^-1)*Psi';




%% compute beampattern
W = 5e9; % bandwidth
f_beam = linspace(-W,1*W,50)*.5;
theta_beam = linspace(-pi/2,pi/2,2000);
beam_pattern = zeros(length(f_beam),length(theta_beam));
e_mod = exp(-1i*2*pi*fc*t);

for ii = 1:length(f_beam)
    ii
    for jj = 1:length(theta_beam)
        tau_jj = x_pos*sin(theta_beam(jj))/c;
        t_jj = n-tau_jj;
        t_jj =t_jj(:);
        %         e_mod  = repmat(exp(1i*2*pi*(fc+f_beam(ii))*(tau_jj)),2*K,1);
        e_jj = exp(1i*2*pi*(fc + f_beam(ii))*t_jj).*e_mod;
        beam_pattern(ii,jj) = norm(V_nyq*(Phi_tilde*(Phi*e_jj)))/sqrt(length(t_nyq));
        %         beam_pattern(ii,jj) = norm(V_nu*(Phi*e_jj))/norm(e_jj);%*norm(h_vec));
    end
end

save('slepian_beam.mat','beam_pattern','theta_beam')

figure(1)
plot(rad2deg(theta_beam),zeros(length(theta_beam),1),'r--')
hold on
plot(rad2deg(theta_beam),db(beam_pattern),'k')
grid on
xlabel('$\theta$ (degrees)','Interpreter','latex')
ylabel('Response(dB)','Interpreter','latex')
legend({'Distortionless response'},'location','northwest','Interpreter','latex')
xlim([-90,90])
ylim([-45,.5])
% saveas(gca,'figs/slepian_mvdr_ula_beampattern_1_bit.png')





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