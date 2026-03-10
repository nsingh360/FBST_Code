T_ap = 1.2;
M = 2^6;
N = 2^5;
MN = M*N;
L = ceil(N);

W = 10; % bandwidth
fs = 2.5*W; % sampling frequency
T = 1/fs;

% mc_iters = 100;
% Sigma = zeros(2*L,2*L);
% 
% for i=1:mc_iters
%     t = unifrnd(0,T_ap);
%     A_nu_basis = zeros(MN,2*L);
%         for kkk=1:2*L
%             f_k = (1/(2*T*L))*(kkk-1-L);
%             A_nu_basis(:,kkk) = exp(1i*2*pi*f_k*t);
%         end
%     Sigma = Sigma + (A_nu_basis.'*conj(A_nu_basis));
% end
% Sigma = Sigma/mc_iters;
% 
% D = eig(Sigma);
% plot(D);
% 
% t_1 = unifrnd(0,T_ap);
% A_nu_basis = zeros(MN,2*L);
% for kkk=1:2*L
%     f_k = (1/(2*T*L))*(kkk-1-L);
%     A_nu_basis(:,kkk) = exp(1i*2*pi*f_k*t_1);
% end
% Sigma_1 = (A_nu_basis.'*conj(A_nu_basis))/MN;
% 
% log(det(Sigma_1)/det(Sigma))


%% bounding Max using Fourier transform
bins = 50;
delta = 2*T_ap/bins;
t = linspace(-T_ap,T_ap, bins);
FT = zeros(bins,1);
freq = linspace(-4*W, 4*W, bins);
n_sinusoids = 100;

a = (randn(n_sinusoids,1) + 1i*randn(n_sinusoids,1))/sqrt(2); % sinusoid amplitudes
f = W*(2*rand(n_sinusoids,1)-1); % sinusoid frequencies

s = sig_gen(a,f,n_sinusoids,t);

for ff=1:bins
    sum_int = 0;
    for tt=1:bins
        sum_int = sum_int + s(tt)*exp(-1i*2*pi*freq(ff)*t(tt))*delta;
    end
    FT(ff,1) = sum_int;
end

%% Max over Bernstein ellipse
rho = (cot(pi/(4*T_ap)))^2;
theta_pts = linspace(-pi, pi, 1000);
berstein_boundary = 0.5*(rho*exp(1i*theta_pts) + (1/rho)*exp(-1i*theta_pts));

c_T = cos(pi/T_ap);
mapped_boundary = T_ap/pi*acos(c_T + 0.5*(1-c_T)*(berstein_boundary+1));
s_complex = sig_gen(a,f,n_sinusoids,mapped_boundary);

function [X] = sig_gen(a,f,n_sinusoids,t)
X = zeros(size(t));
for ii = 1:n_sinusoids % For each antenna
    sig = (a(ii)*exp(1i*2*pi*(f(ii))*(t)));
    X = X + sig;
end
end