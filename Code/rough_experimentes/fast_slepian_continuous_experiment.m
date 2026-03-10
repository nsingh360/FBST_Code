clc;
clear;
close all;

T = 2e-9; % signal defined over [0, T]
W = 5e9; % bandwidth of signal
N = 100; % sample to approximate the inner product over
Ts = T/N; 

%% create time samples
t = zeros(N,1);
for i=1:N
    t(i) = (i-1)*Ts;
end

%% create signal
n_sinusoids = 100;
a = (randn(n_sinusoids,1) + 1i*randn(n_sinusoids,1))/sqrt(2); % sinusoid amplitudes
f = W*(2*rand(n_sinusoids,1)-1); % sinusoid frequencies
s = sig_gen(a,f,n_sinusoids,t);

% %% construct (slow) slepian approximation using discrete least squares
% delta = 0;
% K = ceil(2*W*T)+6; % slepian subspace dimension
% if K > N
%     error("K should be less than N");
% end
% 
% slepian = dpss(N,W*T,K);
% Phi = ((slepian'*slepian + delta*eye(K))^-1)*(slepian');
% alpha_discrete = Phi*s;
% 
% s_estimate_discrete = slepian*alpha_discrete;

%% construct (slow) slepian approximation using continuous least squares 
%% NOTE : standard inner product is assumed over the basis
K = ceil(2*W*T)+6;
slepian = dpss(N,W*T,K);

b = zeros(K,1);
A = zeros(K,K);

for k1=1:K
    for k2=1:K
        b(k1) = slepian(:,k1)'*s;
        A(k2,k1) = slepian(:,k1)'*slepian(:,k2);
    end
end

alpha_slow = inv(A)*b;
s_estimate_slow = slepian*alpha_slow;

%% construct (fast) slepian approximation 
%% Error on individual terms (somewhat close to optimal)
epsilon = 1;
epsilonF = 2*epsilon/3; % Error in approximating B_{N,W}-F_{N,W}F_{N,W}^*
mu = 2*log(8*N-4)/pi^2+1/log(2)+1/log(6/pi);
deltah = (log(8*N-4)/pi)*epsilonF/mu;
deltaa = (1/log(2))*epsilonF/mu;
deltab = (1/log(6/pi))*epsilonF/mu;
epsilonS = epsilon/3; % Error in approximating B_{N,W}-S_KS_K^* 

%% Diagonal matrices (store as column vectors)
w = W*Ts;
R = round(w*(N)-1/2);
w_ = (R+1/2)/N;
vecn = ((-(N-1)/2):((N-1)/2))';
DAcos = cos((2*pi*w_)*vecn);
DAsin = sin((2*pi*w_)*vecn);
DBcos = cos((pi*(w+w_))*vecn);
DBsin = sin((pi*(w+w_))*vecn);

%% Approximate Scaled Hilbert Matrix (1/pi)*H ~~ Z*Z'
rh = ceil(log(8*N-4)*log(4*pi/deltah)/pi^2);
rh=2;

k1 = 1/(2*N-1);
k = sqrt(1-k1^2);
Kk = ellipticK(k^2);
[~,~,p] = ellipj((1:2:(2*rh-1))*Kk/(2*rh),k^2);
p = (N-1/2)*p;

Z = zeros(N,rh);
vecn = ((1/2):(N-1/2))';
Z(:,1) = sqrt(2*p(1)/pi)*ones(N,1)./(vecn+p(1));
for j = 2:rh
    Z(:,j) = sqrt(p(j)/p(j-1))*(1-(p(j)+p(j-1))./(vecn+p(j))).*(Z(:,j-1));
end

%% Approximate A1 ~~ Va*Ca*Va' and B0 ~~ Vb*Cb*Vb'
ra = ceil(log(2/(3*pi*deltaa))/(2*log(2)));
while((2/(3*pi))*sqrt(2/((4*ra-1)*4*ra))/2^(2*ra-2) <= deltaa)
    ra = ra-1;
end
ra=1;

vecZeta = zeta(2:2:(2*ra));
vecCoeff = (2/((N)*pi))*(1-(1-2.^(-(1:2:(2*ra-1)))).*vecZeta);
vecGamma = gamma(1:(2*ra));
Ca = zeros(2*ra,2*ra);
for k = 1:ra
    for l = 0:(2*k-1)
        Ca(l+1,2*k-l) = (-1)^(2*k-1-l)*vecCoeff(k).*vecGamma(2*k)./(vecGamma(l+1)*vecGamma(2*k-l));
    end
end

rb = 0;
while((2/pi)*sqrt(2/((4*rb+1)*(4*rb+2)))*(pi*abs(w-w_)*(N))^(2*rb+1)/gamma(2*rb+2) > deltab)
    rb = rb+1;
end

rb=1;
vecGamma = gamma(1:(2*rb));
vecCoeff = (2/((N)*pi))*(-1).^(0:(rb-1)).*(pi*(w-w_)*(N)).^(1:2:(2*rb-1))./gamma(2:2:(2*rb));
Cb = zeros(2*rb-1,2*rb-1);
for k = 0:(rb-1)
    for l = 0:(2*k)
        Cb(l+1,2*k-l+1) = (-1)^(2*k-l)*vecCoeff(k+1).*vecGamma(2*k+1)./(vecGamma(l+1)*vecGamma(2*k+1-l));
    end
end

% Va and Vb share many columns
V = ones(N,max(2*ra,2*rb-1));
vecn = (0:(1/(N)):(1-1/(N)))';
for k = 1:max(2*ra-1,2*rb-2)
    V(:,k+1) = V(:,k).*vecn;
end

%% Form FST matrices
DAsinZ = bsxfun(@times,DAsin,Z);
DAcosZ = bsxfun(@times,DAcos,Z);
DAsinVa = bsxfun(@times,DAsin,V(:,1:(2*ra)));
DAcosVa = bsxfun(@times,DAcos,V(:,1:(2*ra)));
DBsinVb = bsxfun(@times,DBsin,V(:,1:(2*rb-1)));
DBcosVb = bsxfun(@times,DBcos,V(:,1:(2*rb-1)));

%% Compute DPSSs and store transition region DPSSs
cutoff = K-6;
[S,lambda,lowerindex,upperindex] = transitionDPSS(N,w,epsilonS,cutoff);
if(cutoff >= 1)
    K = cutoff;
    eig_weights = [ones(K-lowerindex+1,1);zeros(upperindex-K,1)]-lambda;
elseif(cutoff > eps)
    eig_weights = (lambda > cutoff)-lambda;
end

L1U1 = [DAsinZ,DAcosZ,DAsinZ(end:-1:1,:),DAcosZ(end:-1:1,:),DAsinVa*Ca,-DAcosVa*Ca,DBcosVb*Cb,DBsinVb*Cb,bsxfun(@times,S,eig_weights')];
L2U2 = [DAcosZ(end:-1:1,:),DAsinZ(end:-1:1,:),DAcosZ,DAsinZ,DAcosVa,DAsinVa,DBcosVb,DBsinVb,S]';

alpha_fast = FST(s,N,R,L2U2);
F_NW = zeros(N,2*R+1);
exp_freq = linspace(-R,R,2*R+1);

for ii=0:N-1
    for j=-R:R
        freq = j/(2*N);
        F_NW(ii+1,j+R+1) = exp(-1i*2*pi*freq*ii);
    end
end

K1 = size(L1U1,2)+2*R+1;
delta = 1e-12;
slepian_fast = [F_NW/sqrt(N) L1U1];
alpha_fast2 = inv(slepian_fast'*slepian_fast + delta*eye(K1))*(slepian_fast'*s);
s_estimate_fast2 = slepian_fast*alpha_fast2;

%% supporting functions
function [X] = sig_gen(a,f,n_sinusoids,t)
X = zeros(size(t));
for ii = 1:n_sinusoids % For each antenna
    sig = (a(ii)*exp(1i*2*pi*(f(ii))*(t)));
    X = X + sig;
end
end

function y = FST(x,N,R,L2U2)
    if(N ~= size(x,1))
        disp('ERROR: Input has incorrect length')
    else
        N = size(x,1);
        Fx = fft(x,N,1);
        y = [Fx([1:(R+1),(N-R+1):N],:);(L2U2)*x];
    end
end