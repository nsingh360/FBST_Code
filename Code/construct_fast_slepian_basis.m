function [A_basis2, V_nyq2, K2] = construct_fast_slepian_basis(W, t, T, L, M, N, T_ap)

MN = M*N;
N_u = ceil(T_ap/T);
K = ceil(2*W*T*(N_u+1));

t_u = [0:(N_u)]'; % uniform sampling points over the interval
t_nu = [N_u*(t-min(t(:)))/T_ap]; % non-uniform sample points, scaled to match interval

ks = reshape(linspace(-ceil(K/2),floor(K/2)-1,K),[K,1]);
F_K = ks*(1/2)*(1/(L*T));

ls = reshape(linspace(-L,L-1,2*L),[2*L,1]);
F_L = ls*(1/2)*(1/(L*T));

A_basis = zeros(2*L,K);

for ll=1:2*L
    for kk=1:K
        A_basis(ll,kk) = sum(exp(1i*2*pi*(F_K(kk) - F_L(ll))*t));
    end
end

fs = 1/T;
t_nyq = [0:N-1]'/fs; % nyquist samples for testing
V_nyq = exp(1i*2*pi*t_nyq*F_K.');

% hilbert matrix approximation
H = zeros(N_u+1,N_u+1);
for i=1:N_u+1
    for j=1:N_u+1
        H(i,j) = 1/(i+j-1);
    end
end
Z = approx_H(N_u+1);

% A1 matrix approximatin using Taylor series
A1 = zeros(N_u+1,N_u+1);
for i=1:N_u+1
    for j=1:N_u+1
        if i==j
            A1(i,j) = 0;
        else
            A1(i,j) = 1/(pi*(i-j)) - 1/((N_u+1)*sin(pi*(i-j)/(N_u+1))) - 1/(pi*(i-j + N_u+1)) - 1/(pi*(i-j - (N_u+1)));
        end
    end
end
[VA,CA] = approx_A1(N_u+1);

% B0 matrix approximation using Taylor series
w = W/fs;
if mod(ceil(2*(N_u+1)*w),2)~=0
    w_ = ceil(2*(N_u+1)*w)/(2*(N_u+1));
else
    w_ = (ceil(2*(N_u+1)*w)+1)/(2*(N_u+1));
end
B0 = zeros(N_u+1,N_u+1);
for i=1:N_u+1
    for j=1:N_u+1
        if i==j 
            B0(i,j) = 2*(w-w_);
        else
            B0(i,j) = 2*sin(pi*(w-w_)*(i-j))/(pi*(i-j));
        end
    end
end

[VB,CB] = approx_B0(B0,N_u+1,w,w_);
DA = diag(exp(1i*2*pi*w_*t_u));
DB = diag(exp(1i*pi*(w+w_)*t_u));

% create exchange matrix
J = zeros(N_u+1,N_u+1);
for i=1:N_u+1
    J(i,N_u-i+2) = 1;
end
L1 = [DA*J*Z DA*Z DA'*J*Z DA'*Z DA*VA*CA' DA'*VA*CA' DB*VB*CB' DB'*VB*CB'];
L2 = [(1/(2*pi*1i))*DA*Z (-1/(2*pi*1i))*DA*J*Z (-1/(2*pi*1i))*DA'*Z (1/(2*pi*1i))*DA'*J*Z (1/(2*1i))*DA*VA (-1/(2*1i))*DA'*VA (1/2)*DB*VB (1/2)*DB'*VB];
L_size = size(L1);
L_dim = L_size(2);

%% Error on individual terms (somewhat close to optimal)
epsilon = 1;
epsilonF = 2*epsilon/3; % Error in approximating B_{N,W}-F_{N,W}F_{N,W}^*
mu = 2*log(8*N_u-4)/pi^2+1/log(2)+1/log(6/pi);
deltah = (log(8*N_u-4)/pi)*epsilonF/mu;
deltaa = (1/log(2))*epsilonF/mu;
deltab = (1/log(6/pi))*epsilonF/mu;
epsilonS = epsilon/3; % Error in approximating B_{N,W}-S_KS_K^* 

%% Diagonal matrices (store as column vectors)
R = round(w*(N_u+1)-1/2-1e-3);
w_ = (R+1/2)/(N_u+1);
vecn = ((-(N_u+1-1)/2):((N_u+1-1)/2))';
DAcos = cos((2*pi*w_)*vecn);
DAsin = sin((2*pi*w_)*vecn);
DBcos = cos((pi*(w+w_))*vecn);
DBsin = sin((pi*(w+w_))*vecn);

%% Approximate Scaled Hilbert Matrix (1/pi)*H ~~ Z*Z'
rh = ceil(log(8*(N_u+1)-4)*log(4*pi/deltah)/pi^2);

k1 = 1/(2*(N_u+1)-1);
k = sqrt(1-k1^2);
Kk = ellipticK(k^2);
[~,~,p] = ellipj((1:2:(2*rh-1))*Kk/(2*rh),k^2);
p = (N_u+1-1/2)*p;

Z = zeros(N_u+1,rh);
vecn = ((1/2):(N_u+1-1/2))';
Z(:,1) = sqrt(2*p(1)/pi)*ones(N_u+1,1)./(vecn+p(1));
for j = 2:rh,
    Z(:,j) = sqrt(p(j)/p(j-1))*(1-(p(j)+p(j-1))./(vecn+p(j))).*(Z(:,j-1));
end

%% Approximate A1 ~~ Va*Ca*Va' and B0 ~~ Vb*Cb*Vb'
ra = ceil(log(2/(3*pi*deltaa))/(2*log(2)));
while((2/(3*pi))*sqrt(2/((4*ra-1)*4*ra))/2^(2*ra-2) <= deltaa)
    ra = ra-1;
end

vecZeta = zeta(2:2:(2*ra));
vecCoeff = (2/((N_u+1)*pi))*(1-(1-2.^(-(1:2:(2*ra-1)))).*vecZeta);
vecGamma = gamma(1:(2*ra));
Ca = zeros(2*ra,2*ra);
for k = 1:ra,
    for l = 0:(2*k-1),
        Ca(l+1,2*k-l) = (-1)^(2*k-1-l)*vecCoeff(k).*vecGamma(2*k)./(vecGamma(l+1)*vecGamma(2*k-l));
    end
end

rb = 0;
while((2/pi)*sqrt(2/((4*rb+1)*(4*rb+2)))*(pi*abs(w-w_)*(MN+1))^(2*rb+1)/gamma(2*rb+2) > deltab)
    rb = rb+1;
end

vecGamma = gamma(1:(2*rb));
vecCoeff = (2/((N_u+1)*pi))*(-1).^(0:(rb-1)).*(pi*(w-w_)*(N_u+1)).^(1:2:(2*rb-1))./gamma(2:2:(2*rb));
Cb = zeros(2*rb-1,2*rb-1);
for k = 0:(rb-1),
    for l = 0:(2*k),
        Cb(l+1,2*k-l+1) = (-1)^(2*k-l)*vecCoeff(k+1).*vecGamma(2*k+1)./(vecGamma(l+1)*vecGamma(2*k+1-l));
    end
end

% Va and Vb share many columns
V = ones(N_u+1,max(2*ra,2*rb-1));
vecn = (0:(1/(N_u+1)):(1-1/(N_u+1)))';
for k = 1:max(2*ra-1,2*rb-2),
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
cutoff = K;
[S,lambda,lowerindex,upperindex] = transitionDPSS(N_u+1,w,epsilonS,cutoff);
if(cutoff >= 1)
    K = cutoff;
    eig_weights = [ones(K-lowerindex+1,1);zeros(upperindex-K,1)]-lambda;
elseif(cutoff > eps)
    eig_weights = (lambda > cutoff)-lambda;
end

L1 = [DAsinZ,DAcosZ,DAsinZ(end:-1:1,:),DAcosZ(end:-1:1,:),DAsinVa*Ca,-DAcosVa*Ca,DBcosVb*Cb,DBsinVb*Cb,bsxfun(@times,S,eig_weights')];
% L1 = [L1, bsxfun(@times,S,eig_weights')];
L_size = size(L1);
L_dim = L_size(2);


L1_nu = interp1(t_u,L1,t_nu,'pchip',0); % non-uniform
exp_basis = exp(-1i*2*pi*(F_L*t.'));
L1_basis = exp_basis*L1_nu;

L1_nyq = interp1(t_u,L1,N_u*(t_nyq-min(t))/T_ap,'pchip',0); % interpolate to nyquist

A_basis2 = [A_basis L1_basis];
V_nyq2 = [V_nyq L1_nyq];
K2 = K + L_dim;

end

function [VB,CB] = approx_B0(B,N,W,W_)
rB = 1;
VB = zeros(N,2*rB-1);
CB = zeros(2*rB-1,2*rB-1);

for n=1:N
    for r=1:2*rB-1
        VB(n,r) = (n-1)^(r-1);
    end
end
VB(1,1) = 0;

for k_=1:rB
    k = k_ - 1;
    for r1 = 1:2*k+1
        l = r1 - 1;
        CB(r1,2*k-l+1) = CB(r1,2*k-l+1) + 2*(1/N)*(1/pi)*((-1)^(3*k-l))*(((pi*(W-W_)*N)^(2*k+1))/factorial(2*k+1))*nchoosek(2*k,l)*((1/N)^(2*k));
    end
end
end

function [VA,CA] = approx_A1(N)
rA = 1;
VA = zeros(N,2*rA);
CA = zeros(2*rA,2*rA);

for n=1:N
    for r=1:2*rA
        VA(n,r) = (n-1)^(r-1);
    end
end
VA(1,1) = 0;

for k=1:rA
    for r1=1:2*k
        l = r1-1;
        CA(r1, 2*k-l) = CA(r1, 2*k-l) + 2*(1/N)*(1/pi)*((1 - (1 - 2^(-(2*k-1)))*zeta(2*k))*nchoosek(2*k-1,l)*((-1)^(2*k-1-l))*((1/N)^(2*k-1)));
    end
end

end

function [Z] = approx_H(N)

A = zeros(N,N);
B = ones(N,1);
for i=1:N
    A(i,i) = i-1 + 1/2;
end

a = 1/2;
b = N - 1/2;
cond_A = 2*N-1;


r = 3;

% select equispaced (in log scale) points 
P = zeros(r,1);
for k=1:r
    P(k,1) = a^((2*k-1)/(2*r))*b^((2*r-2*k+1)/(2*r));
end

Z = zeros(N,r);
for k=1:r
    if k==1
        Z(:,k) = sqrt(2*P(k))*inv(A + P(k)*eye(N))*B;
    else
        Z(:,k) = sqrt(P(k)/P(k-1))*(eye(N) - (P(k) + P(k-1))*inv(A + P(k)*eye(N)))*Z(:,k-1);
    end
end

end