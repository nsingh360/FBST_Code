function [alpha,error,error_bound,M,rho] = fe_bias_error(sig_gen, fe_continuous, a, f, n_sinusoids, T_extended, fs, W, L, t_samples,transform)

min_t = min(t_samples);
max_t = max(t_samples);
T_aperture = max(t_samples) - min(t_samples);

%% calculate continuous time estimate for alpha
Gram_matrix = zeros(2*L,2*L);
f_vec = zeros(2*L,1);
for l1=1:2*L
    for l2=1:2*L
        f_l1 = (fs/(2*L))*(l1-1-L);
        f_l2 = (fs/(2*L))*(l2-1-L);
        Gram_matrix(l1,l2) = integral(@(t)exp(1i*2*pi*(f_l2-f_l1)*t),min(t_samples),max(t_samples),'RelTol', 1e-12, 'AbsTol', 1e-14);
    end
    f_vec(l1) = integral(@(t)sig_gen(a,f,n_sinusoids,t).*exp(-1i*2*pi*f_l1*t),min(t_samples),max(t_samples),'RelTol', 1e-12, 'AbsTol', 1e-14);
end

% Gram_matrix = vpa(Gram_matrix,100);
% f_vec = vpa(f_vec, 100);
alpha = (Gram_matrix)\f_vec; % avoid using inv() for numerical stability.

% %% measure error
error = sqrt(integral(@(t)(abs(sig_gen(a,f,n_sinusoids,t) - fe_continuous(alpha,fs,L,t))).^2,min(t_samples),max(t_samples)));

%% bound error
rho = cot(pi/(4*T_extended))^2;
theta_pts = linspace(-pi, pi, 1000);
berstein_boundary = 0.5*(rho*exp(1i*theta_pts) + (1/rho)*exp(-1i*theta_pts));

c_T = cos(pi/T_extended);
mapped_boundary = T_extended/pi*acos(c_T + 0.5*(1-c_T)*(berstein_boundary+1));
if transform==1
    mapped_boundary = (T_aperture/2)*(mapped_boundary + (max_t+min_t)/T_aperture);
end
s_complex = sig_gen(a,f,n_sinusoids,mapped_boundary);
[~,idx] = max(abs(s_complex));
M = max(abs(s_complex));

c2 = 2*pi*abs(imag(mapped_boundary(idx)));
x = linspace(-1,1,1000);
if transform==1
    x = (T_aperture/2)*(x + (max_t+min_t)/T_aperture);
end
c1 = max(abs(sig_gen(a,f,n_sinusoids,x)));
error_bound = c1*exp(c2*W)*rho^(-L)/(rho-1);

end