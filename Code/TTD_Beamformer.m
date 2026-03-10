function [Y] = TTD_Beamformer(X,delays, fs, Psi, varargin)
% Function generates the filters associated with a given non-integer delay
%
%%%% Required Inputs %%%%%
%   fs - sampling frequnecy
%   Psi - Interpolation kernel. Need to accept MxN matrix inputs
%%%% Optional Inputs %%%%%
%   R - length of filter is 2*R

%%%%%%%%%%%%% Argument Parser %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;

addRequired(p,'X');
addRequired(p,'delays');
addRequired(p,'fs');
addRequired(p,'Psi');
addOptional(p,'R', 4, @isnumeric);

parse(p, X, delays, fs, Psi, varargin{:});
X = p.Results.X;
delays = p.Results.delays;
fs = p.Results.fs;
Psi = p.Results.Psi;
R = p.Results.R;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N,M] = size(X);
N_filter = N+2*R+1;

[out] = Kernel_Filter_Gen(fs, Psi,'R',R);
[h, h_idx] = out.Filter_Gen(delays);
h = squeeze(h);
h_idx = squeeze(h_idx);

wrapN = @(x, n) (1 + mod(x-1, n));

H = zeros(N_filter,M);

for i=1:M
    H(wrapN(h_idx(i,:),N_filter),i) = h(i,:);
end

Hf = fft(H,[],1);
Xf = fft(X,N_filter,1);
Y = ifft(Hf.*Xf,[],1);
Y = Y(1:N,:);

end

