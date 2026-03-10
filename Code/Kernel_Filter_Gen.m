function [out] = Kernel_Filter_Gen(fs, Psi, varargin)
% Function generates the filters associated with a given non-integer delay
%
%%%% Required Inputs %%%%%
%   fs - sampling frequnecy
%   Psi - Interpolation kernel. Need to accept MxN matrix inputs
%%%% Optional Inputs %%%%%
%   R - length of filter is 2*R

%%%%%%%%%%%%% Argument Parser %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;

addRequired(p,'fs');
addRequired(p,'Psi');
addOptional(p,'R', 4, @isnumeric);

parse(p, fs, Psi, varargin{:});
fs = p.Results.fs;
Psi = p.Results.Psi;
R = p.Results.R;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% window of indices
idx_window = [-R+1:R]';
% System of equation matrix
K = zeros(2*R,2*R);
% generate system of equations
for ii = 1:length(idx_window)
    v = Psi(idx_window(ii))*ones(1,2*R-abs(idx_window(ii)));
    K = K + diag(v,idx_window(ii));
end

% Compute inverse
K_inv = K^-1;

% save results
out.K_inv = K_inv;
out.Psi = Psi;
out.fs = fs;
out.R = R;
out.idx_window = idx_window;
out.Filter_Gen = @Filter_Gen;

% function that generates filters and the tap locations
    function [h, h_idx] = Filter_Gen(tau)
        [M,N] = size(tau);
        % by convention f(t-tau) so
        tau = -tau;
        % vectorize
        tau_vec = reshape(tau,M*N,1);
        % integer delay
        tau_int = floor(tau_vec*out.fs);
        % fractional delay
        tau_frac = tau_vec - floor(tau_vec*out.fs)/out.fs;
        % calculate b matrix
        b = out.Psi(out.idx_window' - tau_frac*out.fs);
        % generate filter
        psi_tilde_vec = (conj(b)*out.K_inv);
        % tap locations
        h_idx_vec = out.idx_window'-tau_int;
        % reshape in an M x N x 2R matrix, adjust so that reference is 0
        % not 1
        h_idx = reshape(h_idx_vec,M,N,2*out.R)-1;
        % flip vectors since they are presumabley applied via convolution
        h = reshape(fliplr(psi_tilde_vec),M,N,2*out.R);
    end

end

