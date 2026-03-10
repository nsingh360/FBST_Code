
function broadband_dvm_test()
% BROADBAND_DVM_TEST - Tests DVM over multiple frequencies to simulate broadband behavior

% Parameters
N = 8;                          % Size of input
tau = 0.1;                      % Delay
Fs = 1;                         % Sampling frequency (normalized)
num_freqs = 16;                 % Number of frequency bins
frequencies = linspace(0, Fs/2, num_freqs);  % Frequency range

% Input signal z (same across frequencies)
z = randn(N, 1);

% Store outputs
Y = zeros(N, num_freqs);

% Run DVM at each frequency
for k = 1:num_freqs
    omega = 2 * pi * frequencies(k);
    alpha = exp(-1j * omega * tau);
    Y(:, k) = dvm(N, alpha, z);
end

% Display results
disp('DVM outputs for each frequency bin:');
disp(abs(Y)); % show magnitudes
imagesc(frequencies, 1:N, abs(Y));
xlabel('Frequency (Hz)');
ylabel('Element Index');
title('Magnitude of DVM Output vs Frequency');
colorbar;

end
