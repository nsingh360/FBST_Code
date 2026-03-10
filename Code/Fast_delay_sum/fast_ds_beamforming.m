function [beamformed_signals, Angles] = fast_ds_beamforming(x, d, fs, c, angles, num_stages)
% Implements Fast Delay-and-Sum (Fast DS) beamforming.
% 
% Parameters:
% x : matrix
%     Received signal at each sensor (size: M x N, where M is number of sensors, N is samples).
% d : float
%     Sensor spacing in meters.
% fs : float
%     Sampling frequency in Hz.
% c : float
%     Propagation speed (e.g., speed of sound or light).
% angles : vector
%     Array of steering angles (in degrees).
% num_stages : int
%     Number of recursive beamforming stages.
% 
% Returns:
% energy : vector
%     Beamformed energy at each steering angle.

[M, N] = size(x); % Number of sensors and time samples
energy = zeros(1, length(angles));
Angles = zeros(2^num_stages,1);
beamformed_signals = zeros(N,2^num_stages);

for s = 1:num_stages
    num_sectors = 2^s;
    for sector = 1:num_sectors
        theta = asind(1 - (2 * (sector - 1) + 1) / (2^s)); % Compute steering angle
        Angles(sector,1) = theta;
        delay_samples = (d * sind(theta) / c) * fs; % Delay per sensor in samples
        
        % Apply delay-and-sum beamforming
        delayed_signals = zeros(size(x));
        for m = 1:M
            fractional_delay = (m - 1) * delay_samples;
            int_delay = floor(fractional_delay);
            frac_part = fractional_delay - int_delay;
            
            if int_delay < N - 1
                delayed_signals(m, 1:end-int_delay) = (1 - frac_part) * x(m, int_delay+1:end) + frac_part * x(m, int_delay+2:end);
            end
        end
        
        % Sum across sensors
        beamformed_signal = sum(delayed_signals, 1);
        beamformed_signals(:,sector) = beamformed_signal;
        energy(sector) = sum(abs(beamformed_signal).^2); % Compute energy
    end
end

end