function [rf_pulse_res,t_rf] = rf_pulse_sinc(flip_angle,bandwidth,num_lobes,sampling_rate,phase)

% Derived parameters
pulse_duration = (num_lobes+1) / bandwidth; % Total pulse duration (s)
num_samples = round(pulse_duration * sampling_rate); % Total time steps
t_rf = linspace(-pulse_duration / 2, pulse_duration / 2, num_samples); % Time vector

% Generate sinc pulse
sinc_pulse = sinc(2 * bandwidth * t_rf);

% Scale the sinc pulse for the desired flip angle
gamma = 42.576e6; % Gyromagnetic ratio (Hz/T)
flip_angle_rad = flip_angle * pi / 180; % Convert flip angle to radians
area = trapz(t_rf, sinc_pulse); % Area under the sinc pulse
if area == 0
    area = 1e-12; % Prevent division by zero
end
scale_factor = flip_angle_rad / (gamma * abs(area)); % Scaling factor
sinc_pulse_scaled = scale_factor * sinc_pulse;

% Plot results
%     figure('Name','RF Pulse')
%     subplot(2, 1, 1);
%     plot(t_rf, sinc_pulse, 'LineWidth', 1.5);
%     xlabel('Time (s)');
%     ylabel('Amplitude (arb. units)');
%     title('Unscaled Sinc Pulse');
%     grid on;

%     subplot(2, 1, 2);
plot(t_rf, sinc_pulse_scaled, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('B1 Amplitude (T)');
title(sprintf('Scaled Sinc Pulse for %.1f° Flip Angle', flip_angle));
grid on;

rf_pulse = sinc_pulse_scaled;

rad_phase = phase * pi / 180;
phase_Rot_mat = [cos(rad_phase) -sin((rad_phase));
    sin((rad_phase)) cos((rad_phase))];

rf_pulse_with_phase = phase_Rot_mat * [rf_pulse ; zeros(1,length(rf_pulse))];
rf_pulse_res = rf_pulse_with_phase(1,:) + 1i * rf_pulse_with_phase(2,:);

end