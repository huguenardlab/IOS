% Parameters
function template_out=WhiskerTemplate(total_pulses, pulse_frequency, pulse_duration, delay_duration, rise_time,sample_rate,total_duration,varargin)
pulse_increment=0;
if (size(varargin,1)>0)
    pulse_increment=varargin{1};
end
pulse_period=1/pulse_frequency;


% Time vector for the entire pulse sequence, but can be shorter than the
% total duration
t_total = 0:1/sample_rate:(total_pulses * (pulse_period) + delay_duration) - 1/sample_rate;

% Generate square pulse signal with linear ramp up and down
template_out = zeros(total_duration*sample_rate,1)';
for i = 1:total_pulses
    pulse_amplitude=1+(i*pulse_increment);
    t_pulse_start = delay_duration + (i - 1) *  pulse_period;
    t_pulse_end = t_pulse_start + pulse_duration+rise_time*2; % plateau will be 50 ms
    
    % Linear ramp up
    ramp_up = (t_total - t_pulse_start) / rise_time .* (t_total >= t_pulse_start & t_total <= (t_pulse_start + rise_time));
    
    % Plateau
    plateau = (t_total > (t_pulse_start + rise_time)) & (t_total <= (t_pulse_end - rise_time));
    
    % Linear ramp down
    ramp_down = (1 - (t_total - (t_pulse_end - rise_time)) / rise_time) .* (t_total > (t_pulse_end - rise_time) & t_total <= t_pulse_end);
    
    template_out(1:size(t_total,2)) = template_out(1:size(t_total,2)) + pulse_amplitude*(ramp_up + plateau + ramp_down);
end

