c = 299792458; % m/s
% Adjust spatial discretization according to use-case
space_res_x = 5e-4; % Space discretization
space_res_y = 2.5e-3;
space_res_z = min(space_res_x, space_res_y);

% According to CFL condition
t_step = 1/(c*sqrt(1/space_res_x^2 + 1/space_res_y^2));

% tsim needs to be bigger than the simulated time window. For minimal file
% size, use the length of the simulated time window.
tsim = 1000e-9;
iterations = ceil(tsim/t_step)+1;

f_s = 1/t_step; % Sampling frequency in Hz
t = 0:1/f_s:tsim;

% FCCW generation with FCCW parameters given in 
M = 32; % amount of frequency components
BW = 59.52e6; % bandwidth
delta_f = BW/(M-1); % frequency step
FCCW_frequencies = zeros(1, M);
FCCW_baseband_signal = zeros(1, length(t));
for m = 1:M
    m_frequency = -29.28e6 + (m-1)*delta_f;
    FCCW_frequencies(m) = m_frequency;
    FCCW_baseband_signal = FCCW_baseband_signal + cos(2*pi*m_frequency*t)+ 1i*sin(2*pi*m_frequency*t);
end
clear m m_frequency

% Shift introduced by carrier frequency
frequency_shift = 1.3e9;

% I/Q shift
LO = exp(1i*2*pi*frequency_shift*t);
I_Q_shifted_signal = real(FCCW_baseband_signal).*real(LO) + ...
    imag(FCCW_baseband_signal).*real(LO*1i);

% Introduce a ramp to excitation function
t_ramp = 4e-9;
ramp = (t_step:t_step:t_ramp) ./ t_ramp;
ramp(end+1:length(t)) = 1;

I_Q_shifted_signal = I_Q_shifted_signal .* ramp;

% Save the excitation signal in a text file to use in gprMax
writetable(table(t', I_Q_shifted_signal', 'VariableNames', {'time', 'FCCW_IQ_excitation'}), 'excitations_FCCW_IQ_ramp.txt', 'delimiter', ' ')