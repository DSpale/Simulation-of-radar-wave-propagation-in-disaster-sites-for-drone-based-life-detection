clear
close all

% Adjust image path
image_path = 'D:\Uni\10. Semester\Masterarbeit\images\FCCW\';
if ~exist(image_path, 'dir')
    mkdir(image_path)
end
[~, simulation_directory, ~] = fileparts(pwd);

save_path = 'baseband_signals';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end

set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Times New Roman');

c = 299792458; % m/s

data_info = h5info('FCCW_bioradar1.out').Attributes;

% Values associated to paremeters in the h5info
dx_dy_dz = data_info(5).Value;
dt = data_info(6).Value;

% Characteristics of the FCCW frequency comb
FCCW.BW = 59.52e6; % MHz
FCCW.frequencies_amount = 32;
FCCW.frequency_step = FCCW.BW / (FCCW.frequencies_amount - 1);
FCCW.period_length = 1/FCCW.frequency_step;
FCCW.carrier_frequency = 1.3e9;

FCCW.starting_frequency = -29.28e6; % in the baseband!
FCCW.frequencies = linspace(FCCW.starting_frequency, FCCW.starting_frequency + FCCW.BW, FCCW.frequencies_amount) + FCCW.carrier_frequency;

% Range resolution
delta_R = c/(2*FCCW.BW);

% buffer for computational inaccuries towards the end of an array; adjust
% according to needs
period_sample_length = ceil(1/(FCCW.frequency_step * dt)) + 100;

simulations = length(dir(fullfile(pwd, '*.out')));

f_s_frame = 61.44e6; % MHz, sampling frequency for samples in frame
T_s_frame = 1/f_s_frame; % s, distance between two samples in frame
datapoints_per_frame = 4096; % samples per frame

% Time vector for the sampled times
sample_times = (T_s_frame:T_s_frame:datapoints_per_frame*T_s_frame)';

% Load data
data = cell(1, 6);
for i_sim = 1:simulations
    summed_Ex = zeros(datapoints_per_frame, 1);
    summed_Ey = zeros(datapoints_per_frame, 1);
    summed_Ez = zeros(datapoints_per_frame, 1);
    summed_Hx = zeros(datapoints_per_frame, 1);
    summed_Hy = zeros(datapoints_per_frame, 1);
    summed_Hz = zeros(datapoints_per_frame, 1);

    % Read data and truncate to >= last period. 100 last values are never
    % used in the processing algorithm (to mitigate edge effects from the
    % Hilbert transform)
    data_name = strcat('FCCW_bioradar', string(i_sim), '.out');
    Ex = h5read(data_name, '/rxs/rx1/Ex');
    Ex = hilbert(Ex);
    Ex = Ex(end-period_sample_length:end);

    Ey = h5read(data_name, '/rxs/rx1/Ey');
    Ey = hilbert(Ey);
    Ey = Ey(end-period_sample_length:end);

    Ez = h5read(data_name, '/rxs/rx1/Ez');
    Ez = hilbert(Ez);
    Ez = Ez(end-period_sample_length:end);

    Hx = h5read(data_name, '/rxs/rx1/Hx');
    Hx = hilbert(Hx);
    Hx = Hx(end-period_sample_length:end);

    Hy = h5read(data_name, '/rxs/rx1/Hy');
    Hy = hilbert(Hy);
    Hy = Hy(end-period_sample_length:end);

    Hz = h5read(data_name, '/rxs/rx1/Hz');
    Hz = hilbert(Hz);
    Hz = Hz(end-period_sample_length:end);

    for i_sample = 1:datapoints_per_frame
        % Get real data by interpolating linearly between two nearest
        % samples a and b
        curr_sample_time = sample_times(i_sample);
        curr_sample = mod(curr_sample_time, FCCW.period_length) / dt;
        curr_iteration = floor(curr_sample_time * FCCW.frequency_step);

        sample_a = floor(curr_sample)+1;
        sample_b = sample_a + 1;
        weight_sample_b = mod(curr_sample, 1);
        weight_sample_a = 1 - weight_sample_b;

        % Construct 4096 samples for the signal processing algorithm with
        % phase shift correction for beat period
        summed_Ex(i_sample) = (weight_sample_a * Ex(sample_a) + weight_sample_b * Ex(sample_b))*exp(-1j*pi/3*curr_iteration);
        summed_Ey(i_sample) = (weight_sample_a * Ey(sample_a) + weight_sample_b * Ey(sample_b))*exp(-1j*pi/3*curr_iteration);
        summed_Ez(i_sample) = (weight_sample_a * Ez(sample_a) + weight_sample_b * Ez(sample_b))*exp(-1j*pi/3*curr_iteration);
        summed_Hx(i_sample) = (weight_sample_a * Hx(sample_a) + weight_sample_b * Hx(sample_b))*exp(-1j*pi/3*curr_iteration);
        summed_Hy(i_sample) = (weight_sample_a * Hy(sample_a) + weight_sample_b * Hy(sample_b))*exp(-1j*pi/3*curr_iteration);
        summed_Hz(i_sample) = (weight_sample_a * Hz(sample_a) + weight_sample_b * Hz(sample_b))*exp(-1j*pi/3*curr_iteration);

    end
    % Correction of factor 2 & application of frequency shift property
    % (shifting to baseband)
    summed_Ex = summed_Ex / 2 .* exp(-2j*pi*FCCW.carrier_frequency*sample_times);
    summed_Ey = summed_Ey / 2 .* exp(-2j*pi*FCCW.carrier_frequency*sample_times);
    summed_Ez = summed_Ez / 2 .* exp(-2j*pi*FCCW.carrier_frequency*sample_times);
    summed_Hx = summed_Hx / 2 .* exp(-2j*pi*FCCW.carrier_frequency*sample_times);
    summed_Hy = summed_Hy / 2 .* exp(-2j*pi*FCCW.carrier_frequency*sample_times);
    summed_Hz = summed_Hz / 2 .* exp(-2j*pi*FCCW.carrier_frequency*sample_times);

    % Save intermediate results
    save_data = [summed_Ex, summed_Ey, summed_Ez, summed_Hx, summed_Hy, summed_Hz];
    save(strcat('.\', save_path, '\summed', string(i_sim), '.mat'), 'save_data');
end

%% BIORADAR SIGNAL PROCESSING
frequencies_baseband = FCCW.frequencies - FCCW.carrier_frequency;

% SORTIE FCCW processing
data = cell(simulations, 4);
for i_sim = 1:simulations
    % Time domain data
    loaded_data = load(strcat('.\', save_path, '\summed', string(i_sim), '.mat'));
    loaded_data = loaded_data.save_data;
    loaded_data = loaded_data(:, 3); % for Ez component
    data{i_sim, 1} = loaded_data;
    % Frequency domain data
    data{i_sim, 2} = fftshift(fft(loaded_data)) / datapoints_per_frame;
    % iFFT to spatial domain
    data{i_sim, 3} = zeros(32, 1); % frequency components
    data{i_sim, 4} = zeros(32, 1); % distances
    for i_frequency_component = 1:32 % samples where the different frequencies are
        % If beat period phase correction is not applied: Use 118 instead
        % of 97
        frequency_component = data{i_sim, 2}(97+(i_frequency_component-1)*128);
        data{i_sim, 3}(i_frequency_component) = frequency_component;
        data{i_sim, 4} = data{i_sim, 4} + abs(frequency_component)*exp(1j*2*pi/32*(i_frequency_component-1).*(0:31)');
    end
end

%%
% Construct one cycle, only first 17 samples contain all data
spatial_domain_data = zeros(17, 2*simulations-2);

for i_sim = 1:simulations
    spatial_domain_data(:, i_sim) = real(data{i_sim, 4}(1:17));
end
for i_sim = (simulations+1):(2*(simulations)-2)
    spatial_domain_data(:, i_sim) = real(data{2*simulations - i_sim, 4}(1:17));
end

N_frames = 20e3;
repetitions = 3;
frames_data = repmat(spatial_domain_data, 1, repetitions);
save(strcat(pwd, '\', save_path, '\frames_data', '.mat'), 'frames_data');

%% DATA PLOTTING OF THE INDIVIDUAL SIGNAL PROCESSING STEPS
close all
f_fft = f_s_frame/datapoints_per_frame*(-datapoints_per_frame/2:datapoints_per_frame/2-1);
figure()
hold on
grid on
% Quarter frame
plot((0:1027)/f_s_frame*1e6, real(data{1, 1}(1:1028))/1e3)
plot((0:1027)/f_s_frame*1e6, imag(data{1, 1}(1:1028))/1e3)
legend('Real part', 'Imaginary part')
axis padded; xlim('tight')
xlabel('Time in µs')
ylabel('$E_\mathrm{z}$ in kV/m', 'interpreter', 'latex')
set(gca, 'FontSize', 16)
exportgraphics(gcf,strcat(image_path, simulation_directory, '_time_domain.pdf'), 'ContentType','vector')

% Frequency domain data
figure()
hold on
grid on
plot(f_fft/1e6, real(data{1, 2}))
plot(f_fft/1e6, imag(data{1, 2}))
legend('Real part', 'Imaginary part')
axis padded; xlim('tight')
xlabel('Frequency in MHz', 'interpreter', 'latex')
ylabel('$E_\mathrm{z}$ in V/m', 'interpreter', 'latex')
set(gca, 'FontSize', 16)
exportgraphics(gcf,strcat(image_path, simulation_directory, '_frequency_domain1.pdf'), 'ContentType','vector')

% To only show the 32 relevant components
figure()
hold on
grid on
f_fft_2 = linspace(f_fft(97), f_fft(4065), length(data{1, 3}));
plot(f_fft_2/1e6, real(data{1, 3}), '-x')
plot(f_fft_2/1e6, imag(data{1, 3}), '-x')
plot(f_fft_2/1e6, abs(data{1, 3}), '-x')
legend('Real part', 'Imaginary part', 'Magnitude')
axis padded; xlim('tight')
y_lim = ylim;
ylim([y_lim(1) y_lim(2) + 400])
clear y_lim
xlabel('Frequency in MHz', 'interpreter', 'latex')
ylabel('$E_\mathrm{z}$ in V/m', 'interpreter', 'latex')
set(gca, 'FontSize', 16)
exportgraphics(gcf,strcat(image_path, simulation_directory, '_frequency_domain2.pdf'), 'ContentType','vector')

% Spatial domain
figure()
hold on
grid on
plot((1:31)*delta_R, real(data{1, 4}(2:end)), '-x')
plot((1:31)*delta_R, imag(data{1, 4}(2:end)), '-x')
xline(16*delta_R,'-',{'Axis of','symmetry'}, 'LabelHorizontalAlignment', 'left', 'FontSize', 16, 'FontName', 'CMU Serif');
legend('Real part', 'Imaginary part')
axis padded
xlabel('Range in m', 'interpreter', 'latex')
ylabel('IFFT of radar frequency domain', 'interpreter', 'latex')
set(gca, 'FontSize', 16)
exportgraphics(gcf,strcat(image_path, simulation_directory, '_spatial_domain.pdf'), 'ContentType','vector')

% Range bin plots. Adjust range bins (here: i = 5:7) according to range
% bins (i = 1 - "range bin 0" is not used in the processing algorithm)
figure()
tl = tiledlayout(3, 1);
for i = 5:7
    nexttile
    hold on
    grid on
    box on
    plot(frames_data(i, :) - mean(frames_data(i, :)))
    ylabel(strcat("Range ", string(i-1)), 'interpreter', 'latex')
    axis padded; xlim('tight')
    set(gca, 'FontSize', 16)
end
xlabel('Frame number', 'interpreter', 'latex')
exportgraphics(gcf,strcat(image_path, simulation_directory, '_range_frames.pdf'), 'ContentType','vector')

% Obtain peak-to-peak breathing activities
range_activity = zeros(16, 1);
for i = 1:16
    range_activity(i) = max(frames_data(i+1, :)) - min(frames_data(i+1, :));
end

figure()
hold on; grid on
scatter((1:16)*delta_R, range_activity, 'filled')
xlabel("Range in m", 'interpreter', 'latex')
ylabel("Peak-to-peak breathing activity", 'interpreter', 'latex')
axis padded
set(gca, 'FontSize', 16)
exportgraphics(gcf,strcat(image_path, simulation_directory, '_range_estimate.pdf'), 'ContentType','vector')

% Calculate prominence ratio PR and save it to .txt-file
PR_squared = range_activity.^2;
PR = (max(PR_squared)*(length(PR_squared)-1))/(sum(PR_squared) - max(PR_squared));

% Write prominence ratio to file
txt_file = fopen(strcat('PR.txt'), 'wt');
fprintf(txt_file, 'Prominence ratio: %.4f', PR);
fclose(txt_file);