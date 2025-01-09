clear
close all

% Adjust image path
image_path = 'D:\Uni\10. Semester\Masterarbeit\images\SP_FCCW\';
if ~exist(image_path, 'dir')
    mkdir(image_path)
end
[~, simulation_directory, ~] = fileparts(pwd);

% Save intermediate results for the superpositioned baseband signals
save_path = 'superpositioned_baseband_signals';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end

set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Times New Roman');

% Read the simulation data filenames
simulation_files = dir(fullfile(pwd, '*.out'));

% For each simulation, 32 different frequency components are simulated with
% one output file for each component
num_of_simulations = length(simulation_files) / 32;
sample_data = h5info(strcat(pwd, '\', 'SP_FCCW_radar1.out'));

iterations = double(sample_data.Attributes(3).Value);
dt = sample_data.Attributes(6).Value;
tsim = (iterations - 1) * dt;

c = 299792458; % m/s

% SORTIE FCCW bioradar parameters
BW = 59.52e6; % MHz
frequency_step = BW / 31;
frequencies = cell(32, 3);
carrier_frequency = 1.3e9;

% Range resolution
delta_R = c/(2*BW);

for m = 1:32
    % Frequencies of the frequency comb
    frequencies{m, 1} = carrier_frequency - 29.28e6 + (m-1)*frequency_step;
    % Period length for differing frequencies
    frequencies{m, 2} = 1/(frequencies{m, 1});
    % Phase shift in seconds for "Hilbert transform"
    frequencies{m, 3} = -0.25 * frequencies{m, 2};
end
max_iterations = ceil(frequencies{1, 2}/dt);

f_s_frame = 61.44e6; % MHz, sampling frequency for samples in frame
T_s_frame = 1/f_s_frame; % s, distance between two samples
datapoints = 4096; % samples per frame

% Time vector for the sampled times
sample_times = (T_s_frame:T_s_frame:datapoints*T_s_frame)';

% Load simulation data
data = cell(1, 6);
for i_sim = 1:num_of_simulations
    summed_Ex = zeros(datapoints, 1);
    summed_Ey = zeros(datapoints, 1);
    summed_Ez = zeros(datapoints, 1);
    summed_Hx = zeros(datapoints, 1);
    summed_Hy = zeros(datapoints, 1);
    summed_Hz = zeros(datapoints, 1);
    for i_frequency_component = 1:32
        % Get index of current file
        i_file = (i_sim-1)*32 + i_frequency_component;
        % Read data and truncate to >= last period
        data_name = strcat('SP_FCCW_radar', string(i_file), '.out');
        Ex = h5read(data_name, '/rxs/rx1/Ex');
        Ex = Ex(end-max_iterations:end);
        Ey = h5read(data_name, '/rxs/rx1/Ey');
        Ey = Ey(end-max_iterations:end);
        Ez = h5read(data_name, '/rxs/rx1/Ez');
        Ez = Ez(end-max_iterations:end);
        Hx = h5read(data_name, '/rxs/rx1/Hx');
        Hx = Hx(end-max_iterations:end);
        Hy = h5read(data_name, '/rxs/rx1/Hy');
        Hy = Hy(end-max_iterations:end);
        Hz = h5read(data_name, '/rxs/rx1/Hz');
        Hz = Hz(end-max_iterations:end);

        for i_sample = 1:datapoints
            % Get real data by interpolating linearly between two nearest
            % samples a and b
            curr_sample_time.r = sample_times(i_sample);
            curr_sample.r = mod(curr_sample_time.r, ...
                frequencies{i_frequency_component, 2}) / dt;

            sample_a.r = floor(curr_sample.r)+1;
            sample_b.r = sample_a.r + 1;
            weight_sample_b.r = mod(curr_sample.r, 1);
            weight_sample_a.r = 1 - weight_sample_b.r;

            % Perform same operations for imaginary data
            curr_sample_time.im = sample_times(i_sample) + ...
                frequencies{i_frequency_component, 3};
            curr_sample.im = mod(curr_sample_time.im, ...
                frequencies{i_frequency_component, 2}) / dt;

            sample_a.im = floor(curr_sample.im)+1;
            sample_b.im = sample_a.im + 1;
            weight_sample_b.im = mod(curr_sample.im, 1);
            weight_sample_a.im = 1 - weight_sample_b.im;

            % Accumulate points from different frequncy components
            summed_Ex(i_sample) = summed_Ex(i_sample) ...
                + (weight_sample_a.r * Ex(sample_a.r) + weight_sample_b.r * Ex(sample_b.r)) ...
                + 1j*(weight_sample_a.im * Ex(sample_a.im) + weight_sample_b.im * Ex(sample_b.im));

            summed_Ey(i_sample) = summed_Ey(i_sample) ...
                + (weight_sample_a.r * Ey(sample_a.r) + weight_sample_b.r * Ey(sample_b.r)) ...
                + 1j*(weight_sample_a.im * Ey(sample_a.im) + weight_sample_b.im * Ey(sample_b.im));

            summed_Ez(i_sample) = summed_Ez(i_sample) ...
                + (weight_sample_a.r * Ez(sample_a.r) + weight_sample_b.r * Ez(sample_b.r)) ...
                + 1j*(weight_sample_a.im * Ez(sample_a.im) + weight_sample_b.im * Ez(sample_b.im));

            summed_Hx(i_sample) = summed_Hx(i_sample) ...
                + (weight_sample_a.r * Hx(sample_a.r) + weight_sample_b.r * Hx(sample_b.r)) ...
                + 1j*(weight_sample_a.im * Hx(sample_a.im) + weight_sample_b.im * Hx(sample_b.im));

            summed_Hy(i_sample) = summed_Hy(i_sample) ...
                + (weight_sample_a.r * Hy(sample_a.r) + weight_sample_b.r * Hy(sample_b.r)) ...
                + 1j*(weight_sample_a.im * Hy(sample_a.im) + weight_sample_b.im * Hy(sample_b.im));

            summed_Hz(i_sample) = summed_Hz(i_sample) ...
                + (weight_sample_a.r * Hz(sample_a.r) + weight_sample_b.r * Hz(sample_b.r)) ...
                + 1j*(weight_sample_a.im * Hz(sample_a.im) + weight_sample_b.im * Hz(sample_b.im));
        end
    end
    % Correction of factor 2 & application of frequency shift property
    % (shift to baseband)
    summed_Ex = summed_Ex / 2 .* exp(-2j*pi*carrier_frequency*sample_times);
    summed_Ey = summed_Ey / 2 .* exp(-2j*pi*carrier_frequency*sample_times);
    summed_Ez = summed_Ez / 2 .* exp(-2j*pi*carrier_frequency*sample_times);
    summed_Hx = summed_Hx / 2 .* exp(-2j*pi*carrier_frequency*sample_times);
    summed_Hy = summed_Hy / 2 .* exp(-2j*pi*carrier_frequency*sample_times);
    summed_Hz = summed_Hz / 2 .* exp(-2j*pi*carrier_frequency*sample_times);

    % Save intermediate baseband results
    save_data = [summed_Ex, summed_Ey, summed_Ez, summed_Hx, summed_Hy, summed_Hz];
    save(strcat('.\', save_path, '\summed', string(i_sim), '.mat'), 'save_data');
end

%% BIORADAR SIGNAL PROCESSING
frequencies_baseband = cell2mat(frequencies(:, 1)) - carrier_frequency;

% Load data from the saved intermediate results
data = cell(num_of_simulations, 4);
for i_sim = 1:num_of_simulations
    loaded_data = load(strcat(save_path, '\summed', string(i_sim), '.mat'));
    loaded_data = loaded_data.save_data;
    loaded_data = loaded_data(:, 3); % For Ez component
    % Time domain data
    data{i_sim, 1} = loaded_data;
    % Frequency domain data
    data{i_sim, 2} = fftshift(fft(loaded_data)) / datapoints;
    % iFFT to spatial domain
    data{i_sim, 3} = zeros(32, 1); % frequency components
    data{i_sim, 4} = zeros(32, 1); % distances
    for i_frequency_component = 1:32 % samples where the different frequencies are located
        frequency_component = data{i_sim, 2}(97+(i_frequency_component-1)*128);
        data{i_sim, 3}(i_frequency_component) = frequency_component;
        data{i_sim, 4} = data{i_sim, 4} + abs(frequency_component) * ...
            exp(1j*2*pi/32*(i_frequency_component-1).*(0:31)');
    end
end

%%
% Construct one cycle, only first 17 samples contain all data
spatial_domain_data = zeros(17, 2*num_of_simulations-2);

for i_sim = 1:num_of_simulations
    spatial_domain_data(:, i_sim) = real(data{i_sim, 4}(1:17));
end
for i_sim = (num_of_simulations+1):(2*(num_of_simulations)-2)
    spatial_domain_data(:, i_sim) = real(data{2*num_of_simulations - i_sim, 4}(1:17));
end

N_frames = 20e3;
%repetions = ceil(N_frames / (2*simulations-2));
repetitions = 3;
frames_data = repmat(spatial_domain_data, 1, repetitions);
save(strcat(pwd, '\', save_path, '\frames_data', '.mat'), 'frames_data');

%% DATA PLOTTING OF THE INDIVIDUAL SIGNAL PROCESSING STEPS
close all
f_fft = f_s_frame/datapoints*(-datapoints/2:datapoints/2-1);
figure()
hold on
grid on
% Quarter frame
plot((0:1027)/f_s_frame*1e6, real(data{1, 1}(1:1028))/1e3)
plot((0:1027)/f_s_frame*1e6, imag(data{1, 1}(1:1028))/1e3)
legend('Real part', 'Imaginary part', 'FontName', 'CMU Serif')
axis padded; xlim('tight')
xlabel('Time in µs')
ylabel('$E_\mathrm{z}$ in kV/m', 'interpreter', 'latex')
set(gca, 'FontSize', 16)
exportgraphics(gcf,strcat(image_path, simulation_directory, ...
    '_time_domain.pdf'), 'ContentType','vector')

figure()
hold on
grid on
plot(f_fft/1e6, real(data{1, 2}))
plot(f_fft/1e6, imag(data{1, 2}))
legend('Real part', 'Imaginary part', 'FontName', 'CMU Serif')
axis padded; xlim('tight')
xlabel('Frequency in MHz', 'interpreter', 'latex')
ylabel('$E_\mathrm{z}$ in V/m', 'interpreter', 'latex')
set(gca, 'FontSize', 16)
exportgraphics(gcf,strcat(image_path, simulation_directory, ...
    '_frequency_domain1.pdf'), 'ContentType','vector')

% To only show the 32 relevant components
figure()
hold on
grid on
f_fft_2 = linspace(f_fft(97), f_fft(4065), length(data{1, 3}));
plot(f_fft_2/1e6, real(data{1, 3}), '-x')
plot(f_fft_2/1e6, imag(data{1, 3}), '-x')
plot(f_fft_2/1e6, abs(data{1, 3}), '-x')
legend('Real part', 'Imaginary part', 'Magnitude', 'FontName', 'CMU Serif')
axis padded; xlim('tight')
y_lim = ylim;
ylim([y_lim(1) y_lim(2) + 400])
clear y_lim
xlabel('Frequency in MHz', 'interpreter', 'latex')
ylabel('$E_\mathrm{z}$ in V/m', 'interpreter', 'latex')
set(gca, 'FontSize', 16)
exportgraphics(gcf,strcat(image_path, simulation_directory, ...
    '_frequency_domain2.pdf'), 'ContentType','vector')

% Spatial domain
figure()
hold on
grid on
plot((1:31)*delta_R, real(data{1, 4}(2:end)), '-x')
plot((1:31)*delta_R, imag(data{1, 4}(2:end)), '-x')
xline(16*delta_R,'-',{'Axis of','symmetry'}, ...
    'LabelHorizontalAlignment', 'left', 'FontSize', 16, 'FontName', 'CMU Serif');
legend('Real part', 'Imaginary part', 'FontName', 'CMU Serif')
axis padded
xlabel('Range in m', 'interpreter', 'latex')
ylabel('IFFT of radar frequency domain', 'interpreter', 'latex')
set(gca, 'FontSize', 16)
exportgraphics(gcf,strcat(image_path, simulation_directory, ...
    '_spatial_domain.pdf'), 'ContentType','vector')

% Range bin plots
figure()
tl = tiledlayout(3, 1);
for i = 2:4
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
exportgraphics(gcf,strcat(image_path, simulation_directory, ...
    '_range_frames.pdf'), 'ContentType','vector')


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
exportgraphics(gcf,strcat(image_path, simulation_directory, ...
    '_range_estimate.pdf'), 'ContentType','vector')

% Calculate prominence ratio PR and save it to .txt-file
PR_squared = range_activity.^2;
PR = (max(PR_squared)*(length(PR_squared)-1))/(sum(PR_squared) - max(PR_squared));

% Write prominence ratio to file
txt_file = fopen('PR.txt', 'wt');
fprintf(txt_file, 'Prominence ratio: %.4f', PR);
fclose(txt_file);