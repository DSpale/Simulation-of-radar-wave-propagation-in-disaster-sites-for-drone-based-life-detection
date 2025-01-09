clear
close all

% Adjust image path
image_path = 'D:\Uni\10. Semester\Masterarbeit\images\MCCW\';
if ~exist(image_path, 'dir')
    mkdir(image_path)
end
[~, simulation_directory, ~] = fileparts(pwd);

set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Times New Roman');

space_res = 5e-4;

%% Get metadata of the measurement
% MACROS
Ex = 2;
Ey = 3;
Ez = 4;
Hx = 5;
Hy = 6;
Hz = 7;

% Read the simulation data filenames
simulation_files = dir(fullfile(pwd, '*.out'));

num_of_simulations = length(simulation_files);
data = cell(1, 7);
% Data format: Info, Ex-Ez, Hx-Hz
for i = 1:num_of_simulations
    data_name = strcat(pwd, '\', 'MCCW_radar', string(i), '.out');
    data{i, 1} = h5info(data_name);
    data{i, 2} = h5read(data_name, '/rxs/rx1/Ex');
    data{i, 3} = h5read(data_name, '/rxs/rx1/Ey');
    data{i, 4} = h5read(data_name, '/rxs/rx1/Ez');
    data{i, 5} = h5read(data_name, '/rxs/rx1/Hx');
    data{i, 6} = h5read(data_name, '/rxs/rx1/Hy');
    data{i, 7} = h5read(data_name, '/rxs/rx1/Hz');
end

% Determine number of iterations, time step and simulated time window
iterations = double(data{1,1}.Attributes(3).Value);
dt = data{1,1}.Attributes(6).Value;
tsim = (iterations - 1) * dt;

%% Plotting
zero_values = false; % Plot only non-zero field quantities
num_shown_displacements = 4; % Displacements depicted in field quantities plot
% First row displacements: Index of simulation file
% Second row displacements: Physical displacements in mm
displacements = round(linspace(1, 101, num_shown_displacements));
displacements = [displacements; (displacements-1)*space_res*1e3];
t = linspace(0, tsim, iterations);

% Plot the field quantities
figure('Position', [100, 100, 600, 600])
if zero_values
    tiledlayout(3, 2)
else
    tiledlayout(3, 1)
end

if zero_values
    nexttile
    title('$E_\mathrm{x}$ component', 'interpreter', 'latex')
    hold on
    grid on
    box on
    for i = 1:num_shown_displacements
        plot(t*1e9, data{displacements(1, i), Ex} / 1000)
    end
    axis padded; xlim('tight')
    ylabel('kV/m')
    set(gca, 'FontSize', 14)
end

nexttile
title('$H_\mathrm{x}$ component', 'interpreter', 'latex')
hold on
grid on
box on
for i = 1:num_shown_displacements
    plot(t*1e9, data{displacements(1, i), Hx})
end
axis padded; xlim('tight')
ylabel('A/m', 'interpreter', 'latex')
set(gca, 'FontSize', 14)

if zero_values
    nexttile
    title('$E_\mathrm{y}$ component', 'interpreter', 'latex')
    hold on
    grid on
    box on
    for i = 1:num_shown_displacements
        plot(t*1e9, data{displacements(1, i), Ey} / 1000)
    end
    axis padded; xlim('tight')
    ylabel('kV/m', 'interpreter', 'latex')
    set(gca, 'FontSize', 14)
end

nexttile
title('$H_\mathrm{y}$ component', 'interpreter', 'latex')
hold on
grid on
box on
for i = 1:num_shown_displacements
    plot(t*1e9, data{displacements(1, i), Hy})
end
axis padded; xlim('tight')
ylabel('A/m', 'interpreter', 'latex')
set(gca, 'FontSize', 14)

nexttile
title('$E_\mathrm{z}$ component', 'interpreter', 'latex')
hold on
grid on
box on
for i = 1:num_shown_displacements
    plot(t*1e9, data{displacements(1, i), Ez} / 1000)
end
axis padded; xlim('tight')
xlabel('Time in ns', 'interpreter', 'latex')
ylabel('kV/m', 'interpreter', 'latex')
set(gca, 'FontSize', 14)

if zero_values
    nexttile
    title('$H_\mathrm{z}$ component', 'interpreter', 'latex')
    hold on
    grid on
    box on
    for i = 1:num_shown_displacements
        plot(t*1e9, data{displacements(1, i), Hz})
    end
    axis padded; xlim('tight')
    xlabel('Time in s', 'interpreter', 'latex')
    ylabel('A/m', 'interpreter', 'latex')
    i
end

radii_label = strcat(string(displacements(2,:)), ' mm');
lg = legend(radii_label, 'Orientation', 'horizontal', 'NumColumns', 2);
lg.Layout.Tile = 'south';
exportgraphics(gcf,strcat(image_path, simulation_directory, ...
    '_field_quantities.pdf'), 'ContentType','vector')


%% Generate figure depicting phase change due to breathing
excitation_frequency = 1e9;
t_corr = 0:dt:(1/excitation_frequency);
cos_corr = cos(2*pi*excitation_frequency*t_corr);

% Correlate all the data
% Layout of data:
%   Row 1: Truncated data
%   Row 2: Correlated values
%   Row 3: Phase with maximum correlation aka shift

data_trunc_corr = cell(num_of_simulations, 2);
starting_time = 18e-9;
for i = 1:num_of_simulations
    % Regard the data from 18 ns onwards
    data_trunc_corr{i, 1} = data{i, Ez}(round(starting_time/dt):end);
end

% Determine the shift for maximal correlation
for i = 1:num_of_simulations
    [c, lags] = xcorr(data_trunc_corr{i, 1}, cos_corr);
    data_trunc_corr{i, 2} = c;
    [~, I_prox] = max(c);
    lags_max = lags(I_prox);
    data_trunc_corr{i, 3} = mod(lags_max, ...
        length(t_corr))/length(t_corr)*360;
end
clear c I_prox lags_max

% Plot the phase shift and amplitude as determined by cross-correlation
figure
hold on
grid on
box on
phases = cell2mat(data_trunc_corr(:, 3));
amplitudes = cellfun(@max, data_trunc_corr(:,1));
displacements = (0:(num_of_simulations-1))*space_res;
scatter(displacements*1e2, phases, 'x')
xlabel('Displacement in cm', 'interpreter', 'latex')
ylabel('Phase in Degrees', 'interpreter', 'latex')
axis padded
set(gca, 'FontSize', 26)
exportgraphics(gcf,strcat(image_path, simulation_directory, ...
    '_phases.pdf'), 'ContentType','vector')

figure
hold on
grid on
box on
plot(displacements*1e2, amplitudes / 1000, 'x')
xlabel('Displacement in cm', 'interpreter', 'latex')
ylabel('Amplitude in kV/m', 'interpreter', 'latex')
axis padded
set(gca, 'FontSize', 26)
exportgraphics(gcf,strcat(image_path, simulation_directory, ...
    '_amplitudes.pdf'), 'ContentType','vector')

x_comps = zeros(1, num_of_simulations);
y_comps = zeros(1, num_of_simulations);
for i = 1:num_of_simulations
    x_comps(i) = amplitudes(i)*cos(phases(i)/360*2*pi);
    y_comps(i) = amplitudes(i)*sin(phases(i)/360*2*pi);
end

% Fit circle to data with Levenberg-Marquardt-optimization
[fitted_center, fitted_radius] = levenberg_marquardt_circle(x_comps, y_comps);

% Display the results
center_magnitude = sqrt(fitted_center(1)^2+fitted_center(2)^2);
phase1 = atan2(y_comps(end)-fitted_center(2), x_comps(end)-fitted_center(1));
phase2 = atan2(y_comps(1)-fitted_center(2), x_comps(1)-fitted_center(1));
fprintf('Fitted center: (%.4f, %.4f), Magnitude: %.4f\n', ...
    fitted_center(1), fitted_center(2), center_magnitude);
fprintf('Fitted radius: %.4f\n', fitted_radius);
fprintf('Phase difference: %.4f°\n', abs(phase1-phase2)/(2*pi)*360);

% Write results to text file
txt_file = fopen('circle_results.txt', 'wt');
fprintf(txt_file, 'Fitted center: (%.4f, %.4f), Magnitude: %.4f\n', ...
    fitted_center(1), fitted_center(2), center_magnitude);
fprintf(txt_file, 'Fitted radius: %.4f\n', fitted_radius);
fprintf(txt_file, 'Phase difference: %.4f°\n', abs(phase1-phase2)/(2*pi)*360);
fclose(txt_file);

% Plot the results
figure
scatter(x_comps, y_comps, 'DisplayName', 'Data Points');
hold on;
scatter(fitted_center(1), fitted_center(2), 'r', 'filled', ...
    'DisplayName', 'Static center')
theta_circle = linspace(0, 2*pi, 100);
circle_x = fitted_center(1) + fitted_radius * cos(theta_circle);
circle_y = fitted_center(2) + fitted_radius * sin(theta_circle);
plot(circle_x, circle_y, 'r--', 'LineWidth', 1, 'DisplayName', 'Fitted Circle');
xlabel('Real $E_\mathrm{z}$ in V/m', 'interpreter', 'latex');
ylabel('Imaginary $E_\mathrm{z}$ in V/m', 'interpreter', 'latex');
axis padded
axis equal
grid on
box on
set(gca, 'FontSize', 14)
exportgraphics(gcf,strcat(image_path, simulation_directory, ...
    '_circle.pdf'), 'ContentType','vector')

%% COSINE FITTING
t = ((0:(length(cell2mat(data_trunc_corr(1,1)))-1))*dt)';

% Solve with linear regression
% A*cos(x - y) = A*cos(x)cos(y) + A*sin(x)sin(y)

for ii = 1:num_of_simulations
    data_LVO = cell2mat(data_trunc_corr(ii,1));
    B = [arrayfun(@(x) cos(x), 2*pi*1e9*t) arrayfun(@(x) sin(x), 2*pi*1e9*t)];
    A = pinv(B' * B) * B' * data_LVO;
    disp(strcat("Simulation ", string(ii), " optimized"))
    amplitudes_LinReg(ii) = sqrt(A(1)^2 + A(2)^2);
    phases_LinReg(ii) = mod(angle(A(1) + 1j*A(2))/(2*pi)*360, 360);
end

% Representation in degrees
phases = phases / (2*pi)*360;

% Plot the results
figure
hold on
grid on
box on
displacements = (0:(num_of_simulations-1))*space_res;
scatter(displacements*1e2, phases_LinReg, 'x')
xlabel('Displacement in cm', 'interpreter', 'latex')
ylabel('Phase in Degrees', 'interpreter', 'latex')
axis padded
set(gca, 'FontSize', 26)
exportgraphics(gcf,strcat(image_path, simulation_directory, ...
    '_phases_LinReg.pdf'), 'ContentType','vector')

figure
hold on
grid on
box on
plot(displacements*1e2, amplitudes_LinReg / 1000, 'x')
xlabel('Displacement in cm', 'interpreter', 'latex')
ylabel('Amplitude in kV/m', 'interpreter', 'latex')
axis padded
set(gca, 'FontSize', 26)
exportgraphics(gcf,strcat(image_path, simulation_directory, ...
    '_amplitudes_LinReg.pdf'), 'ContentType','vector')

x_comps = zeros(1, num_of_simulations);
y_comps = zeros(1, num_of_simulations);
for i = 1:num_of_simulations
    x_comps(i) = amplitudes_LinReg(i)*cos(phases_LinReg(i)/360*2*pi);
    y_comps(i) = amplitudes_LinReg(i)*sin(phases_LinReg(i)/360*2*pi);
end

% Construct circle
[fitted_center, fitted_radius] = levenberg_marquardt_circle(x_comps, y_comps);

% Display the results
center_magnitude = sqrt(fitted_center(1)^2+fitted_center(2)^2);
phase1 = atan2(y_comps(end)-fitted_center(2), x_comps(end)-fitted_center(1));
phase2 = atan2(y_comps(1)-fitted_center(2), x_comps(1)-fitted_center(1));
fprintf('Fitted center: (%.4f, %.4f), Magnitude: %.4f\n', ...
    fitted_center(1), fitted_center(2), center_magnitude);
fprintf('Fitted radius: %.4f\n', fitted_radius);
fprintf('Phase difference: %.4f°\n', abs(phase1-phase2)/(2*pi)*360);

% Write once again to text file
txt_file = fopen('circle_results_LinReg.txt', 'wt');
fprintf(txt_file, 'Fitted center: (%.4f, %.4f), Magnitude: %.4f\n', ...
    fitted_center(1), fitted_center(2), center_magnitude);
fprintf(txt_file, 'Fitted radius: %.4f\n', fitted_radius);
fprintf(txt_file, 'Phase difference: %.4f°\n', abs(phase1-phase2)/(2*pi)*360);
fclose(txt_file);

% Plot the results
figure
scatter(x_comps, y_comps, 'DisplayName', 'Data Points');
hold on;
scatter(fitted_center(1), fitted_center(2), 'r', 'filled', ...
    'DisplayName', 'Static center')
theta_circle = linspace(0, 2*pi, 100);
circle_x = fitted_center(1) + fitted_radius * cos(theta_circle);
circle_y = fitted_center(2) + fitted_radius * sin(theta_circle);
plot(circle_x, circle_y, 'r--', 'LineWidth', 1, 'DisplayName', 'Fitted Circle');
xlabel('Real $E_\mathrm{z}$ in V/m', 'interpreter', 'latex');
ylabel('Imaginary $E_\mathrm{z}$ in V/m', 'interpreter', 'latex');
axis padded
axis equal
grid on
box on
set(gca, 'FontSize', 14)
exportgraphics(gcf,strcat(image_path, simulation_directory, ...
    '_circle_LinReg.pdf'), 'ContentType','vector')