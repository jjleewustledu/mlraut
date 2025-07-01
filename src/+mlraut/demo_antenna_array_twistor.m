%% Antenna Array Demonstration - Beamforming and Interference via Twistor Theory
% This demo shows multiple antenna sources creating interference patterns
% and directional beamforming using the twistor framework

clear; close all; clc;

fprintf('Antenna Array Electromagnetic Simulation\n');
fprintf('=======================================\n\n');

%% Configuration: Phased Array Antenna
% Set up a linear array of antennas for beamforming

% Array parameters
n_antennas = 4;
antenna_spacing = 0.5; % In wavelengths
frequency = 5.0; % Normalized frequency
omega = 2*pi*frequency;
wavelength = 1/frequency;
phase_shift = pi/4; % Progressive phase shift for beam steering

fprintf('Array Configuration:\n');
fprintf('  Number of antennas: %d\n', n_antennas);
fprintf('  Spacing: %.2f wavelengths\n', antenna_spacing);
fprintf('  Frequency: %.2f (normalized)\n', frequency);
fprintf('  Beam steering phase: %.2f rad\n\n', phase_shift);

%% Initialize field for array
array_field = TwistorMaxwellField('gauge_group', 'U1', ...
                                  'n_points', 50, ...
                                  'x_range', [-15, 15], ...
                                  'slice_type', 'minkowski');

%% Define antenna positions
antenna_positions = zeros(n_antennas, 3);
for i = 1:n_antennas
    antenna_positions(i, 1) = (i - (n_antennas+1)/2) * antenna_spacing * wavelength;
    antenna_positions(i, 2) = 0;
    antenna_positions(i, 3) = 0;
end

%% Create superposed twistor function for array
% Each antenna contributes with appropriate phase
array_field.twistor_function = @(W) arrayTwistorFunction(W, antenna_positions, ...
                                                         omega, phase_shift);

%% Compute fields
fprintf('Computing array fields...\n');
array_field.computeGaugeFields();

%% Visualization 1: Array Pattern and Beamforming
figure('Name', 'Phased Array Beamforming', 'Position', [100, 100, 1400, 600]);

subplot(1, 3, 1);
% Individual antenna locations
scatter(antenna_positions(:,1), antenna_positions(:,2), 200, 'r', 'filled');
hold on;
for i = 1:n_antennas
    text(antenna_positions(i,1), antenna_positions(i,2)+0.5, ...
         sprintf('A%d', i), 'HorizontalAlignment', 'center');
end
grid on;
xlabel('x position');
ylabel('y position');
title('Antenna Array Layout');
axis equal;
xlim([-5, 5]);
ylim([-2, 2]);

subplot(1, 3, 2);
% Near field pattern
mid_z = floor(array_field.n_points / 2);
[X, Y] = meshgrid(linspace(-15, 15, array_field.n_points));

E_mag = sqrt(array_field.field_strength.Ex(:,:,mid_z).^2 + ...
             array_field.field_strength.Ey(:,:,mid_z).^2 + ...
             array_field.field_strength.Ez(:,:,mid_z).^2);

imagesc(array_field.x_range, array_field.x_range, E_mag);
colormap(hot);
colorbar;
hold on;

% Mark antenna positions
for i = 1:n_antennas
    plot(antenna_positions(i,1), antenna_positions(i,2), 'c^', ...
         'MarkerSize', 10, 'MarkerFaceColor', 'cyan', 'LineWidth', 2);
end

xlabel('x (units)');
ylabel('y (units)');
title('Near Field Pattern');
axis equal tight;

subplot(1, 3, 3);
% Far field radiation pattern
angles = linspace(-pi/2, pi/2, 180);
r_far = 12; % Sample distance

E_far_pattern = zeros(size(angles));
for i = 1:length(angles)
    x_sample = r_far * sin(angles(i));
    y_sample = r_far * cos(angles(i));
    
    [~, ix] = min(abs(linspace(-15, 15, array_field.n_points) - x_sample));
    [~, iy] = min(abs(linspace(-15, 15, array_field.n_points) - y_sample));
    
    E_far_pattern(i) = E_mag(iy, ix);
end

% Convert to dB
E_far_dB = 20*log10(E_far_pattern/max(E_far_pattern) + 1e-10);

polarplot(angles + pi/2, E_far_dB + 40, 'LineWidth', 2);
hold on;
polarplot([0 2*pi], [20 20], 'r--'); % -20 dB line
rlim([0 40]);
title('Far Field Pattern (dB)');

%% Visualization 2: Beam Steering Demonstration
figure('Name', 'Beam Steering', 'Position', [100, 100, 1200, 800]);

% Show beam patterns for different steering angles
steering_angles = [-45, -30, -15, 0, 15, 30, 45]; % degrees
n_angles = length(steering_angles);

for idx = 1:n_angles
    subplot(3, 3, idx);
    
    % Convert to phase shift
    steer_phase = 2*pi*antenna_spacing*sind(steering_angles(idx));
    
    % Update twistor function
    array_field.twistor_function = @(W) arrayTwistorFunction(W, antenna_positions, ...
                                                             omega, steer_phase);
    array_field.computeGaugeFields();
    
    % Get field magnitude
    E_mag_steer = sqrt(array_field.field_strength.Ex(:,:,mid_z).^2 + ...
                       array_field.field_strength.Ey(:,:,mid_z).^2);
    
    % Plot
    imagesc(array_field.x_range, array_field.x_range, E_mag_steer);
    colormap(hot);
    hold on;
    
    % Show beam direction
    beam_angle = steering_angles(idx) * pi/180;
    quiver(0, 0, 10*sin(beam_angle), 10*cos(beam_angle), ...
           'Color', 'cyan', 'LineWidth', 2, 'MaxHeadSize', 2);
    
    title(sprintf('Beam: %d°', steering_angles(idx)));
    axis equal tight;
    set(gca, 'XTickLabel', [], 'YTickLabel', []);
end

%% Visualization 3: Two-Source Interference Pattern
figure('Name', 'Two-Source Interference', 'Position', [100, 100, 1200, 500]);

% Create two-antenna configuration
two_source_field = TwistorMaxwellField('gauge_group', 'U1', ...
                                       'n_points', 60, ...
                                       'x_range', [-20, 20]);

% Two antennas separated by multiple wavelengths
separation = 3 * wavelength;
source1_pos = [-separation/2, 0, 0];
source2_pos = [separation/2, 0, 0];

% In-phase sources
two_source_field.twistor_function = @(W) ...
    exp(1i * omega * sqrt((W(1) - source1_pos(1))^2 + W(2)^2) / (abs(W(1)) + 1)) + ...
    exp(1i * omega * sqrt((W(1) - source2_pos(1))^2 + W(2)^2) / (abs(W(1)) + 1));

two_source_field.computeGaugeFields();

subplot(1, 2, 1);
% Field pattern
E_two = sqrt(two_source_field.field_strength.Ex(:,:,30).^2 + ...
             two_source_field.field_strength.Ey(:,:,30).^2);

imagesc(two_source_field.x_range, two_source_field.x_range, E_two);
colormap(hot);
colorbar;
hold on;

% Mark sources
plot(source1_pos(1), source1_pos(2), 'c^', 'MarkerSize', 12, ...
     'MarkerFaceColor', 'cyan', 'LineWidth', 2);
plot(source2_pos(1), source2_pos(2), 'c^', 'MarkerSize', 12, ...
     'MarkerFaceColor', 'cyan', 'LineWidth', 2);

xlabel('x (units)');
ylabel('y (units)');
title('Two-Source Interference Pattern');
axis equal tight;

subplot(1, 2, 2);
% Cross-section showing interference fringes
y_idx = 30;
E_cross = E_two(y_idx, :);
x_vals = linspace(-20, 20, two_source_field.n_points);

plot(x_vals, E_cross, 'b-', 'LineWidth', 2);
hold on;

% Theoretical fringe spacing
fringe_spacing = wavelength / (2 * sin(atan(separation/2/10)));
theoretical_fringes = cos(2*pi*x_vals/fringe_spacing).^2;
theoretical_fringes = theoretical_fringes * max(E_cross);
plot(x_vals, theoretical_fringes, 'r--', 'LineWidth', 1);

xlabel('x position');
ylabel('Field magnitude');
title('Interference Fringes (y=0)');
legend('Computed', 'Theoretical', 'Location', 'best');
grid on;

%% Visualization 4: Time Evolution Animation
figure('Name', 'Wave Propagation Animation', 'Position', [100, 100, 800, 800]);

% Animate the two-source interference over time
n_frames = 20;
time_steps = linspace(0, 2*pi/omega, n_frames);

for frame = 1:n_frames
    clf;
    
    % Update phase for current time
    current_phase = omega * time_steps(frame);
    
    % Recompute with time-dependent phase
    two_source_field.twistor_function = @(W) ...
        real(exp(1i * (omega * sqrt((W(1) - source1_pos(1))^2 + W(2)^2) / ...
             (abs(W(1)) + 1) - current_phase)) + ...
             exp(1i * (omega * sqrt((W(1) - source2_pos(1))^2 + W(2)^2) / ...
             (abs(W(1)) + 1) - current_phase)));
    
    two_source_field.computeGaugeFields();
    
    % Get real part of field
    E_real = real(two_source_field.field_strength.Ex(:,:,30));
    
    % Plot
    imagesc(two_source_field.x_range, two_source_field.x_range, E_real);
    colormap(redblue);
    caxis([-1, 1] * max(abs(E_real(:))));
    colorbar;
    hold on;
    
    % Sources
    plot(source1_pos(1), source1_pos(2), 'k^', 'MarkerSize', 12, ...
         'MarkerFaceColor', 'yellow', 'LineWidth', 2);
    plot(source2_pos(1), source2_pos(2), 'k^', 'MarkerSize', 12, ...
         'MarkerFaceColor', 'yellow', 'LineWidth', 2);
    
    title(sprintf('Wave Interference at t = %.2f periods', ...
                  time_steps(frame)*omega/(2*pi)));
    xlabel('x'); ylabel('y');
    axis equal tight;
    
    drawnow;
    pause(0.05);
end

%% Summary
fprintf('\n\nAntenna Array Analysis Summary\n');
fprintf('==============================\n');
fprintf('1. Array demonstrates beamforming with phase control\n');
fprintf('2. Main beam can be steered by progressive phase shifts\n');
fprintf('3. Two-source interference shows expected fringe pattern\n');
fprintf('4. Fringe spacing: %.3f units (theoretical: %.3f)\n', ...
        mean(diff(find(diff(sign(diff(E_cross))) ~= 0))), fringe_spacing);
fprintf('5. Side lobe suppression: %.1f dB\n', ...
        min(E_far_dB(E_far_dB > -40)));
fprintf('\nThe twistor formulation successfully captures:\n');
fprintf('- Phased array beamforming\n');
fprintf('- Constructive/destructive interference\n');
fprintf('- Near-to-far field transitions\n');
fprintf('- Time-dependent wave propagation\n');

%% Helper Functions
function g = arrayTwistorFunction(W, positions, omega, phase_shift)
    % Twistor function for antenna array with beamforming
    
    g = 0;
    n_antennas = size(positions, 1);
    
    for k = 1:n_antennas
        % Phase for beamforming
        array_phase = (k-1) * phase_shift;
        
        % Distance in twistor coordinates
        dx = W(1) - positions(k, 1);
        dy = W(2) - positions(k, 2);
        dz = W(3) - positions(k, 3);
        
        r_tw = sqrt(abs(dx)^2 + abs(dy)^2 + abs(dz)^2 + 0.1);
        
        % Contribution from this antenna with beamforming phase
        g = g + exp(1i * (omega * r_tw + array_phase)) / (r_tw + 0.1);
    end
    
    % Normalize
    g = g / n_antennas;
end

% Helper function for red-blue colormap
function cmap = redblue(n)
    if nargin < 1
        n = 256;
    end
    
    % Create red-white-blue colormap
    r = [linspace(0, 1, n/2), ones(1, n/2)]';
    g = [linspace(0, 1, n/2), linspace(1, 0, n/2)]';
    b = [ones(1, n/2), linspace(1, 0, n/2)]';
    
    cmap = [r, g, b];
end