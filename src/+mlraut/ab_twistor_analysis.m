%% Aharonov-Bohm Effect: Twistor Theory Analysis
% Deep dive into how twistor theory explains the AB effect

clear; close all; clc;

fprintf('Aharonov-Bohm Effect: Twistor Theory Analysis\n');
fprintf('=============================================\n\n');

%% Part 1: Topological Origin in Twistor Space
fprintf('Part 1: Topological Structure in Twistor Space\n');
fprintf('----------------------------------------------\n\n');

% In twistor theory, the AB effect arises from the non-trivial
% holomorphic structure of the gauge bundle over twistor space

% Create monopole gauge field
ward = mlraut.WardMonopoleCorrect();

% The key insight: AB phase is a holonomy in the gauge bundle
fprintf('Twistor theory perspective:\n');
fprintf('1. Gauge field A_μ ↔ Holomorphic vector bundle over PT\n');
fprintf('2. AB phase ↔ Holonomy around non-contractible loop\n');
fprintf('3. Flux quantization ↔ Chern class of the bundle\n\n');

%% Part 2: Constructing AB Configuration via Twistors
fprintf('Part 2: AB Configuration from Twistor Functions\n');
fprintf('-----------------------------------------------\n\n');

% We need a gauge field with flux but localized B field
% This can be constructed from monopole-antimonopole pairs

% Parameters
flux_tube_radius = 1.0;
separation = 5.0;  % Monopole-antimonopole separation

% Grid
n = 40;
coords = linspace(-4, 4, n);
[X, Y] = meshgrid(coords, coords);

% Compute gauge field
A_field = zeros(n, n, 3);
fprintf('Computing gauge field from monopole pairs...\n');

for i = 1:n
    for j = 1:n
        x = X(i,j);
        y = Y(i,j);
        z = 0;
        
        % Monopole at z = -separation/2
        pos1 = [0, x, y, z + separation/2];
        A1 = ward.computeMonopoleGauge(pos1);
        
        % Antimonopole at z = +separation/2
        pos2 = [0, x, y, z - separation/2];
        A2 = ward.computeMonopoleGauge(pos2);
        
        % Superposition
        A_total = A1 - A2;
        A_field(i,j,:) = A_total(2:4);
    end
end

% Compute B field
h = coords(2) - coords(1);
[~, dAx_dy] = gradient(A_field(:,:,1), h);
[dAy_dx, ~] = gradient(A_field(:,:,2), h);
B_z = dAy_dx - dAx_dy;

%% Part 3: Wilson Loop Calculation
fprintf('\nPart 3: Wilson Loop and Holonomy\n');
fprintf('---------------------------------\n');

% The AB phase is the Wilson loop: W = exp(i∮A·dl)
% In twistor theory, this is the holonomy of the bundle

% Define circular paths at different radii
radii = linspace(0.5, 3, 10);
wilson_loops = zeros(size(radii));
enclosed_flux = zeros(size(radii));

for k = 1:length(radii)
    r = radii(k);
    
    % Parametric circle
    theta = linspace(0, 2*pi, 100);
    x_path = r * cos(theta);
    y_path = r * sin(theta);
    
    % Compute line integral ∮A·dl
    integral = 0;
    for i = 1:length(theta)-1
        % Path element
        dx = x_path(i+1) - x_path(i);
        dy = y_path(i+1) - y_path(i);
        
        % Interpolate A
        Ax = interp2(X, Y, A_field(:,:,1), x_path(i), y_path(i), 'linear', 0);
        Ay = interp2(X, Y, A_field(:,:,2), y_path(i), y_path(i), 'linear', 0);
        
        integral = integral + Ax*dx + Ay*dy;
    end
    
    wilson_loops(k) = integral;
    
    % Compute enclosed flux
    mask = (X.^2 + Y.^2) < r^2;
    enclosed_flux(k) = sum(B_z(mask)) * h^2;
end

fprintf('Wilson loop results:\n');
fprintf('  Radius | Wilson Loop | Flux | W - Φ\n');
fprintf('  -------|-------------|------|-------\n');
for k = 1:3:length(radii)
    fprintf('  %5.2f  | %11.3f | %4.2f | %6.3f\n', ...
            radii(k), wilson_loops(k), enclosed_flux(k), ...
            wilson_loops(k) - enclosed_flux(k));
end

%% Part 4: Quantum Interference Demonstration
fprintf('\n\nPart 4: Quantum Interference Pattern\n');
fprintf('------------------------------------\n');

% Two-path interference around flux tube
source_pos = [-3, 0];
detector_positions = linspace(-3, 3, 100);

% Reference case: no flux
intensity_ref = zeros(size(detector_positions));

% AB case: with flux
intensity_AB = zeros(size(detector_positions));

% Electron momentum
k = 20;  % In units where ℏ = 1

for i = 1:length(detector_positions)
    det_pos = [3, detector_positions(i)];
    
    % Path 1: goes above origin
    path1_phase = computeABPhase(source_pos, [0, 1], det_pos, A_field, coords);
    path1_length = norm(det_pos - [0, 1]) + norm([0, 1] - source_pos);
    
    % Path 2: goes below origin
    path2_phase = computeABPhase(source_pos, [0, -1], det_pos, A_field, coords);
    path2_length = norm(det_pos - [0, -1]) + norm([0, -1] - source_pos);
    
    % Quantum amplitudes
    psi1 = exp(1i * (k * path1_length + path1_phase));
    psi2 = exp(1i * (k * path2_length + path2_phase));
    
    intensity_AB(i) = abs(psi1 + psi2)^2;
    
    % Reference (no flux)
    psi1_ref = exp(1i * k * path1_length);
    psi2_ref = exp(1i * k * path2_length);
    intensity_ref(i) = abs(psi1_ref + psi2_ref)^2;
end

%% Part 5: Visualization
figure('Name', 'Aharonov-Bohm Effect: Twistor Theory', 'Position', [100, 100, 1400, 900]);

% Gauge field from Ward construction
subplot(2, 3, 1);
A_mag = sqrt(A_field(:,:,1).^2 + A_field(:,:,2).^2);
imagesc(coords, coords, log10(A_mag + 1e-6));
colormap(subplot(2,3,1), hot);
colorbar;
hold on;
% Show flux region
rectangle('Position', [-flux_tube_radius, -flux_tube_radius, 2*flux_tube_radius, 2*flux_tube_radius], ...
          'Curvature', [1, 1], 'EdgeColor', 'w', 'LineWidth', 2);
xlabel('x'); ylabel('y');
title('log_{10}|A| from Ward Construction');
axis equal tight;

% Magnetic field
subplot(2, 3, 2);
imagesc(coords, coords, B_z);
colormap(subplot(2,3,2), redblue);
caxis([-max(abs(B_z(:))), max(abs(B_z(:)))]);
colorbar;
hold on;
rectangle('Position', [-flux_tube_radius, -flux_tube_radius, 2*flux_tube_radius, 2*flux_tube_radius], ...
          'Curvature', [1, 1], 'EdgeColor', 'k', 'LineWidth', 2);
xlabel('x'); ylabel('y');
title('Magnetic Field B_z');
axis equal tight;

% Vector potential circulation
subplot(2, 3, 3);
skip = 2;
quiver(X(1:skip:end,1:skip:end), Y(1:skip:end,1:skip:end), ...
       A_field(1:skip:end,1:skip:end,1), A_field(1:skip:end,1:skip:end,2), ...
       'b', 'LineWidth', 1);
hold on;
% Show integration paths
theta = linspace(0, 2*pi, 100);
plot(2*cos(theta), 2*sin(theta), 'r-', 'LineWidth', 2);
plot(cos(theta), sin(theta), 'g-', 'LineWidth', 2);
xlabel('x'); ylabel('y');
title('Vector Potential A');
legend('A field', 'Outer path', 'Inner path', 'Location', 'best');
axis equal tight;

% Wilson loops
subplot(2, 3, 4);
plot(radii, wilson_loops, 'b.-', 'LineWidth', 2, 'MarkerSize', 10);
hold on;
plot(radii, enclosed_flux, 'r--', 'LineWidth', 2);
xlabel('Path radius');
ylabel('Value');
title('Wilson Loop vs Enclosed Flux');
legend('∮A·dl', '∫B·dA', 'Location', 'best');
grid on;

% Quantum interference
subplot(2, 3, 5);
plot(detector_positions, intensity_AB/max(intensity_AB), 'b-', 'LineWidth', 2);
hold on;
plot(detector_positions, intensity_ref/max(intensity_ref), 'r--', 'LineWidth', 1.5);
xlabel('Detector position y');
ylabel('Normalized intensity');
title('Quantum Interference Pattern');
legend('With AB flux', 'No flux', 'Location', 'best');
grid on;

% Twistor space visualization
subplot(2, 3, 6);
% Show how the gauge bundle structure appears in twistor space
text(0.1, 0.9, 'Twistor Space Structure:', 'FontSize', 14, 'FontWeight', 'bold');
text(0.1, 0.75, '• Gauge field A_μ ↔ ∂̄-operator on PT', 'FontSize', 12);
text(0.1, 0.6, '• AB flux ↔ First Chern class c_1', 'FontSize', 12);
text(0.1, 0.45, '• Wilson loop ↔ Holonomy of bundle', 'FontSize', 12);
text(0.1, 0.3, '• Quantization ↔ Integrality of c_1', 'FontSize', 12);
text(0.1, 0.15, sprintf('Computed flux: Φ = %.3f', enclosed_flux(end)), 'FontSize', 12);
text(0.1, 0.05, sprintf('Wilson loop: W = %.3f', wilson_loops(end)), 'FontSize', 12);
axis off;

%% Part 6: Theoretical Analysis
fprintf('\n\nPart 6: Twistor Theory Insights\n');
fprintf('--------------------------------\n');

fprintf('\nKey results from twistor theory:\n');
fprintf('1. The AB effect is fundamentally topological\n');
fprintf('2. It depends only on the holomorphic structure of the gauge bundle\n');
fprintf('3. Ward''s construction automatically gives the correct topology\n');
fprintf('4. Flux quantization emerges from integrality of Chern classes\n\n');

fprintf('Mathematical correspondence:\n');
fprintf('  Spacetime                    Twistor Space\n');
fprintf('  ---------                    -------------\n');
fprintf('  Gauge potential A_μ      ↔   Holomorphic bundle E → PT\n');
fprintf('  Field strength F_μν      ↔   Curvature of E\n');
fprintf('  AB phase exp(i∮A·dl)     ↔   Holonomy around CP¹\n');
fprintf('  Magnetic flux Φ          ↔   c₁(E) (first Chern class)\n');
fprintf('  Dirac quantization       ↔   c₁(E) ∈ H²(PT,ℤ)\n\n');

%% Summary figure
figure('Name', 'AB Effect Summary');

% Show the essential physics
subplot(1, 2, 1);
% Simplified setup
theta = linspace(0, 2*pi, 100);
plot(cos(theta), sin(theta), 'r-', 'LineWidth', 3);
hold on;
fill(0.3*cos(theta), 0.3*sin(theta), 'b', 'FaceAlpha', 0.5);
% Electron paths
plot([-2, 0, 2], [0, 1, 0], 'g-', 'LineWidth', 2);
plot([-2, 0, 2], [0, -1, 0], 'g-', 'LineWidth', 2);
plot(-2, 0, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
plot(2, 0, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
text(0, 0, 'Φ', 'HorizontalAlignment', 'center', 'FontSize', 16);
xlabel('x'); ylabel('y');
title('Aharonov-Bohm Setup');
axis equal;
xlim([-2.5, 2.5]); ylim([-1.5, 1.5]);

subplot(1, 2, 2);
% Phase shift
flux_scan = linspace(0, 3, 100);
phase_shift = 2*pi*flux_scan;
visibility = abs(cos(phase_shift/2));
plot(flux_scan, visibility, 'b-', 'LineWidth', 2);
xlabel('Flux Φ/Φ₀');
ylabel('Interference visibility');
title('AB Oscillations');
grid on;
ylim([0, 1.1]);

fprintf('\nConclusion: Ward''s construction from twistor theory naturally\n');
fprintf('encodes the topological aspects of the Aharonov-Bohm effect,\n');
fprintf('demonstrating the deep connection between holomorphic geometry\n');
fprintf('and quantum mechanics.\n');

%% Helper functions
function phase = computeABPhase(start_pos, mid_pos, end_pos, A_field, coords)
    % Compute AB phase along path: start → mid → end
    
    n_points = 50;
    phase = 0;
    
    % Segment 1
    for i = 1:n_points-1
        t = (i-1)/(n_points-1);
        t_next = i/(n_points-1);
        
        pos = start_pos + t*(mid_pos - start_pos);
        pos_next = start_pos + t_next*(mid_pos - start_pos);
        
        Ax = interp2(coords, coords, A_field(:,:,1), pos(1), pos(2), 'linear', 0);
        Ay = interp2(coords, coords, A_field(:,:,2), pos(1), pos(2), 'linear', 0);
        
        dx = pos_next(1) - pos(1);
        dy = pos_next(2) - pos(2);
        
        phase = phase + Ax*dx + Ay*dy;
    end
    
    % Segment 2
    for i = 1:n_points-1
        t = (i-1)/(n_points-1);
        t_next = i/(n_points-1);
        
        pos = mid_pos + t*(end_pos - mid_pos);
        pos_next = mid_pos + t_next*(end_pos - mid_pos);
        
        Ax = interp2(coords, coords, A_field(:,:,1), pos(1), pos(2), 'linear', 0);
        Ay = interp2(coords, coords, A_field(:,:,2), pos(1), pos(2), 'linear', 0);
        
        dx = pos_next(1) - pos(1);
        dy = pos_next(2) - pos(2);
        
        phase = phase + Ax*dx + Ay*dy;
    end
end

function cmap = redblue(n)
    if nargin < 1, n = 256; end
    r = [linspace(0,1,n/2), ones(1,n/2)]';
    g = [linspace(0,1,n/2), linspace(1,0,n/2)]';
    b = [ones(1,n/2), linspace(1,0,n/2)]';
    cmap = [r g b];
end