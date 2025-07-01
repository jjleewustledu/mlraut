classdef RiemannSphere < handle
    % RIEMANNSPHERE Class for plotting complex time-series on the Riemann sphere
    %
    % The RiemannSphere class provides visualization of complex-valued time-series
    % data mapped onto the Riemann sphere using stereographic projection.
    %
    % KEY FEATURES:
    %   - Stereographic Projection: Maps complex numbers to points on unit sphere
    %   - Handles Infinity: Maps infinite values to the north pole of the sphere  
    %   - Time-based Coloring: Colors trajectory based on time progression
    %   - Customizable Visualization: Semi-transparent sphere, equator/meridians
    %   - Start/End Markers: Clear indication of trajectory beginning and end
    %   - Multiple Examples: Built-in demo with various trajectory types
    %
    % STEREOGRAPHIC PROJECTION FORMULA:
    %   For complex number z = x + iy:
    %   X = 2x/(1 + |z|²)
    %   Y = 2y/(1 + |z|²)  
    %   Z = (|z|² - 1)/(1 + |z|²)
    %   Infinite values map to north pole (0, 0, 1)
    %
    % PROPERTIES:
    %   z        - Complex time-series data (vector)
    %   t        - Time vector (defaults to 1:length(z))
    %   options  - Structure with visualization options
    %
    % METHODS:
    %   RiemannSphere(z, t, options) - Constructor
    %   plot()                       - Create Riemann sphere visualization
    %   demo()                       - Run demonstration with example trajectories
    %   setOptions(options)          - Update visualization options
    %   getProjectedCoordinates()    - Get 3D coordinates from stereographic projection
    %
    % VISUALIZATION OPTIONS:
    %   sphere_alpha     - Transparency of sphere surface (default: 0.3)
    %   show_equator     - Display equator line (default: true)
    %   show_meridians   - Display meridian lines (default: true)
    %   colormap_name    - Colormap for time-based coloring (default: 'viridis')
    %   line_width       - Width of trajectory line (default: 2)
    %   point_size       - Size of trajectory points (default: 30)
    %
    % EXAMPLE:
    %   t = linspace(0, 4*pi, 200);
    %   z = exp(1i*t) + 0.5*exp(2i*t);
    %   rs = RiemannSphere(z, t);
    %   rs.plot();
    %
    %   % With custom options
    %   options.colormap_name = 'jet';
    %   options.line_width = 3;
    %   rs = RiemannSphere(z, t, options);
    %   rs.plot();
    %
    %   % Run demo
    %   RiemannSphere.demo();
    %
    % See also: STEREOGRAPHIC_PROJECTION, COMPLEX, SPHERE, PLOT3
    
    properties (Access = public)
        z        % Complex time-series data
        t        % Time vector
        options  % Visualization options structure
    end
    
    properties (Access = private)
        default_options  % Default visualization options
    end
    
    methods (Access = public)
        
        function obj = RiemannSphere(z, t, options)
            % RIEMANNSPHERE Constructor for RiemannSphere class
            %
            % Syntax:
            %   rs = RiemannSphere(z)
            %   rs = RiemannSphere(z, t)
            %   rs = RiemannSphere(z, t, options)
            %
            % Inputs:
            %   z       - Complex time-series data (vector)
            %   t       - Time vector (optional, defaults to 1:length(z))
            %   options - Structure with visualization options (optional)
            %
            % Output:
            %   obj - RiemannSphere object
            
            if nargin < 1
                error('RiemannSphere:InvalidInput', 'Complex data z is required');
            end
            
            % Set default options
            obj.default_options = struct();
            obj.default_options.sphere_alpha = 0.9;
            obj.default_options.show_equator = false;
            obj.default_options.show_meridians = false;
            obj.default_options.colormap_name = 'viridis';
            obj.default_options.line_width = 2;
            obj.default_options.point_size = 30;
            
            % Validate and set complex data
            obj.z = z(:);  % Ensure column vector
            
            % Set time vector
            if nargin < 2 || isempty(t)
                obj.t = (1:length(obj.z))';
            else
                obj.t = t(:);  % Ensure column vector
                if length(obj.t) ~= length(obj.z)
                    error('RiemannSphere:SizeMismatch', ...
                          'Time vector t must have same length as complex data z');
                end
            end
            
            % Set options
            if nargin < 3 || isempty(options)
                obj.options = obj.default_options;
            else
                obj.setOptions(options);
            end
        end
        
        function setOptions(obj, new_options)
            % SETOPTIONS Update visualization options
            %
            % Syntax:
            %   rs.setOptions(new_options)
            %
            % Input:
            %   new_options - Structure with new option values
            %
            % Updates existing options with new values, keeping defaults
            % for unspecified fields.
            
            obj.options = obj.default_options;
            
            if isstruct(new_options)
                fields = fieldnames(new_options);
                for i = 1:length(fields)
                    if isfield(obj.default_options, fields{i})
                        obj.options.(fields{i}) = new_options.(fields{i});
                    else
                        warning('RiemannSphere:UnknownOption', ...
                               'Unknown option: %s', fields{i});
                    end
                end
            end
        end
        
        function [X, Y, Z] = getProjectedCoordinates(obj)
            % GETPROJECTEDCOORDINATES Get 3D coordinates from stereographic projection
            %
            % Syntax:
            %   [X, Y, Z] = rs.getProjectedCoordinates()
            %
            % Outputs:
            %   X, Y, Z - 3D coordinates on unit sphere
            %
            % Maps complex numbers to Riemann sphere using stereographic projection.
            % Infinite values are mapped to the north pole (0, 0, 1).
            
            % Handle infinite values (map to north pole)
            is_finite = isfinite(obj.z);
            z_finite = obj.z(is_finite);
            
            % Compute stereographic projection for finite values
            z_abs_sq = abs(z_finite).^2;
            denom = 1 + z_abs_sq;
            
            X_finite = 2 * real(z_finite) ./ denom;
            Y_finite = 2 * imag(z_finite) ./ denom;
            Z_finite = (z_abs_sq - 1) ./ denom;
            
            % Initialize full coordinate arrays
            X = zeros(size(obj.z));
            Y = zeros(size(obj.z));
            Z = zeros(size(obj.z));
            
            % Fill finite values
            X(is_finite) = X_finite;
            Y(is_finite) = Y_finite;
            Z(is_finite) = Z_finite;
            
            % Infinite values go to north pole
            X(~is_finite) = 0;
            Y(~is_finite) = 0;
            Z(~is_finite) = 1;
        end
        
        function fig_handle = plot(obj, fig_handle)
            % PLOT Create Riemann sphere visualization
            %
            % Syntax:
            %   rs.plot()
            %   rs.plot(fig_handle)
            %   fig_handle = rs.plot()
            %
            % Input:
            %   fig_handle - Optional figure handle (creates new if not provided)
            %
            % Output:
            %   fig_handle - Handle to figure containing the plot
            %
            % Creates a 3D visualization of the complex time-series on the
            % Riemann sphere with customizable appearance options.
            
            % Create or use existing figure
            if nargin < 2 || isempty(fig_handle)
                fig_handle = figure('Position', [100, 100, 800, 600]);
            else
                figure(fig_handle);
                clf;
            end
            
            % Get projected coordinates
            [X, Y, Z] = obj.getProjectedCoordinates();
            
            % Create unit sphere
            [sphere_x, sphere_y, sphere_z] = sphere(50);
            
            % Plot the sphere
            surf(sphere_x, sphere_y, sphere_z, 'FaceAlpha', obj.options.sphere_alpha, ...
                 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.9]);
            hold on;
            
            % Add equator and meridians if requested
            if obj.options.show_equator
                theta = linspace(0, 2*pi, 100);
                plot3(cos(theta), sin(theta), zeros(size(theta)), 'k--', 'LineWidth', 1);
            end
            
            if obj.options.show_meridians
                phi = linspace(-pi/2, pi/2, 100);
                for theta_mer = [0, pi/2, pi, 3*pi/2]
                    x_mer = cos(phi) * cos(theta_mer);
                    y_mer = cos(phi) * sin(theta_mer);
                    z_mer = sin(phi);
                    plot3(x_mer, y_mer, z_mer, 'k--', 'LineWidth', 0.5, ...
                          'Color', [0.5, 0.5, 0.5]);
                end
            end
            
            % Color the trajectory by time
            if length(obj.t) == length(obj.z)
                % Use time for coloring
                color_data = obj.t;
                colormap(obj.options.colormap_name);
                
                % Plot trajectory as colored line
                for i = 1:length(X)-1
                    line([X(i), X(i+1)], [Y(i), Y(i+1)], [Z(i), Z(i+1)], ...
                         'Color', obj.getColorFromTime(obj.t(i), min(obj.t), max(obj.t)), ...
                         'LineWidth', obj.options.line_width);
                end
                
                % Add colored scatter points
                scatter3(X, Y, Z, obj.options.point_size, color_data, 'filled');
                colorbar;
                if length(obj.t) > 1
                    clabel = colorbar;
                    clabel.Label.String = 'Time';
                end
            else
                % Default coloring
                plot3(X, Y, Z, 'b-', 'LineWidth', obj.options.line_width);
                scatter3(X, Y, Z, obj.options.point_size, 'r', 'filled');
            end
            
            % Mark start and end points
            if length(obj.z) > 1
                scatter3(X(1), Y(1), Z(1), 100, 'g', 'filled', ...
                        'MarkerEdgeColor', 'k', 'LineWidth', 2);
                scatter3(X(end), Y(end), Z(end), 100, 'r', 'filled', ...
                        'MarkerEdgeColor', 'k', 'LineWidth', 2);
                
                % Add text labels
                text(X(1), Y(1), Z(1), '  Start', 'FontSize', 12, 'FontWeight', 'bold');
                text(X(end), Y(end), Z(end), '  End', 'FontSize', 12, 'FontWeight', 'bold');
            end
            
            % Set axis properties
            axis equal;
            grid on;
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            title('Complex Time-Series on Riemann Sphere');
            
            % Set viewing angle
            view(45, 30);
            
            % Add legend
            if length(obj.z) > 1
                legend('Sphere', 'Trajectory', 'Points', 'Start', 'End', 'Location', 'best');
            end
            
            hold off;
        end
        
    end
    
    methods (Access = private)
        
        function color = getColorFromTime(obj, t_val, t_min, t_max)
            % GETCOLORFROMTIME Get color from colormap based on normalized time
            %
            % Private method to map time values to colors using the current colormap
            
            if t_max == t_min
                norm_t = 0.5;
            else
                norm_t = (t_val - t_min) / (t_max - t_min);
            end
            
            % Get colormap
            cmap = colormap(obj.options.colormap_name);
            cmap_size = size(cmap, 1);
            
            % Get color index
            color_idx = max(1, min(cmap_size, round(norm_t * (cmap_size - 1)) + 1));
            color = cmap(color_idx, :);
        end
        
    end
    
    methods (Static)
        
        function demo()
            % DEMO Run demonstration with example trajectories
            %
            % Syntax:
            %   RiemannSphere.demo()
            %
            % Creates a figure with four subplots showing different types
            % of complex trajectories on the Riemann sphere:
            %   1. Spiral trajectory
            %   2. Figure-8 pattern  
            %   3. Complex random walk
            %   4. Trajectory with singularity
            %
            % This static method demonstrates the capabilities of the
            % RiemannSphere class with various interesting trajectories.
            
            fprintf('Riemann Sphere Visualization Demo\n');
            fprintf('=================================\n\n');
            
            % Create figure for demo
            fig = figure('Position', [50, 50, 1200, 900]);
            
            % Example 1: Simple spiral
            fprintf('Example 1: Spiral trajectory\n');
            t1 = linspace(0, 4*pi, 200);
            z1 = (0.1 + 0.8*t1/(4*pi)) .* exp(1i*t1);
            
            subplot(2,2,1);
            rs1 = RiemannSphere(z1, t1);
            rs1.plot(fig);
            title('Spiral Trajectory');
            
            % Example 2: Figure-8 pattern
            fprintf('Example 2: Figure-8 pattern\n');
            t2 = linspace(0, 2*pi, 200);
            z2 = exp(1i*t2) + 0.5*exp(-2i*t2);
            
            subplot(2,2,2);
            rs2 = RiemannSphere(z2, t2);
            rs2.plot(fig);
            title('Figure-8 Pattern');
            
            % Example 3: Random walk
            fprintf('Example 3: Complex random walk\n');
            n_steps = 100;
            rng(42); % For reproducibility
            z3 = cumsum(0.1*(randn(n_steps,1) + 1i*randn(n_steps,1)));
            
            subplot(2,2,3);
            rs3 = RiemannSphere(z3);
            rs3.plot(fig);
            title('Random Walk');
            
            % Example 4: Trajectory approaching infinity
            fprintf('Example 4: Trajectory with singularity\n');
            t4 = linspace(0.1, 2, 100);
            z4 = 1./(t4 - 2.1) + 0.5i;  % Approaches infinity as t approaches 2.1
            
            subplot(2,2,4);
            rs4 = RiemannSphere(z4, t4);
            rs4.plot(fig);
            title('Trajectory with Singularity');
            
            % Add overall title
            sgtitle('RiemannSphere Class Demonstration', 'FontSize', 16, 'FontWeight', 'bold');
            
            fprintf('\nDemo complete! Each subplot shows a different type of complex trajectory.\n');
            fprintf('Try creating your own RiemannSphere objects with different complex data.\n');
        end
        
    end
    
end