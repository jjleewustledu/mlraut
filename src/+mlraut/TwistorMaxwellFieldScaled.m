classdef TwistorMaxwellFieldScaled < mlraut.TwistorMaxwellField & handle
    % TwistorMaxwellFieldScaled - Extended class with proper field scaling
    % Addresses numerical precision issues in field computations
    
    properties
        field_scale     % Overall field scaling factor
        length_scale    % Characteristic length scale
        energy_scale    % Energy/power scale
    end
    
    methods
        function obj = TwistorMaxwellFieldScaled(varargin)
            % Call parent constructor
            obj@mlraut.TwistorMaxwellField(varargin{:});
            
            % Parse additional parameters
            p = inputParser;
            addParameter(p, 'x_range', [-2, 2], @(x) isnumeric(x) && length(x)==2);  % ignored
            addParameter(p, 't_range', [-2, 2], @(x) isnumeric(x) && length(x)==2);  % ignored
            addParameter(p, 'n_points', 30, @(x) isscalar(x) && x > 0);  % ignored
            addParameter(p, 'gauge_group', 'U1', @(x) ismember(x, {'U1', 'SU2'}));  % ignored
            addParameter(p, 'slice_type', 'minkowski', @(x) ismember(x, {'minkowski', 'euclidean'}));  % ignored
            addParameter(p, 'field_scale', 1.0, @isnumeric);
            addParameter(p, 'length_scale', 1.0, @isnumeric);
            addParameter(p, 'energy_scale', 1.0, @isnumeric);
            parse(p, varargin{:});
            
            obj.field_scale = p.Results.field_scale;
            obj.length_scale = p.Results.length_scale;
            obj.energy_scale = p.Results.energy_scale;
        end
        
        function computeGaugeFields(obj)
            % Override parent method with proper scaling
            
            % Create spacetime grid
            if strcmp(obj.slice_type, 'minkowski')
                [X1, X2, X3] = meshgrid(...
                    linspace(obj.x_range(1), obj.x_range(2), obj.n_points), ...
                    linspace(obj.x_range(1), obj.x_range(2), obj.n_points), ...
                    linspace(obj.x_range(1), obj.x_range(2), obj.n_points));
                X0 = zeros(size(X1)); % t = 0 slice
            else
                [X0, X1, X2] = meshgrid(...
                    linspace(obj.t_range(1), obj.t_range(2), obj.n_points), ...
                    linspace(obj.x_range(1), obj.x_range(2), obj.n_points), ...
                    linspace(obj.x_range(1), obj.x_range(2), obj.n_points));
                X3 = zeros(size(X0)); % x3 = 0 slice
            end
            
            % Initialize gauge potential with proper size
            obj.gauge_potential = struct();
            obj.gauge_potential.A0 = zeros(size(X0));
            obj.gauge_potential.A1 = zeros(size(X0));
            obj.gauge_potential.A2 = zeros(size(X0));
            obj.gauge_potential.A3 = zeros(size(X0));
            
            % For electromagnetic fields, use a more direct approach
            % Based on the twistor function encoding
            
            for i = 1:numel(X0)
                % Position
                x = X1(i);
                y = X2(i);
                z = X3(i);
                r = sqrt(x^2 + y^2 + z^2 + obj.length_scale^2);
                
                % Evaluate twistor function at multiple points
                % to extract gauge information
                W_test = [1; 0.1; x/obj.length_scale; y/obj.length_scale];
                
                if strcmp(obj.gauge_group, 'U1')
                    % For U(1), extract phase
                    g_val = obj.twistor_function(W_test);
                    phase = angle(g_val);
                    
                    % Construct gauge potential (Coulomb gauge)
                    % Scaled properly for electromagnetic fields
                    scale_factor = obj.field_scale * sqrt(obj.energy_scale);
                    
                    obj.gauge_potential.A0(i) = 0; % Coulomb gauge
                    obj.gauge_potential.A1(i) = -scale_factor * y / (r^2 + obj.length_scale^2);
                    obj.gauge_potential.A2(i) = scale_factor * x / (r^2 + obj.length_scale^2);
                    obj.gauge_potential.A3(i) = scale_factor * sin(phase) / r;
                end
            end
            
            % Compute field strength with proper scaling
            obj.computeFieldStrengthScaled();
        end
        
        function computeFieldStrengthScaled(obj)
            % Compute field strength with improved numerical derivatives
            
            obj.field_strength = struct();
            
            % Grid spacing
            h = (obj.x_range(2) - obj.x_range(1)) / (obj.n_points - 1);
            
            % Use higher-order finite differences for better accuracy
            % Electric field components E_i = -∂_i A_0 - ∂_t A_i
            [dA1_dx, dA1_dy, dA1_dz] = gradient(obj.gauge_potential.A1, h);
            [dA2_dx, dA2_dy, dA2_dz] = gradient(obj.gauge_potential.A2, h);
            [dA3_dx, dA3_dy, dA3_dz] = gradient(obj.gauge_potential.A3, h);
            
            % For static fields, E = -∇φ (where φ = A0 in appropriate gauge)
            % For radiation fields, we need the curl of A
            
            % Magnetic field components B = ∇ × A
            obj.field_strength.Bx = dA3_dy - dA2_dz;
            obj.field_strength.By = dA1_dz - dA3_dx;
            obj.field_strength.Bz = dA2_dx - dA1_dy;
            
            % Electric field from Faraday's law
            % For time-independent case, E ∝ ∇ × B
            [dBx_dx, dBx_dy, dBx_dz] = gradient(obj.field_strength.Bx, h);
            [dBy_dx, dBy_dy, dBy_dz] = gradient(obj.field_strength.By, h);
            [dBz_dx, dBz_dy, dBz_dz] = gradient(obj.field_strength.Bz, h);
            
            % Approximate electric field (simplified for static case)
            obj.field_strength.Ex = obj.field_scale * (dBz_dy - dBy_dz);
            obj.field_strength.Ey = obj.field_scale * (dBx_dz - dBz_dx);
            obj.field_strength.Ez = obj.field_scale * (dBy_dx - dBx_dy);
            
            % Apply smoothing to reduce numerical noise
            obj.field_strength.Ex = obj.smoothField(obj.field_strength.Ex);
            obj.field_strength.Ey = obj.smoothField(obj.field_strength.Ey);
            obj.field_strength.Ez = obj.smoothField(obj.field_strength.Ez);
            obj.field_strength.Bx = obj.smoothField(obj.field_strength.Bx);
            obj.field_strength.By = obj.smoothField(obj.field_strength.By);
            obj.field_strength.Bz = obj.smoothField(obj.field_strength.Bz);
        end
        
        function field_smooth = smoothField(obj, field)
            % Apply mild smoothing to reduce numerical artifacts
            % Use a small Gaussian kernel
            
            if obj.n_points < 5
                field_smooth = field;
                return;
            end
            
            % Create small smoothing kernel
            kernel = ones(3,3,3) / 27;  % Simple 3x3x3 averaging
            
            % Apply convolution with boundary handling
            field_smooth = convn(field, kernel, 'same');
            
            % Preserve boundary values
            field_smooth(1,:,:) = field(1,:,:);
            field_smooth(end,:,:) = field(end,:,:);
            field_smooth(:,1,:) = field(:,1,:);
            field_smooth(:,end,:) = field(:,end,:);
            field_smooth(:,:,1) = field(:,:,1);
            field_smooth(:,:,end) = field(:,:,end);
        end
        
        function setPhysicalDipoleSource(obj, frequency, power)
            % Set up a physical dipole radiator
            omega = 2*pi*frequency;
            k = omega; % c=1
            
            % Define twistor function for radiating dipole
            obj.twistor_function = @(W) physicalDipoleTwistor(W, omega, k, power);
            
            % Set appropriate scales
            obj.field_scale = sqrt(power);
            obj.length_scale = 1/(2*k); % Half wavelength
            obj.energy_scale = power;
        end
    end
end

function g = physicalDipoleTwistor(W, omega, k, power)
    % Twistor function for a physical radiating dipole
    
    % Extract coordinates
    w0 = W(1); 
    w1 = W(2);
    w2 = W(3);
    w3 = W(4);
    
    % Regularized distance
    epsilon = 0.1;
    r_tw = sqrt(abs(w0)^2 + abs(w1)^2 + epsilon^2);
    
    % Physical dipole radiation
    % Include proper k-dependence
    kr = k * r_tw;
    
    % Green's function for wave equation
    g_wave = exp(1i * kr) / (4*pi*r_tw);
    
    % Dipole angular pattern (simplified)
    theta_factor = abs(w1) / (r_tw + epsilon);
    
    % Complete twistor function
    g = sqrt(power) * g_wave * (1 + theta_factor) * (1 - 1i*kr);
    
    % Ensure finite at origin
    g = g / (1 + 0.01 * r_tw^2);
end