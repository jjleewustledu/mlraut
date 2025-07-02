classdef WardConstructionImproved < mlraut.WardConstructionFull & handle
    % Improved Ward construction with better scaling and numerics
    
    properties
        scale_factor    % Global scale for fields
        use_residues    % Whether to use residue theorem
    end
    
    methods
        function obj = WardConstructionImproved(n_points, box_size)
            % Call parent constructor
            obj@mlraut.WardConstructionFull(n_points, box_size);
            
            obj.scale_factor = 1.0;
            obj.use_residues = true;
        end
        
        function A_mu = computeGaugeAtPoint(obj, x, twistor_func, gauge_group)
            % Override with improved computation
            
            % Convert to 2-spinor notation
            x_spinor = obj.spacetimeToSpinor(x);
            
            % Initialize result
            A_mu = zeros(4, 1);
            
            % Get scale from twistor function
            % Sample at a few reference points
            ref_points = [
                1, 0, 1, 0;
                0, 1, 0, 1;
                1, 1, 1, 1;
                1, 0, 0, 1
            ];
            
            g_scales = zeros(size(ref_points, 1), 1);
            for i = 1:size(ref_points, 1)
                Z_ref = ref_points(i, :)' / norm(ref_points(i, :));
                g_scales(i) = abs(twistor_func(Z_ref));
            end
            
            % Use median to avoid outliers
            typical_scale = median(g_scales(g_scales > 0));
            if isempty(typical_scale) || typical_scale < 1e-10
                typical_scale = 1;
            end
            
            % Check if we can use residue theorem
            if obj.use_residues
                % Try to find poles analytically
                poles = obj.findPoles(twistor_func, x_spinor);
                if ~isempty(poles)
                    A_mu = obj.computeViaResidues(x_spinor, twistor_func, poles);
                    return;
                end
            end
            
            % Fall back to contour integration with better numerics
            contributions = zeros(obj.contour_points.n_total, 4);
            
            for k = 1:obj.contour_points.n_total
                zeta = obj.contour_points.zeta(k);
                weight = obj.contour_points.weights(k);
                
                % Improved parametrization
                if abs(zeta) <= 1
                    % Standard patch
                    pi_spinor = [1; zeta];
                    jacobian = 1;
                else
                    % Alternate patch for |ζ| > 1
                    pi_spinor = [1/zeta; 1];
                    jacobian = 1/abs(zeta)^2;
                end
                
                % Don't normalize pi_spinor yet
                pi_norm = sqrt(abs(pi_spinor(1))^2 + abs(pi_spinor(2))^2);
                
                % Incidence relation
                omega_spinor = 1i * x_spinor * pi_spinor;
                Z = [omega_spinor; pi_spinor];
                
                % Evaluate function
                g = twistor_func(Z);
                
                if abs(g) > obj.epsilon * typical_scale
                    % Improved derivative computation
                    dlog_g_dz = obj.computeLogarithmicDerivative(twistor_func, Z, k, zeta, x_spinor);
                    
                    % Ward's formula with proper factors
                    % A_μ = (i/2π) ∮ π_{A'} d(log g)/dζ (∂x^μ/∂ω^A) dζ
                    
                    % The key insight: ∂x^μ/∂ω^A gives the gauge index
                    % For Euclidean signature with our conventions:
                    factor = weight * jacobian / (2*pi * pi_norm^2);
                    
                    % Extract components (this is the delicate index calculation)
                    contributions(k, 1) = factor * real(pi_spinor(1) * conj(dlog_g_dz));
                    contributions(k, 2) = factor * real(pi_spinor(2) * conj(dlog_g_dz));
                    contributions(k, 3) = factor * imag(pi_spinor(1) * conj(dlog_g_dz));
                    contributions(k, 4) = factor * imag(pi_spinor(2) * conj(dlog_g_dz));
                end
            end
            
            % Sum contributions
            A_mu = sum(contributions, 1)';
            
            % Apply physical scaling
            r = norm(x);
            if r > obj.epsilon
                % Scale like 1/r for monopole-like solutions
                A_mu = A_mu * typical_scale * obj.scale_factor / sqrt(1 + r^2);
            else
                A_mu = A_mu * typical_scale * obj.scale_factor;
            end
            
            % Ensure real for Euclidean theory
            A_mu = real(A_mu);
        end
        
        function dlog_g_dz = computeLogarithmicDerivative(obj, twistor_func, Z, k, zeta, x_spinor)
            % Improved computation of d(log g)/dζ
            
            h = 1e-8;  % Smaller step for better accuracy
            
            % Use complex derivative
            % d/dζ = (1/2)(∂/∂ξ - i∂/∂η) where ζ = ξ + iη
            
            % Four-point stencil for better accuracy
            Z_plus = Z;
            Z_minus = Z;
            Z_iplus = Z;
            Z_iminus = Z;
            
            % Vary ζ in the appropriate component
            if abs(zeta) <= 1
                % Varying π' = (1, ζ)
                Z_plus(4) = Z(4) + h;
                Z_minus(4) = Z(4) - h;
                Z_iplus(4) = Z(4) + 1i*h;
                Z_iminus(4) = Z(4) - 1i*h;
                
                % Update ω consistently for shifted π
                if norm(x_spinor) > 0
                    pi_plus = [1; zeta + h];
                    pi_minus = [1; zeta - h];
                    pi_iplus = [1; zeta + 1i*h];
                    pi_iminus = [1; zeta - 1i*h];
                    
                    Z_plus(1:2) = 1i * x_spinor * pi_plus;
                    Z_minus(1:2) = 1i * x_spinor * pi_minus;
                    Z_iplus(1:2) = 1i * x_spinor * pi_iplus;
                    Z_iminus(1:2) = 1i * x_spinor * pi_iminus;
                end
            else
                % Varying π' = (ζ, 1) → (1/ζ, 1)
                % d/dζ acts on 1/ζ component
                inv_zeta = 1/zeta;
                Z_plus(3) = inv_zeta - h/zeta^2;
                Z_minus(3) = inv_zeta + h/zeta^2;
                Z_iplus(3) = inv_zeta - 1i*h/zeta^2;
                Z_iminus(3) = inv_zeta + 1i*h/zeta^2;
                
                % Update ω consistently
                if norm(x_spinor) > 0
                    pi_plus = [inv_zeta - h/zeta^2; 1];
                    pi_minus = [inv_zeta + h/zeta^2; 1];
                    pi_iplus = [inv_zeta - 1i*h/zeta^2; 1];
                    pi_iminus = [inv_zeta + 1i*h/zeta^2; 1];
                    
                    Z_plus(1:2) = 1i * x_spinor * pi_plus;
                    Z_minus(1:2) = 1i * x_spinor * pi_minus;
                    Z_iplus(1:2) = 1i * x_spinor * pi_iplus;
                    Z_iminus(1:2) = 1i * x_spinor * pi_iminus;
                end
            end
            
            % Evaluate function at shifted points
            g = twistor_func(Z);
            g_plus = twistor_func(Z_plus);
            g_minus = twistor_func(Z_minus);
            g_iplus = twistor_func(Z_iplus);
            g_iminus = twistor_func(Z_iminus);
            
            % Complex derivative of log
            if abs(g) > obj.epsilon
                d_real = (log(abs(g_plus)) - log(abs(g_minus))) / (2*h);
                d_imag = (log(abs(g_iplus)) - log(abs(g_iminus))) / (2*h);
                
                % Phase contribution
                phase_diff = angle(g_plus/g_minus);
                if abs(phase_diff) > pi
                    phase_diff = phase_diff - 2*pi*sign(phase_diff);
                end
                d_phase = phase_diff / (2*h);
                
                dlog_g_dz = d_real + 1i*(d_imag + d_phase);
            else
                dlog_g_dz = 0;
            end
            
            % Account for coordinate transformation if needed
            if abs(zeta) > 1
                dlog_g_dz = -dlog_g_dz / zeta^2;
            end
        end
        
        function poles = findPoles(obj, twistor_func, x_spinor)
            % Try to identify poles of the twistor function
            % This is problem-specific
            
            poles = [];  % Empty for now - could implement pole detection
            
            % For known functions, we could hard-code pole locations
            % e.g., for g(Z) = 1/(Z·pole + ε), extract pole
        end
        
        function A_mu = computeViaResidues(obj, x_spinor, twistor_func, poles)
            % Use residue theorem for known poles
            
            A_mu = zeros(4, 1);
            
            % For each pole, compute residue contribution
            for p = 1:size(poles, 1)
                pole = poles(p, :)';
                
                % Residue of π_{A'} d(log g)/dζ at the pole
                % This requires knowing the order and residue
                
                % Simplified for simple poles
                residue = 2*pi*1i;  % Standard simple pole residue
                
                % Contribution to gauge potential
                % Details depend on the specific pole structure
            end
        end
    end
end