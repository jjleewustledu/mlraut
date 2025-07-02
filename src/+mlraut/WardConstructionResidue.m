classdef WardConstructionResidue < handle
    % Ward construction using residue theorem for proper pole handling
    
    properties
        epsilon
        scale_factor
    end
    
    methods
        function obj = WardConstructionResidue()
            obj.epsilon = 1e-10;
            obj.scale_factor = 1.0;

            warning("mlraut:DeprecationError", ...
                "mapping from twistor space residues to spacetime gauge potentials requires the proper index structure");
        end
        
        function A_mu = computeMonopoleGauge(obj, x, pole_location)
            % Compute gauge potential for monopole using residue method
            % pole_location is in twistor space: Y = (Y_A, Y^{A'})
            
            % Convert to spinor notation
            x_spinor = obj.spacetimeToSpinor(x);
            
            % Find poles of <Z,Y> = 0 as function of ζ
            % For π' = (π₀', π₁') with different parametrizations
            
            % Try parametrization 1: π' = (1, ζ)
            [poles1, residues1] = obj.findPoles_param1(x_spinor, pole_location);
            
            % Try parametrization 2: π' = (ζ, 1)  
            [poles2, residues2] = obj.findPoles_param2(x_spinor, pole_location);
            
            % Initialize gauge potential
            A_mu = zeros(4, 1);
            
            % Sum residue contributions
            for k = 1:length(poles1)
                if abs(poles1(k)) < 1/obj.epsilon  % Pole in patch 1
                    A_mu = A_mu + obj.residueContribution(x_spinor, poles1(k), residues1(k), 1);
                end
            end
            
            for k = 1:length(poles2)
                if abs(poles2(k)) > obj.epsilon  % Pole in patch 2
                    A_mu = A_mu + obj.residueContribution(x_spinor, poles2(k), residues2(k), 2);
                end
            end
            
            % Apply physical normalization
            A_mu = real(A_mu) * obj.scale_factor;
        end
        
        function [poles, residues] = findPoles_param1(obj, x_spinor, Y)
            % Find poles for parametrization π' = (1, ζ)
            % Solve <Z,Y> = 0 where Z = (ix^{AA'}π_{A'}, π^{A'})
            
            % <Z,Y> = ω^A Y_A - Y^{A'} π_{A'}
            % With ω^A = ix^{AA'}π_{A'} and π' = (1, ζ):
            % <Z,Y> = i(x₀₀' + x₀₁'ζ)Y₀ + i(x₁₀' + x₁₁'ζ)Y₁ - Y⁰' - Y¹'ζ
            
            % Rearrange as aζ + b = 0
            a = 1i*x_spinor(1,2)*Y(1) + 1i*x_spinor(2,2)*Y(2) - Y(4);
            b = 1i*x_spinor(1,1)*Y(1) + 1i*x_spinor(2,1)*Y(2) - Y(3);
            
            if abs(a) > obj.epsilon
                poles = -b/a;
                
                % Residue of 1/<Z,Y> at the pole
                % Res = lim_{ζ→pole} (ζ-pole)/<Z,Y>
                % For simple pole: Res = 1/(d<Z,Y>/dζ) = 1/a
                residues = 1/a;
            else
                poles = [];
                residues = [];
            end
        end
        
        function [poles, residues] = findPoles_param2(obj, x_spinor, Y)
            % Find poles for parametrization π' = (ζ, 1)
            
            % Similar calculation with roles swapped
            a = 1i*x_spinor(1,1)*Y(1) + 1i*x_spinor(2,1)*Y(2) - Y(3);
            b = 1i*x_spinor(1,2)*Y(1) + 1i*x_spinor(2,2)*Y(2) - Y(4);
            
            if abs(a) > obj.epsilon
                poles = -b/a;
                residues = 1/a;
            else
                poles = [];
                residues = [];
            end
        end
        
        function A_contrib = residueContribution(obj, x_spinor, pole, residue, patch)
            % Calculate gauge potential contribution from a residue
            
            % At pole ζ₀, the residue theorem gives:
            % A_μ = (i/2π) × 2πi × Res[π_{A'} ∂log(g)/∂ζ dx^μ/dω^A]
            
            % Construct π_{A'} at the pole
            if patch == 1
                pi_at_pole = [1; pole];
            else
                pi_at_pole = [pole; 1];
            end
            
            % Normalize
            pi_norm = sqrt(abs(pi_at_pole(1))^2 + abs(pi_at_pole(2))^2);
            pi_normalized = pi_at_pole / pi_norm;
            
            % The gauge potential components come from:
            % A_μ dx^μ = -i π_{A'} d(log g)
            
            % For monopole g = 1/<Z,Y>, we have d(log g) = -d<Z,Y>/<Z,Y>
            % The residue theorem gives the coefficient
            
            A_contrib = zeros(4, 1);
            
            % Components related to π_{A'} structure
            % This is the key formula from Ward's construction
            factor = -residue;  % includes 2πi from residue theorem
            
            % Euclidean signature mapping
            A_contrib(1) = real(factor * pi_normalized(1) * conj(pi_normalized(1)));  % A_t
            A_contrib(2) = real(factor * pi_normalized(2) * conj(pi_normalized(1)));  % A_x  
            A_contrib(3) = imag(factor * pi_normalized(2) * conj(pi_normalized(1)));  % A_y
            A_contrib(4) = real(factor * pi_normalized(2) * conj(pi_normalized(2)));  % A_z
        end
        
        function x_spinor = spacetimeToSpinor(obj, x)
            % Convert spacetime to spinor (same as parent)
            x_spinor = [x(1) + x(2), x(3) + 1i*x(4);
                       x(3) - 1i*x(4), x(1) - x(2)] / sqrt(2);
        end
        
        function A_analytic = analyticMonopole(obj, x)
            % Direct analytic formula for comparison
            r = sqrt(x(2)^2 + x(3)^2 + x(4)^2);
            
            if r > obj.epsilon
                A_analytic = zeros(4, 1);
                A_analytic(2) = -x(3) / (r * (r + x(4)));  % A_x
                A_analytic(3) = x(2) / (r * (r + x(4)));   % A_y
            else
                A_analytic = zeros(4, 1);
            end
        end
    end
end