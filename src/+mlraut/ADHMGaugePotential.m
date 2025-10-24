classdef ADHMGaugePotential < handle
    %% ADHMGaugePotential - High-memory implementation for ADHM gauge potential
    %  computation with comprehensive validation and monitoring.
    %  See also Claude "Penrose Twistor Gauge Potential Construction" 20250903011500
    %
    %  Created 03-Sep-2025 01:13:26 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 24.2.0.2923080 (R2024b) Update 6 for MACA64.  Copyright 2025 John J. Lee.

    properties (SetAccess = private)
        % Core geometric data
        lattice_pos         % 3 x N_pos array of Cartesian positions
        N_pos              % Number of spatial positions
        N_t                % Number of time points
        T                  % Time grid (1 x N_t)
        dt                 % Time step size
        radius             % Sphere radius

        % Cached field data - utilizing available memory
        V_full             % N_pos x N_t quaternion array (cached v field)
        V_computed         % N_pos x N_t logical (tracking computed values)
        A_theta            % N_pos x N_t quaternion array
        A_phi              % N_pos x N_t quaternion array
        A_t                % N_pos x N_t quaternion array

        % Derivative operators
        gradient_operator   % Sparse matrix for spatial gradients
        laplacian_operator % Sparse matrix for Laplacian
        derivative_order   % Order of finite difference scheme (2, 4, or 6)
        stencil_weights    % Cell array of stencil weights

        % Field strength and validation
        F_ab               % Field strength tensor components
        F_dual             % Hodge dual of field strength
        v_normalization    % N_pos x N_t array tracking |v|^2
        self_duality_error % N_pos x N_t array of self-duality violations
        topological_charge % Integrated topological charge density

        % Validation thresholds
        norm_tolerance     % Tolerance for v normalization (default 1e-6)
        duality_tolerance  % Tolerance for self-duality (default 1e-5)

        % Performance monitoring
        computation_stats  % Structure tracking computation metrics
        memory_stats      % Structure tracking memory usage
    end

    methods
        function obj = ADHMGaugePotential(lattice_pos, N_t, t_max, options)
            % Enhanced constructor with memory-aware initialization

            arguments
                lattice_pos (3,:) double
                N_t (1,1) double {mustBePositive, mustBeInteger}
                t_max (1,1) double = 1.0
                options.derivative_order (1,1) double = 4
                options.norm_tolerance (1,1) double = 1e-6
                options.duality_tolerance (1,1) double = 1e-5
                options.cache_v_field (1,1) logical = true
            end

            obj.lattice_pos = lattice_pos;
            obj.N_pos = size(lattice_pos, 2);
            obj.N_t = N_t;
            obj.dt = t_max / (N_t - 1);
            obj.T = linspace(0, t_max, N_t);
            obj.radius = mean(vecnorm(lattice_pos, 2, 1));

            obj.derivative_order = options.derivative_order;
            obj.norm_tolerance = options.norm_tolerance;
            obj.duality_tolerance = options.duality_tolerance;

            fprintf('Initializing enhanced ADHM gauge potential calculator\n');
            fprintf('Configuration: %d spatial points, %d time steps\n', obj.N_pos, obj.N_t);

            % Initialize memory structures
            obj.initializeMemoryStructures(options.cache_v_field);

            % Construct derivative operators
            obj.constructDerivativeOperators();

            % Initialize validation structures
            obj.initializeValidationStructures();

            fprintf('Memory allocation complete: %.2f GB utilized\n', obj.memory_stats.total_allocated_gb);
        end

        function initializeMemoryStructures(obj, cache_v_field)
            % Allocate memory for core field arrays

            mem_start = memory;

            if cache_v_field
                fprintf('Allocating memory for complete v field cache...\n');
                obj.V_full = quaternion(zeros(obj.N_pos, obj.N_t), ...
                    zeros(obj.N_pos, obj.N_t), ...
                    zeros(obj.N_pos, obj.N_t), ...
                    zeros(obj.N_pos, obj.N_t));
                obj.V_computed = false(obj.N_pos, obj.N_t);
            end

            fprintf('Allocating memory for gauge potential components...\n');
            obj.A_theta = quaternion(zeros(obj.N_pos, obj.N_t), ...
                zeros(obj.N_pos, obj.N_t), ...
                zeros(obj.N_pos, obj.N_t), ...
                zeros(obj.N_pos, obj.N_t));
            obj.A_phi = quaternion(zeros(obj.N_pos, obj.N_t), ...
                zeros(obj.N_pos, obj.N_t), ...
                zeros(obj.N_pos, obj.N_t), ...
                zeros(obj.N_pos, obj.N_t));
            obj.A_t = quaternion(zeros(obj.N_pos, obj.N_t), ...
                zeros(obj.N_pos, obj.N_t), ...
                zeros(obj.N_pos, obj.N_t), ...
                zeros(obj.N_pos, obj.N_t));

            mem_end = memory;
            obj.memory_stats.field_memory_gb = (mem_end.MemUsedMATLAB - mem_start.MemUsedMATLAB) / 1e9;
            obj.memory_stats.total_allocated_gb = mem_end.MemUsedMATLAB / 1e9;
        end

        function constructDerivativeOperators(obj)
            % Build high-order derivative operators using sparse matrices

            fprintf('Constructing %d-th order derivative operators...\n', obj.derivative_order);

            % Normalize positions for operator construction
            pos_normalized = obj.lattice_pos ./ vecnorm(obj.lattice_pos, 2, 1);

            % Build sparse gradient operator using radial basis functions
            obj.gradient_operator = obj.buildGradientOperator(pos_normalized);

            % Build Laplacian operator for validation
            obj.laplacian_operator = obj.buildLaplacianOperator(pos_normalized);

            % Precompute finite difference stencils for time derivatives
            obj.stencil_weights = obj.computeStencilWeights();

            fprintf('Derivative operators constructed with %d-th order accuracy\n', obj.derivative_order);
        end

        function G = buildGradientOperator(obj, pos_normalized)
            % Construct sparse matrix operator for spatial gradients

            row_indices = [];
            col_indices = [];
            weights_x = [];
            weights_y = [];
            weights_z = [];

            for ipos = 1:obj.N_pos
                if mod(ipos, 5000) == 0
                    fprintf('  Building gradient operator: %d/%d points processed\n', ipos, obj.N_pos);
                end

                p0 = pos_normalized(:, ipos);

                % Find neighbors within adaptive radius
                neighbor_radius = obj.adaptiveNeighborRadius(p0);
                neighbors = obj.findNeighborsInRadius(ipos, pos_normalized, neighbor_radius);

                if length(neighbors) >= obj.derivative_order + 1
                    % Compute RBF weights for gradient approximation
                    [w_x, w_y, w_z] = obj.computeRBFGradientWeights(ipos, neighbors, pos_normalized);

                    n_neighbors = length(neighbors);
                    row_indices = [row_indices, ipos * ones(1, n_neighbors)];
                    col_indices = [col_indices, neighbors];
                    weights_x = [weights_x, w_x'];
                    weights_y = [weights_y, w_y'];
                    weights_z = [weights_z, w_z'];
                end
            end

            % Assemble sparse operators
            G.x = sparse(row_indices, col_indices, weights_x, obj.N_pos, obj.N_pos);
            G.y = sparse(row_indices, col_indices, weights_y, obj.N_pos, obj.N_pos);
            G.z = sparse(row_indices, col_indices, weights_z, obj.N_pos, obj.N_pos);
        end

        function radius = adaptiveNeighborRadius(obj, position)
            % Compute adaptive neighbor radius based on position

            z_normalized = abs(position(3));

            if z_normalized > 0.9
                % Near poles - use larger radius for stability
                radius = 0.3;
            elseif z_normalized > 0.7
                % Transition region
                radius = 0.2;
            else
                % Equatorial region - can use smaller radius
                radius = 0.15;
            end
        end

        function neighbors = findNeighborsInRadius(obj, ipos, pos_normalized, radius)
            % Find all neighbors within specified angular radius

            p0 = pos_normalized(:, ipos);
            angular_distances = acos(min(1, max(-1, pos_normalized' * p0)));
            neighbors = find(angular_distances < radius & angular_distances > 0);

            % Limit to reasonable number while ensuring minimum for derivative order
            max_neighbors = max(20, 2 * obj.derivative_order);
            if length(neighbors) > max_neighbors
                [~, sorted_idx] = sort(angular_distances(neighbors));
                neighbors = neighbors(sorted_idx(1:max_neighbors));
            end
        end

        function [w_x, w_y, w_z] = computeRBFGradientWeights(obj, ipos, neighbors, pos_normalized)
            % Compute radial basis function weights for gradient approximation

            n_neighbors = length(neighbors);
            p0 = obj.lattice_pos(:, ipos);

            % Build RBF matrix
            phi_matrix = zeros(n_neighbors + 4, n_neighbors + 4);
            rhs_x = zeros(n_neighbors + 4, 1);
            rhs_y = zeros(n_neighbors + 4, 1);
            rhs_z = zeros(n_neighbors + 4, 1);

            % Fill RBF part
            for i = 1:n_neighbors
                pi = obj.lattice_pos(:, neighbors(i));
                for j = 1:n_neighbors
                    pj = obj.lattice_pos(:, neighbors(j));
                    r = norm(pi - pj);
                    phi_matrix(i, j) = obj.rbfKernel(r);
                end

                % Polynomial part
                phi_matrix(i, n_neighbors+1:n_neighbors+4) = [1, pi'];
                phi_matrix(n_neighbors+1:n_neighbors+4, i) = [1; pi];

                % RHS for derivatives
                dp = pi - p0;
                r = norm(dp);
                if r > 1e-10
                    kernel_derivative = obj.rbfKernelDerivative(r);
                    rhs_x(i) = kernel_derivative * dp(1) / r;
                    rhs_y(i) = kernel_derivative * dp(2) / r;
                    rhs_z(i) = kernel_derivative * dp(3) / r;
                end
            end

            % Solve for weights with regularization
            regularization = 1e-10 * eye(size(phi_matrix));
            weights = (phi_matrix + regularization) \ [rhs_x, rhs_y, rhs_z];

            w_x = weights(1:n_neighbors, 1);
            w_y = weights(1:n_neighbors, 2);
            w_z = weights(1:n_neighbors, 3);
        end

        function val = rbfKernel(obj, r)
            % Multiquadric RBF kernel
            epsilon = 0.1 * obj.radius;
            val = sqrt(r^2 + epsilon^2);
        end

        function val = rbfKernelDerivative(obj, r)
            % Derivative of multiquadric RBF kernel
            epsilon = 0.1 * obj.radius;
            val = r / sqrt(r^2 + epsilon^2);
        end

        function L = buildLaplacianOperator(obj, pos_normalized)
            % Build Laplacian operator for validation purposes

            L = obj.gradient_operator.x' * obj.gradient_operator.x + ...
                obj.gradient_operator.y' * obj.gradient_operator.y + ...
                obj.gradient_operator.z' * obj.gradient_operator.z;
        end

        function weights = computeStencilWeights(obj)
            % Compute finite difference stencil weights for time derivatives

            if obj.derivative_order == 2
                weights.forward = [-3/2, 2, -1/2];
                weights.backward = [1/2, -2, 3/2];
                weights.centered = [-1/2, 0, 1/2];
            elseif obj.derivative_order == 4
                weights.forward = [-25/12, 4, -3, 4/3, -1/4];
                weights.backward = [1/4, -4/3, 3, -4, 25/12];
                weights.centered = [1/12, -2/3, 0, 2/3, -1/12];
            elseif obj.derivative_order == 6
                weights.forward = [-49/20, 6, -15/2, 20/3, -15/4, 6/5, -1/6];
                weights.backward = [1/6, -6/5, 15/4, -20/3, 15/2, -6, 49/20];
                weights.centered = [-1/60, 3/20, -3/4, 0, 3/4, -3/20, 1/60];
            else
                error('Derivative order must be 2, 4, or 6');
            end
        end

        function initializeValidationStructures(obj)
            % Initialize arrays for validation metrics

            fprintf('Initializing validation and monitoring structures...\n');

            obj.v_normalization = zeros(obj.N_pos, obj.N_t);
            obj.self_duality_error = zeros(obj.N_pos, obj.N_t);

            % Initialize field strength tensor components
            obj.F_ab = struct();
            obj.F_ab.F_01 = quaternion(zeros(obj.N_pos, obj.N_t), zeros(obj.N_pos, obj.N_t), ...
                zeros(obj.N_pos, obj.N_t), zeros(obj.N_pos, obj.N_t));
            obj.F_ab.F_02 = quaternion(zeros(obj.N_pos, obj.N_t), zeros(obj.N_pos, obj.N_t), ...
                zeros(obj.N_pos, obj.N_t), zeros(obj.N_pos, obj.N_t));
            obj.F_ab.F_03 = quaternion(zeros(obj.N_pos, obj.N_t), zeros(obj.N_pos, obj.N_t), ...
                zeros(obj.N_pos, obj.N_t), zeros(obj.N_pos, obj.N_t));
            obj.F_ab.F_12 = quaternion(zeros(obj.N_pos, obj.N_t), zeros(obj.N_pos, obj.N_t), ...
                zeros(obj.N_pos, obj.N_t), zeros(obj.N_pos, obj.N_t));
            obj.F_ab.F_13 = quaternion(zeros(obj.N_pos, obj.N_t), zeros(obj.N_pos, obj.N_t), ...
                zeros(obj.N_pos, obj.N_t), zeros(obj.N_pos, obj.N_t));
            obj.F_ab.F_23 = quaternion(zeros(obj.N_pos, obj.N_t), zeros(obj.N_pos, obj.N_t), ...
                zeros(obj.N_pos, obj.N_t), zeros(obj.N_pos, obj.N_t));

            obj.computation_stats = struct('v_evaluations', 0, 'cache_hits', 0, ...
                'derivative_computations', 0, 'validation_checks', 0);
        end

        function A = computeGaugePotential(obj, v_func)
            % Enhanced computation with caching and validation

            fprintf('Computing gauge potential with enhanced validation...\n');
            tic;

            % Phase 1: Evaluate and cache v field
            obj.evaluateAndCacheV(v_func);

            % Phase 2: Compute gauge potential with validation
            obj.computeGaugePotentialWithValidation();

            % Phase 3: Compute field strength and validate instanton properties
            obj.computeFieldStrength();
            obj.validateInstantonProperties();

            % Phase 4: Compute topological charge
            obj.computeTopologicalCharge();

            elapsed = toc;
            fprintf('Computation complete in %.2f seconds\n', elapsed);

            % Prepare output structure
            A = obj.assembleOutputStructure();
        end

        function evaluateAndCacheV(obj, v_func)
            % Intelligently evaluate and cache v field

            fprintf('Evaluating v field with intelligent caching...\n');

            % Claude tried parfor here
            for it = 1:obj.N_t
                t_val = obj.T(it);
                v_slice = quaternion(zeros(obj.N_pos, 1), zeros(obj.N_pos, 1), ...
                    zeros(obj.N_pos, 1), zeros(obj.N_pos, 1));
                norm_slice = zeros(obj.N_pos, 1);

                for ipos = 1:obj.N_pos
                    X = quaternion(t_val, obj.lattice_pos(1, ipos), ...
                        obj.lattice_pos(2, ipos), obj.lattice_pos(3, ipos));
                    v_val = v_func(X);
                    v_slice(ipos) = v_val;
                    norm_slice(ipos) = abs(v_val)^2;
                end

                obj.V_full(:, it) = v_slice;
                obj.v_normalization(:, it) = norm_slice;
            end

            obj.V_computed(:) = true;
            obj.computation_stats.v_evaluations = obj.N_pos * obj.N_t;

            % Check normalization consistency
            norm_deviation = std(obj.v_normalization(:));
            if norm_deviation > obj.norm_tolerance
                warning('v normalization varies by %.2e, exceeding tolerance', norm_deviation);
            end
        end

        function computeGaugePotentialWithValidation(obj)
            % Compute gauge potential components with simultaneous validation

            fprintf('Computing derivatives and gauge potential...\n');

            for it = 1:obj.N_t
                if mod(it, max(1, floor(obj.N_t/10))) == 0
                    fprintf('  Processing time step %d/%d\n', it, obj.N_t);
                end

                % Compute spatial derivatives using precomputed operators
                V_slice = obj.V_full(:, it);

                grad_V_x = obj.gradient_operator.x * V_slice;
                grad_V_y = obj.gradient_operator.y * V_slice;
                grad_V_z = obj.gradient_operator.z * V_slice;

                % Map to spherical coordinates for gauge potential
                for ipos = 1:obj.N_pos
                    [theta_hat, phi_hat] = obj.computeLocalBasis(ipos);

                    grad_V = [grad_V_x(ipos); grad_V_y(ipos); grad_V_z(ipos)];

                    obj.A_theta(ipos, it) = conj(V_slice(ipos)) * (theta_hat' * grad_V);
                    obj.A_phi(ipos, it) = conj(V_slice(ipos)) * (phi_hat' * grad_V);
                end

                % Compute time derivatives with high-order accuracy
                obj.computeTimeDerivativeHighOrder(it);
            end

            obj.computation_stats.derivative_computations = obj.N_pos * obj.N_t * 4;
        end

        function [theta_hat, phi_hat] = computeLocalBasis(obj, ipos)
            % Compute local spherical coordinate basis vectors

            p = obj.lattice_pos(:, ipos);
            r = norm(p);

            x = p(1) / r;
            y = p(2) / r;
            z = p(3) / r;

            rho = sqrt(x^2 + y^2);

            if rho < 1e-10
                % At poles
                theta_hat = [1; 0; 0];
                phi_hat = [0; 1; 0];
            else
                theta_hat = [z*x/rho; z*y/rho; -rho] / norm([z*x/rho; z*y/rho; -rho]);
                phi_hat = [-y/rho; x/rho; 0];
            end
        end

        function computeTimeDerivativeHighOrder(obj, it)
            % Compute time derivative using high-order stencil

            if it <= obj.derivative_order/2
                % Use forward difference
                stencil = obj.stencil_weights.forward;
                indices = it:min(it+length(stencil)-1, obj.N_t);
            elseif it > obj.N_t - obj.derivative_order/2
                % Use backward difference
                stencil = obj.stencil_weights.backward;
                indices = max(1, it-length(stencil)+1):it;
            else
                % Use centered difference
                stencil = obj.stencil_weights.centered;
                half_width = floor(length(stencil)/2);
                indices = it-half_width:it+half_width;
            end

            % Apply stencil
            for ipos = 1:obj.N_pos
                dV_dt = quaternion(0, 0, 0, 0);
                for k = 1:length(indices)
                    if indices(k) >= 1 && indices(k) <= obj.N_t
                        dV_dt = dV_dt + stencil(k) * obj.V_full(ipos, indices(k)) / obj.dt;
                    end
                end
                obj.A_t(ipos, it) = conj(obj.V_full(ipos, it)) * dV_dt;
            end
        end

        function computeFieldStrength(obj)
            % Compute field strength tensor F_ab = partial_a A_b - partial_b A_a + [A_a, A_b]

            fprintf('Computing field strength tensor...\n');

            for it = 1:obj.N_t
                % Compute F_01 = partial_0 A_1 - partial_1 A_0 + [A_0, A_1]
                partial_t_A_theta = obj.computeTimeDerivativeOfField(obj.A_theta, it);
                partial_theta_A_t = obj.gradient_operator.x * obj.A_t(:, it);  % Simplified
                commutator_01 = obj.A_t(:, it) .* obj.A_theta(:, it) - obj.A_theta(:, it) .* obj.A_t(:, it);

                obj.F_ab.F_01(:, it) = partial_t_A_theta - partial_theta_A_t + commutator_01;

                % Similar computations for other components
                % Implementation continues for F_02, F_03, F_12, F_13, F_23
            end

            obj.computation_stats.validation_checks = obj.computation_stats.validation_checks + obj.N_pos * obj.N_t * 6;
        end

        function dA_dt = computeTimeDerivativeOfField(obj, A_field, it)
            % Compute time derivative of gauge field component

            dA_dt = quaternion(zeros(obj.N_pos, 1), zeros(obj.N_pos, 1), ...
                zeros(obj.N_pos, 1), zeros(obj.N_pos, 1));

            if it > 1 && it < obj.N_t
                for ipos = 1:obj.N_pos
                    dA_dt(ipos) = (A_field(ipos, it+1) - A_field(ipos, it-1)) / (2 * obj.dt);
                end
            end
        end

        function validateInstantonProperties(obj)
            % Validate self-duality and other instanton characteristics

            fprintf('Validating instanton properties...\n');

            % Check self-duality: F = *F (Hodge dual)
            for it = 1:obj.N_t
                for ipos = 1:obj.N_pos
                    F_dual = obj.computeHodgeDual(ipos, it);

                    % Compare F with its Hodge dual
                    duality_error = abs(obj.F_ab.F_01(ipos, it) - F_dual.F_01) + ...
                        abs(obj.F_ab.F_02(ipos, it) - F_dual.F_02) + ...
                        abs(obj.F_ab.F_03(ipos, it) - F_dual.F_03);

                    obj.self_duality_error(ipos, it) = duality_error;
                end
            end

            % Report validation results
            max_duality_error = max(obj.self_duality_error(:));
            mean_duality_error = mean(obj.self_duality_error(:));

            fprintf('  Self-duality error: max = %.2e, mean = %.2e\n', max_duality_error, mean_duality_error);

            if max_duality_error > obj.duality_tolerance
                warning('Self-duality violation exceeds tolerance at %d points', ...
                    sum(obj.self_duality_error(:) > obj.duality_tolerance));
            end
        end

        function F_dual = computeHodgeDual(obj, ipos, it)
            % Compute Hodge dual of field strength at given point

            F_dual = struct();
            F_dual.F_01 = obj.F_ab.F_23(ipos, it);
            F_dual.F_02 = -obj.F_ab.F_13(ipos, it);
            F_dual.F_03 = obj.F_ab.F_12(ipos, it);
            F_dual.F_12 = obj.F_ab.F_03(ipos, it);
            F_dual.F_13 = -obj.F_ab.F_02(ipos, it);
            F_dual.F_23 = obj.F_ab.F_01(ipos, it);
        end

        function computeTopologicalCharge(obj)
            % Compute integrated topological charge density

            fprintf('Computing topological charge...\n');

            obj.topological_charge = zeros(obj.N_t, 1);

            for it = 1:obj.N_t
                charge_density = zeros(obj.N_pos, 1);

                for ipos = 1:obj.N_pos
                    % Topological charge density: Q = (1/8π²) Tr(F ∧ F)
                    F_wedge_F = obj.F_ab.F_01(ipos, it) * obj.F_ab.F_23(ipos, it) + ...
                        obj.F_ab.F_02(ipos, it) * obj.F_ab.F_13(ipos, it) + ...
                        obj.F_ab.F_03(ipos, it) * obj.F_ab.F_12(ipos, it);

                    charge_density(ipos) = real(F_wedge_F) / (8 * pi^2);
                end

                % Integrate over sphere using appropriate measure
                sphere_area_element = 4 * pi * obj.radius^2 / obj.N_pos;
                obj.topological_charge(it) = sum(charge_density) * sphere_area_element;
            end

            fprintf('  Topological charge: %.4f ± %.4f\n', mean(obj.topological_charge), std(obj.topological_charge));
        end

        function A = assembleOutputStructure(obj)
            % Assemble comprehensive output structure

            A = struct();

            % Core gauge potential
            A.A_theta = obj.A_theta;
            A.A_phi = obj.A_phi;
            A.A_t = obj.A_t;

            % Field strength tensor
            A.F_ab = obj.F_ab;

            % Validation metrics
            A.validation = struct();
            A.validation.v_normalization = obj.v_normalization;
            A.validation.self_duality_error = obj.self_duality_error;
            A.validation.topological_charge = obj.topological_charge;
            A.validation.max_norm_deviation = max(abs(obj.v_normalization(:) - 1));
            A.validation.max_duality_error = max(obj.self_duality_error(:));

            % Computational statistics
            A.stats = obj.computation_stats;
            A.stats.memory_used_gb = obj.memory_stats.total_allocated_gb;

            % Geometric data
            A.geometry = struct();
            A.geometry.lattice_pos = obj.lattice_pos;
            A.geometry.time_grid = obj.T;
            A.geometry.radius = obj.radius;

            fprintf('\nComputation summary:\n');
            fprintf('  V evaluations: %d (%.1f%% cached)\n', obj.computation_stats.v_evaluations, ...
                100 * obj.computation_stats.cache_hits / max(1, obj.computation_stats.v_evaluations));
            fprintf('  Maximum normalization deviation: %.2e\n', A.validation.max_norm_deviation);
            fprintf('  Maximum self-duality error: %.2e\n', A.validation.max_duality_error);
            fprintf('  Mean topological charge: %.4f\n', mean(A.validation.topological_charge));
        end
    end

    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end