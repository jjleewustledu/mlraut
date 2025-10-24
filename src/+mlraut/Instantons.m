classdef Instantons < handle & mlsystem.IHandle
    %% See Ward & Wells "Twistor Geometry and Field Theory" (1990), esp. ch. 8, and AHDM (1978)
    %  
    %  Created 20-Aug-2025 16:31:40 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 24.2.0.2923080 (R2024b) Update 6 for MACA64.  Copyright 2025 John J. Lee.
    
    properties 
    end

    properties (Dependent)
        k  % instanton number (Pointrjagin index) := 1, 2, 3, ...
        lambda_sub  % \mathbb{R}^4
        x_sub  % \mathbb{H}^4
    end

    methods  %% GET/SET
        function g = get.k(this)
            g = this.k_;
        end

        function g = get.lambda_sub(this)
            g = this.lambda_sub_;
        end

        function set.lambda_sub(this, s)
            assert(isreal(s))
            assert(this.k == length(s))
            this.lambda_sub_ = s;
        end

        function g = get.x_sub(this)
            g = this.x_sub_;
        end

        function set.x_sub(this, s)
            assert(isa(s, "quaternion"))
            assert(this.k == length(s))
            this.x_sub_ = s;
        end
    end

    methods
        function A = gauge_potential(this, x)
            %% returns v^* \nabla_a v.

            arguments
                this mlraut.Instantons
                x quaternion
            end

            % The most straightforward approach employs centered differences when possible:
            % ∇_a v(x) ≈ [v(x + ε_a) - v(x - ε_a)] / (2Δx_a),
            % where ε_a represents a unit displacement in the a-th direction and Δx_a is your lattice spacing. 
            % At boundary points, you will need to use forward or backward differences as appropriate.

        end

        function G = greens_function(this, x, y)
            %% satisfies D^2 G(x, y) = -\delta(x - y)

            G = conj(this.v(x)) * this.v(y);
            G = G / (4 * pi^2 * norm(x - y)^2);
        end

        function phi_ = phi(this, x)
            %% returns 1 + \sum_{i=1}^k \lambda_i^2 \norm{x_i - x}^{-2} \in \mathbb{R}

            arguments
                this mlraut.Instantons
                x quaternion
            end

            phi_ = 1;
            for i = 1:this.k
                phi_ = phi_ + this.lambda_sub(i)^2 / norm(this.x_sub(i) - x)^2;
            end
        end

        function v_ = v(this, x, opts)
            %% returns \varphi^{-1/2} \qty[ 1 \frac{\lambda_1}{x - x_1} \cdots \frac{\lambda_k}{x - x_k} ]^* \in \mathbb{H}

            arguments
                this mlraut.Instantons
                x quaternion
                opts.return_v_ast logical = false
            end

            square_bracket = quaternion([1, 0, 0, 0]);
            for i = 1:this.k
                square_bracket = [square_bracket, this.lambda_sub(i) / (x - this.x_sub(i))]; %#ok<AGROW>
            end
            if opts.return_v_ast
                v_ = square_bracket / sqrt(this.phi(x));
            else
                v_ = conj(square_bracket) / sqrt(this.phi(x));
            end
        end

        function x_ = x(~, x_sup_a)
            %% returns \mathbb{H} := \mathbb{E}^4
            %  s.t. \sqrt{2} x = x^0 - x^1\mathbf{i} - x^2\mathbf{j} - x^3\mathbf{k}.

            x_sup_0 = x_sup_a(1);
            x_sup_1 = x_sup_a(2);
            x_sup_2 = x_sup_a(3);
            x_sup_3 = x_sup_a(4);

            x_ = quaternion(x_sup_0, -x_sup_1, -x_sup_2, -x_sup_3) / sqrt(2);
        end

        function this = Instantons(ihcp, k, opts)
            arguments
                ihcp mlraut.HCP = []
                k {mustBeInteger} = 1
                opts.lambda_sub {mustBeReal} = 1
                opts.x_sub quaternion = quaternion([1, 0, 0, 0])
            end

            this.ihcp_ = ihcp;
            this.k_ = k;
            this.lambda_sub = opts.lambda_sub;
            this.x_sub_ = opts.x_sub;
        end
    end

    %% PRIVATE

    properties (Access=private)
        ihcp_  % mlraut.HCP
        k_
        lambda_sub_
        x_sub_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
