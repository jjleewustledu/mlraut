classdef NetworkCollections < handle & mlraut.NetworkData
    %% Provides weights ctx, cbm, str, thal to efficiently calculate total anatomical BOLD signals as vec
    %  from HCP_signals.  See also mlraut.FultzMulti.
    %  
    %  Created 06-May-2025 14:07:01 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 25.1.0.2897550 (R2025a) Prerelease Update 5 for MACA64.  Copyright 2025 John J. Lee.
    
    properties
        anat_weights
    end

    properties (Dependent)
        anat_mask
        anat_weights_fqfn
    end

    methods %% GET
        function g = get.anat_mask(~)
            error("mlraut:NotImplementedError", stackstr())
        end
        function g = get.anat_weights_fqfn(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "AnalyticSignalHCP", "NetworkCollections_anat_weights.mat");
        end
    end

    methods
        function this = NetworkCollections(varargin)
            %% 
            %  ihcp mlraut.AnalyticSignal {mustBeNonempty}
            %  psi {mustBeNumeric} : supply []

            this = this@mlraut.NetworkData(varargin{:});

            if isfile(this.anat_weights_fqfn)
                ld = load(this.anat_weights_fqfn);
                this.anat_weights = ld.anat_weights;
            end

            if isempty(this.anat_weights)
                ctx_data = mlraut.CorticalData(varargin{:});
                anat_weights.ctx = this.build_weights(ctx_data);

                cbm_data = mlraut.CerebellarData(varargin{:});
                anat_weights.cbm = this.build_weights(cbm_data);

                str_data = mlraut.StriatalData(varargin{:});
                anat_weights.str = this.build_weights(str_data);

                thal_data = mlraut.ThalamicData(varargin{:});
                anat_weights.thal = this.build_weights(thal_data);

                this.anat_weights = anat_weights;
                save(this.anat_weights_fqfn, "anat_weights");
            end
        end
    end

    %% PROTECTED

    methods (Access = protected)
        function w = build_weights(this, data)
            anat = logical(data.anat_mask);
            anat_mass = sum(anat, "all");
            Nyeo = sum(~contains(this.NETWORKS_YEO_NAMES, 'task'));
            w = zeros(1, numel(this.NETWORKS_YEO_NAMES));
            for y = 1:Nyeo
                w(y) = sum(this.networks_Yeo == y & anat, "all") / anat_mass;
            end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
