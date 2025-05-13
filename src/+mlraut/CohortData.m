classdef CohortData < handle & mlsystem.IHandle
    %% line1
    %  line2
    %  
    %  Created 08-Feb-2024 18:43:56 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    

    properties (Abstract)
        json_fqfn
        num_frames_to_trim
        out_dir
        root_dir
        stats_fqfn
        task_dtseries_fqfn
        task_niigz_fqfn
        task_ref_niigz_fqfn
        task_ref_dscalar_fqfn
        thickness_dscalar_fqfn
        tr
        t1w_fqfn        
        wmparc_fqfn
    end

    properties (Dependent)
        is_7T
        json
        mninonlinear_dir
        sub
        task
        task_dir
    end

    methods %% GET, SET
        function g = get.is_7T(this)
            g = contains(this.task, '7T');
        end
        function g = get.json(this)
            if ~isempty(this.json_)
                g = this.json_;
                return
            end

            this.json_ = jsonread(this.json_fqfn);
            g = this.json_;
        end
        function     set.json(this, s)
            assert(isstruct(s), stackstr())
            this.json_ = s;
        end
        function g = get.mninonlinear_dir(this)
            g = fullfile(this.root_dir, this.sub, "MNINonLinear");
        end
        function g = get.sub(this)
            g = this.ihcp_.current_subject;
        end
        function g = get.task(this)
            g = this.ihcp_.current_task;
        end
        function g = get.task_dir(this)
            g = fullfile(this.root_dir, this.sub, "MNINonLinear", "Results", this.task);
        end
    end

    methods
        function g = surf_gii_fqfn(~, varargin)
            error("mlraut:NotImplementedError", stackstr());
        end
    end

    methods (Static)
        function this = create(ihcp)
            arguments
                ihcp mlraut.HCP
            end

            switch class(ihcp)
                case 'mlraut.HCP'
                    this = mlraut.HCPYoungAdultData(ihcp);
                case {'mlraut.AnalyticSignalHCP', 'mlraut.AnalyticSignalHCPPar'}
                    this = mlraut.HCPYoungAdultData(ihcp);
                case {'mlraut.AnalyticSignalGBM', 'mlraut.AnalyticSignalGBMPar'}
                    this = mlraut.GBMCiftifyData2(ihcp);
                case {'mlraut.AnalyticSignalHCPAging', 'mlraut.AnalyticSignalHCPAgingPar'}
                    this = mlraut.HCPAgingData(ihcp);
                otherwise
                    error("mlraut:ValueError", "%s: received an %s object.", stackstr(), class(ihcp))
            end
        end
    end

    %% PROTECTED

    properties (Access = protected)
        ihcp_
        json_
    end

    methods (Access = protected)
        function this = CohortData(ihcp)
            arguments
                ihcp mlraut.HCP
            end

            this.ihcp_ = ihcp;
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
