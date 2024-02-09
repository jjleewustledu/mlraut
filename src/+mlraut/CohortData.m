classdef CohortData < handle & mlsystem.IHandle
    %% line1
    %  line2
    %  
    %  Created 08-Feb-2024 18:43:56 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    

    properties (Abstract)
        out_dir
        root_dir
        task_dir
        task_dtseries_fqfn
        task_niigz_fqfn
        task_signal_reference_fqfn  
    end

    properties (Dependent)
        sub
        task
    end

    methods %% GET
        function g = get.sub(this)
            g = this.ihcp_.current_sub;
        end
        function g = get.task(this)
            g = this.ihcp_.current_task;
        end
    end

    methods (Static)
        function this = create(ihcp)
            arguments
                ihcp mlraut.HCP
            end

            switch class(ihcp)
                case 'mlraut.AnalyticSignal'
                    this = mlraut.HCPYoungAdultData(ihcp);
                case 'mlraut.AnalyticSignalYoungAdult'
                    this = mlraut.HCPYoungAdultData(ihcp);
                case 'mlraut.AnalyticSignalGBM'
                    this = mlraut.GBMCiftifyData(ihcp);
                case 'mlraut.AnalyticSignalHCPAging'
                    this = mlraut.HCPAgingData(ihcp);
                otherwise
                    error("mlraut:ValueError", "%s: received an %s object.", stackstr(), class(ihcp))
            end
        end
    end

    %% PROTECTED

    properties (Access = protected)
        ihcp_
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
