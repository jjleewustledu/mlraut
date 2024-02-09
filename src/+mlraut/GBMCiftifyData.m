classdef GBMCiftifyData < handle & mlraut.CohortData
    %% line1
    %  line2
    %  
    %  Created 08-Feb-2024 18:48:01 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    

    properties (Dependent)
        out_dir
        root_dir
        task_dir
        task_dtseries_fqfn
        task_niigz_fqfn
        task_signal_reference_fqfn        
    end

    methods %% GET
        function g = get.out_dir(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "AnalyticSignalGBM", "analytic_signal", "matlabout");
            % /home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalGBM/analytic_signal/matlabout
        end
        function g = get.root_dir(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "AnalyticSignalGBM", "analytic_signal", "dockerout", "ciftify");
            % /home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify
        end
        function g = get.task_dir(this)
            g = fullfile(this.root_dir, this.sub, "MNINonLinear", "Results", this.task);
        end
        function g = get.task_dtseries_fqfn(this)
            mg = mglob(fullfile(this.task_dir, sprintf("%s*_Atlas_s0.dtseries.nii", this.task)));
            assert(~isemptytext(mg), stackstr())
            g = mg(1);
        end
        function g = get.task_niigz_fqfn(this)
            mg = mglob(fullfile(this.task_dir, sprintf("%s*.nii.gz", this.task)));
            assert(~isemptytext(mg), stackstr())
            g = mg(1);
        end
        function g = get.task_signal_reference_fqfn(this)
            mg = mglob(fullfile(this.task_dir, sprintf("%s*.nii.gz", this.task)));
            assert(~isemptytext(mg), stackstr())
            g = mg(1);
        end
    end

    methods
        function this = GBMCiftifyData(varargin)            
            this = this@mlraut.CohortData(varargin{:});
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
