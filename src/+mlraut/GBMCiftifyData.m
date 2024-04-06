classdef GBMCiftifyData < handle & mlraut.CohortData
    %% line1
    %  line2
    %  
    %  Created 08-Feb-2024 18:48:01 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    

    properties        
    end

    properties (Dependent)
        atlas_fqfn
        json_fqfn
        out_dir
        root_dir
        task_dtseries_fqfn
        task_niigz_fqfn
        task_signal_reference_fqfn   
        t1w_fqfn
        wmparc_fqfn     

        CE_fqfn
        num_frames_to_trim
        tr
        WT_fqfn
    end

    methods %% GET
        function g = get.atlas_fqfn(this)
            g = this.task_dtseries_fqfn;
        end
        function g = get.json_fqfn(this)
            g = fullfile(this.root_dir, this.sub, "gbm.json");  % mm voxels
            assert(isfile(g), stackstr())
        end
        function g = get.out_dir(this)
            if ~isemptytext(this.out_dir_)
                g = this.out_dir_;
                return
            end

            g = fullfile(getenv("SINGULARITY_HOME"), "AnalyticSignalGBM", "analytic_signal", "matlabout");
            assert(isfolder(g))
            this.out_dir_ = g;
            % /home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalGBM/analytic_signal/matlabout
        end
        function     set.out_dir(this, s)
            assert(istext(s))
            %ensuredir(s);
            this.out_dir_ = s;
        end
        function g = get.root_dir(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "AnalyticSignalGBM", "analytic_signal", "dockerout", "ciftify");
            assert(isfolder(g))
            % /home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify
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
        function g = get.t1w_fqfn(this)
            g = fullfile(this.mninonlinear_dir, "T1w.nii.gz");  % mm voxels
            assert(isfile(g), stackstr())
        end
        function g = get.wmparc_fqfn(this)
            g = fullfile(this.mninonlinear_dir, "wmparc.nii.gz");  % mm voxels
            assert(isfile(g), stackstr())
        end

        function g = get.CE_fqfn(this)
            g = fullfile(this.mninonlinear_dir, "CE_on_T1w.nii.gz");  % mm voxels
            assert(isfile(g), stackstr())
        end
        function g = get.num_frames_to_trim(this)
            g = 0;
            % g = this.json.num_frames_to_trim;
        end
        function g = get.tr(this)
            g = this.json.tr;
        end
        function g = get.WT_fqfn(this)
            g = fullfile(this.mninonlinear_dir, "WT_on_T1w.nii.gz");  % mm voxels
            assert(isfile(g), stackstr())
        end
    end

    methods
        function this = GBMCiftifyData(varargin)            
            this = this@mlraut.CohortData(varargin{:});
        end
    end

    %% PROTECTED

    properties (Access = protected)
        out_dir_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
