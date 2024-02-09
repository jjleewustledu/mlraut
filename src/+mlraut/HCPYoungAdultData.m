classdef HCPYoungAdultData < handle & mlraut.CohortData
    %% line1
    %  line2
    %  
    %  Created 08-Feb-2024 18:47:21 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
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
            if contains(computer, "MAC")
                g = "/Volumes/PrecunealSSD2/AnalyticSignal";
                assert(isfolder(g));
                return
            end
            if contains(computer, "GLNXA64")
                g = fullfile(getenv("SINGULARITY_HOME"), "AnalyticSignal");
                assert(isfolder(g));
                return
            end
            if contains(hostname, "cluster")
                g = fullfile(getenv("SINGULARITY_HOME"), "AnalyticSignal");
                assert(isfolder(g));
                return
            end
            error("mlraut:NotImplementedError", stackstr());
        end
        function g = get.root_dir(~)
            if contains(computer, "MAC")
                g = "/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200";
                assert(isfolder(g));
                return
            end
            if contains(hostname, "vglab") || contains(hostname, "linux")
                g = "/home/usr/jjlee/mnt/CHPC_hcpdb/packages/unzip/HCP_1200";
                assert(isfolder(g));
                return
            end
            if contains(hostname, "cluster")
                g = "/ceph/hcpdb/packages/unzip/HCP_1200";
                assert(isfolder(g));
                return
            end
            error("mlraut:NotImplementedError", stackstr());
        end
        function g = get.task_dir(this)
            g = fullfile(this.root_dir, this.sub, "MNINonLinear", "Results", this.task);
        end
        function g = get.task_dtseries_fqfn(this)
            mg = mglob(fullfile(this.task_dir, this.task + "_Atlas_MSMAll_hp2000_clean.dtseries.nii"));
            % rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii
            if isemptytext(mg)
                mg = mglob(fullfile(this.task_dir, this.task + "_Atlas_hp2000_clean.dtseries.nii"));
                % rfMRI_REST1_RL_Atlas_hp2000_clean.dtseries.nii
            end
            assert(~isemptytext(mg), stackstr())
            g = mg(end);
        end
        function g = get.task_niigz_fqfn(this)
            mg = mglob(fullfile(this.task_dir, this.task + "_hp*_clean.nii.gz"));
            % rfMRI_REST1_RL_hp2000_clean.nii.gz
            if isemptytext(mg)
                mg = mglob(fullfile(this.task_dir, this.task + ".nii.gz"));
                % rfMRI_REST1_RL.nii.gz
            end
            assert(~isemptytext(mg), stackstr())
            g = mg(1);
        end
        function g = get.task_signal_reference_fqfn(this)
            mg = mglob(fullfile(this.task_dir, this.task + "_SBRef*.nii.gz"));
            assert(~isemptytext(mg), stackstr())
            g = mg(1);
        end
    end
    
    methods
        function this = HCPYoungAdultData(varargin)
            this = this@mlraut.CohortData(varargin{:});
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
