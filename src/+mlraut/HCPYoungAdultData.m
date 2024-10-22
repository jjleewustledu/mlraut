classdef HCPYoungAdultData < handle & mlraut.CohortData
    %% line1
    %  line2
    %  
    %  Created 08-Feb-2024 18:47:21 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    

    properties
        num_frames_to_trim = 4
        tr = 0.72
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
    end

    methods %% GET, SET
        function g = get.atlas_fqfn(this)
            mg = mglob(fullfile(this.task_dir, this.task + "_Atlas_hp2000_clean_vn.dscalar.nii"));
            g = mg(1);
        end
        function g = get.json_fqfn(this)
            g = fullfile(this.out_dir, this.sub + ".json");  % mm voxels
            %ensuredir(myfileparts(g))
        end
        function g = get.out_dir(this)
            if ~isemptytext(this.out_dir_)
                g = this.out_dir_;
                return
            end

            if contains(computer, "MAC")
                g = "/Volumes/PrecunealSSD2/AnalyticSignalHCP";
                assert(isfolder(g));
                this.out_dir_ = g;
                return
            end
            g = fullfile(getenv("SINGULARITY_HOME"), "AnalyticSignalHCP");
            assert(isfolder(g));
            this.out_dir_ = g;
        end
        function     set.out_dir(this, s)
            if isempty(s) || isemptytext(s)
                return
            end
            assert(istext(s))
            %ensuredir(s);
            this.out_dir_ = s;
        end
        function g = get.root_dir(~)
            if contains(computer, "MAC")
                g = "/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200";
                assert(isfolder(g));
                return
            end
            if contains(hostname, "vglab") || contains(hostname, "linux") || contains(hostname, "pascal")
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
        function g = get.task_dtseries_fqfn(this)
            if this.is_7T
                mg = mglob(fullfile(this.task_dir, this.task + "_Atlas_1.6mm_MSMAll_hp2000_clean.dtseries.nii"));
                % rfMRI_REST1_7T_PA_Atlas_1.6mm_MSMAll_hp2000_clean.dtseries.nii
                if isemptytext(mg)
                    mg = mglob(fullfile(this.task_dir, this.task + "_Atlas_1.6mm_hp2000_clean.dtseries.nii"));
                    % rfMRI_REST1_7T_PA_Atlas_1.6mm_hp2000_clean.dtseries.nii
                end
                assert(~isemptytext(mg), stackstr())
                g = mg(end);
                return
            end

            mg = mglob(fullfile(this.task_dir, this.task + "_Atlas_MSMAll_hp2000_clean.dtseries.nii"));
            % rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii
            if isemptytext(mg)
                mg = mglob(fullfile(this.task_dir, this.task + "_Atlas_hp2000_clean.dtseries.nii"));
                % rfMRI_REST1_RL_Atlas_hp2000_clean.dtseries.nii
            end
            if isemptytext(mg)
                mg = mglob(fullfile(this.task_dir, this.task + "_Atlas_MSMAll.dtseries.nii"));
                % rfMRI_REST1_RL_Atlas_hp2000_clean.dtseries.nii
            end
            assert(~isemptytext(mg), stackstr())
            g = mg(end);
        end
        function g = get.task_niigz_fqfn(this)
            if this.is_7T
                mg = mglob(fullfile(this.task_dir, this.task + "_hp*_clean.nii.gz"));
                % rfMRI_REST1_7T_PA_hp2000_clean.nii.gz
                if isemptytext(mg)
                    mg = mglob(fullfile(this.task_dir, this.task + ".nii.gz"));
                    % rfMRI_REST1_7T_PA.nii.gz
                end
                assert(~isemptytext(mg), stackstr())
                g = mg(end);
                return
            end

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
        function g = get.t1w_fqfn(this)
            if this.is_7T
                g = fullfile(this.mninonlinear_dir, "T1w_restore.1.60.nii.gz");  % mm voxels
                assert(isfile(g), stackstr())
                return
            end

            g = fullfile(this.mninonlinear_dir, "T1w_restore.2.nii.gz");  % mm voxels
            assert(isfile(g), stackstr())
        end
        function g = get.wmparc_fqfn(this)
            if this.is_7T
                g = fullfile(this.mninonlinear_dir, "ROIs", "wmparc.1.60.nii.gz");  % mm voxels
                assert(isfile(g), stackstr())
                return
            end

            g = fullfile(this.mninonlinear_dir, "ROIs", "wmparc.2.nii.gz");  % mm voxels
            assert(isfile(g), stackstr())
        end
    end
    
    methods
        function this = HCPYoungAdultData(varargin)
            this = this@mlraut.CohortData(varargin{:});
        end
    end

    %% PROTECTED

    properties (Access = protected)
        out_dir_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
