classdef HCPYoungAdult7TData < handle & mlraut.CohortData
    %% line1
    %  line2
    %  
    %  Created 09-Jun-2025 00:03:40 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 24.2.0.2923080 (R2024b) Update 6 for MACA64.  Copyright 2025 John J. Lee.
    

    properties
        num_frames_to_trim = 10
        tr = 1
    end

    properties (Dependent)
        extended_task
        extended_task_dir
        json_fqfn
        out_dir
        root_dir
        stats_fqfn
        task_dtseries_fqfn
        task_niigz_fqfn
        task_ref_niigz_fqfn
        task_ref_dscalar_fqfn
        thickness_dscalar_fqfn
        t1w_fqfn
        wmparc_fqfn
    end

    methods %% GET, SET
        function g = get.extended_task(this)
            g = this.task;
        end
        function g = get.extended_task_dir(this)
            g = this.task_dir;
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
            
            g = fullfile(getenv("SINGULARITY_HOME"), "AnalyticSignalHCP");
            assert(isfolder(g));
            this.out_dir_ = g;
        end
        function     set.out_dir(this, s)
            if isemptytext(s)
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
        function g = get.stats_fqfn(this)
            mg = mglob(fullfile(this.task_dir, this.task + "_Atlas_stats.dscalar.nii"));
            % rfMRI_REST1_7T_PA_Atlas_stats.dscalar.nii, rfMRI_REST1_RL_Atlas_stats.dscalar.nii
            g = mg(end);
        end
        function g = get.task_dtseries_fqfn(this)
            mg = mglob(fullfile(this.task_dir, this.task + "_Atlas_1.6mm_MSMAll_hp2000_clean.dtseries.nii"));
            % rfMRI_REST1_7T_PA_Atlas_1.6mm_MSMAll_hp2000_clean.dtseries.nii
            if isemptytext(mg)
                mg = mglob(fullfile(this.task_dir, this.task + "_Atlas_1.6mm_hp2000_clean.dtseries.nii"));
                % rfMRI_REST1_7T_PA_Atlas_1.6mm_hp2000_clean.dtseries.nii
            end
            assert(~isemptytext(mg), stackstr())
            g = mg(end);
        end
        function g = get.task_niigz_fqfn(this)
            mg = mglob(fullfile(this.task_dir, this.task + "_hp*_clean.nii.gz"));
            % rfMRI_REST1_7T_PA_hp2000_clean.nii.gz
            if isemptytext(mg)
                mg = mglob(fullfile(this.task_dir, this.task + ".nii.gz"));
                % rfMRI_REST1_7T_PA.nii.gz
            end
            assert(~isemptytext(mg), stackstr())
            g = mg(end);
        end
        function g = get.task_ref_niigz_fqfn(this)
            mg = mglob(fullfile(this.task_dir, this.task + "_SBRef*.nii.gz"));
            assert(~isemptytext(mg), stackstr())
            g = mg(1);
        end
        function g = get.task_ref_dscalar_fqfn(this)
            mg = mglob(fullfile(this.task_dir, this.task + "_Atlas_hp2000_clean_vn.dscalar.nii"));
            if isemptytext(mg)
                mg = mglob(fullfile(this.task_dir, this.task + "_Atlas_hp2000_clean_bias.dscalar.nii"));
            end
            assert(~isemptytext(mg), stackstr())
            g = mg(1);
        end
        function g = get.thickness_dscalar_fqfn(this)
            g = fullfile(this.mninonlinear_dir, "fsaverage_LR32k", this.sub + ".thickness_MSMAll.32k_fs_LR.dscalar.nii");
            assert(isfile(g), stackstr())
        end
        function g = get.t1w_fqfn(this)
            g = fullfile(this.mninonlinear_dir, "T1w_restore.1.60.nii.gz");  % mm voxels
            assert(isfile(g), stackstr())
        end
        function g = get.wmparc_fqfn(this)
            g = fullfile(this.mninonlinear_dir, "ROIs", "wmparc.1.60.nii.gz");  % mm voxels
            assert(isfile(g), stackstr())
        end
    end
    
    methods
        function this = HCPYoungAdult7TData(ihcp, out_dir)
            arguments
                ihcp mlraut.HCP
                out_dir {mustBeTextScalar} = ""
            end

            this = this@mlraut.CohortData(ihcp);

            if isemptytext(out_dir)
                out_dir = fullfile(getenv("SINGULARITY_HOME"), "AnalyticSignalHCP");
            end
            this.out_dir_ = out_dir;
        end

        function g = surf_gii_fqfn(this, hemis)
            arguments
                this mlraut.HCPYoungAdult7TData
                hemis {mustBeTextScalar} = "L"
            end

            if startsWith(hemis, "L", IgnoreCase=true)
                hemis = "L";
            elseif startsWith(hemis, "R", IgnoreCase=true)
                hemis = "R";
            end

            mg = mglob(fullfile(this.mninonlinear_dir, "*."+hemis+".sphere_MSMall.164k_fs_LR.surf.gii"));
            % 995174.L.midthickness_MSMAll.164k_fs_LR.surf.gii
            if isemptytext(mg)
                mg = mglob(fullfile(this.mninonlinear_dir, "*."+hemis+".sphere.164k_fs_LR.surf.gii"));
                % 995174.L.midthickness.164k_fs_LR.surf.gii
            end
            assert(~isemptytext(mg), stackstr())
            g = mg(end);
        end
    end

    %% PROTECTED

    properties (Access = protected)
        out_dir_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
