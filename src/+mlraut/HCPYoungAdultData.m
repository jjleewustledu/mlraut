classdef HCPYoungAdultData < handle & mlraut.CohortData
    %% line1
    %  line2
    %  
    %  Created 08-Feb-2024 18:47:21 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    

    properties (Dependent)
        extended_task
        extended_task_dir
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
        function g = get.num_frames_to_trim(this)
            if this.is_7T
                g = 10;
                return
            end
            g = 14;
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
                g = "/Users/jjlee/mnt/CHPC_hcpdb/packages/unzip/HCP_1200";
                assert(isfolder(g));
                return
            end
            if contains(hostname, "vglab") || contains(hostname, "linux") || contains(hostname, "pascal")
                g = "/home/usr/jjlee/mnt/CHPC_hcpdb/packages/unzip/HCP_1200";
                assert(isfolder(g));
                return
            end
            if isInParallelWorker()
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
        function g = get.tr(this)
            if this.is_7T
                g = 1;
                return
            end
            g = 0.72;
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
        function this = HCPYoungAdultData(ihcp, out_dir)
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
                this mlraut.HCPYoungAdultData
                hemis {mustBeTextScalar} = "L"
            end

            if startsWith(hemis, "L", IgnoreCase=true)
                hemis = "L";
            elseif startsWith(hemis, "R", IgnoreCase=true)
                hemis = "R";
            end

            if this.is_7T
                mg = mglob(fullfile(this.mninonlinear_dir, "*."+hemis+".sphere_MSMall.164k_fs_LR.surf.gii"));
                % 995174.L.midthickness_MSMAll.164k_fs_LR.surf.gii
                if isemptytext(mg)
                    mg = mglob(fullfile(this.mninonlinear_dir, "*."+hemis+".sphere.164k_fs_LR.surf.gii"));
                    % 995174.L.midthickness.164k_fs_LR.surf.gii
                end
                assert(~isemptytext(mg), stackstr())
                g = mg(end);
                return
            end

            mg = mglob(fullfile(this.mninonlinear_dir, "fsaverage_LR32k", "*."+hemis+".sphere_MSMall.32k_fs_LR.surf.gii"));
            % 995174.L.midthickness_MSMAll.32k_fs_LR.surf.gii
            if isemptytext(mg)
                mg = mglob(fullfile(this.mninonlinear_dir, "fsaverage_LR32k", "*."+hemis+".sphere.32k_fs_LR.surf.gii"));
                % 995174.L.midthickness.32k_fs_LR.surf.gii
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
