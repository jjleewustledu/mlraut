classdef HCPAgingData < handle & mlraut.CohortData
    %% line1
    %  line2
    %  
    %  Created 08-Feb-2024 18:47:29 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    

    properties
        num_frames_to_trim = 4
        tr = 0.8
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

    methods %% GET
        function g = get.extended_task(this)
            if strcmp(this.task, "rfMRI_REST")
                g = "fMRI_CONCAT_ALL";
                return
            end
            g = this.task;
        end
        function g = get.extended_task_dir(this)
            ext_root_dir = strrep(this.root_dir, "HCPAgingRec", "rfMRIExtended");
            g = fullfile(ext_root_dir, this.sub, "MNINonLinear", "Results", this.extended_task);
        end
        function g = get.json_fqfn(this)
            g = fullfile(this.out_dir, this.sub, this.sub + ".json");  % mm voxels
            ensuredir(myfileparts(g))
        end
        function g = get.out_dir(this)
            if ~isemptytext(this.out_dir_)
                g = this.out_dir_;
                return
            end

            if contains(computer, "MAC")
                g = "/Volumes/PrecunealSSD2/AnalyticSignalHCPAging";
                assert(isfolder(g));
                this.out_dir_ = g;
                return
            end
            if contains(computer, "GLNXA64")
                g = fullfile(getenv("SINGULARITY_HOME"), "AnalyticSignalHCPAging");
                assert(isfolder(g));
                this.out_dir_ = g;
                return
            end
            if contains(hostname, "cluster")
                g = fullfile(getenv("SINGULARITY_HOME"), "AnalyticSignalHCPAging");
                assert(isfolder(g));
                this.out_dir_ = g;
                return
            end
            error("mlraut:NotImplementedError", stackstr());
        end
        function     set.out_dir(this, s)
            assert(istext(s))
            %ensuredir(s);
            this.out_dir_ = s;
        end
        function g = get.root_dir(~)
            if contains(computer, "MAC")
                g = "/Volumes/PrecunealSSD2/HCPAging/HCPAgingRec/fmriresults01";
                assert(isfolder(g));
                return
            end
            if contains(computer, "GLNXA64")
                g = fullfile(getenv("SINGULARITY_HOME"), "HCPAging", "HCPAgingRec", "fmriresults01");
                assert(isfolder(g));
                return
            end
            if contains(hostname, "cluster")
                g = fullfile(getenv("SINGULARITY_HOME"), "HCPAging", "HCPAgingRec", "fmriresults01");
                assert(isfolder(g));
                return
            end
            error("mlraut:NotImplementedError", stackstr());
        end
        function g = get.stats_fqfn(this)
            mg = mglob(fullfile(this.extended_task_dir, this.extended_task + "_Atlas_MSMAll_mean.dscalar.nii"));
            % rfMRI_REST1_AP_Atlas_MSMAll_mean.dscalar.nii
            if isemptytext(mg)
                mg = mglob(fullfile(this.extended_task_dir, this.extended_task + "_Atlas_mean.dscalar.nii"));
                % rfMRI_REST1_AP_Atlas_mean.dscalar.nii
            end
            assert(~isemptytext(mg), stackstr())
            g = mg(1);
        end
        function g = get.task_dtseries_fqfn(this)
            mg = mglob(fullfile(this.extended_task_dir, this.extended_task + "_Atlas_MSMAll_hp0_clean.dtseries.nii"));
            % fMRI_CONCAT_ALL_Atlas_MSMAll_hp0_clean.dtseries.nii
            if isemptytext(mg)
                mg = mglob(fullfile(this.extended_task_dir, this.extended_task + "_Atlas_hp0_clean.dtseries.nii"));
                % fMRI_CONCAT_ALL_Atlas_hp0_clean.dtseries.nii
            end
            if isemptytext(mg)
                mg = mglob(fullfile(this.task_dir, this.task + "_Atlas_MSMAll_hp0_clean.dtseries.nii"));
                % rfMRI_REST_Atlas_MSMAll_hp0_clean.dtseries.nii
                % rfMRI_REST1_AP_Atlas_MSMAll_hp0_clean.dtseries.nii
            end
            if isemptytext(mg)
                mg = mglob(fullfile(this.task_dir, this.task + "_Atlas_hp0_clean.dtseries.nii"));
                % rfMRI_REST_Atlas_MSMAll_hp0_clean.dtseries.nii
                % rfMRI_REST1_AP_Atlas_MSMAll_hp0_clean.dtseries.nii
            end
            assert(~isemptytext(mg), stackstr())
            g = mg(1);
        end
        function g = get.task_niigz_fqfn(this)
            mg = mglob(fullfile(this.extended_task_dir, this.extended_task + "_hp*_clean.nii.gz"));
            assert(~isemptytext(mg), stackstr())
            g = mg(1);
        end
        function g = get.task_ref_niigz_fqfn(this)
            mg = mglob(fullfile(this.extended_task_dir, this.extended_task + "_SBRef.nii.gz"));
            if isemptytext(mg)                
                mg = mglob(fullfile(this.extended_task_dir, this.extended_task + "_mean.nii.gz"));
            end
            assert(~isemptytext(mg), stackstr())
            g = mg(1);
        end
        function g = get.task_ref_dscalar_fqfn(this)
            mg = mglob(fullfile(this.extended_task_dir, this.extended_task + "_Atlas_MSMAll_hp0_clean_vn.dscalar.nii"));
            if isemptytext(mg)
                mg = mglob(fullfile(this.extended_task_dir, this.extended_task + "_Atlas_MSMAll_hp0_vn.dscalar.nii"));
            end
            if isemptytext(mg)
                mg = mglob(fullfile(this.extended_task_dir, this.extended_task + "_Atlas_MSMAll_mean.dscalar.nii"));
            end
            if isemptytext(mg)
                mg = mglob(fullfile(this.extended_task_dir, this.extended_task + "_Atlas_hp0_clean_vn.dscalar.nii"));
            end
            if isemptytext(mg)
                mg = mglob(fullfile(this.extended_task_dir, this.extended_task + "_Atlas_hp0_vn.dscalar.nii"));
            end
            if isemptytext(mg)
                mg = mglob(fullfile(this.extended_task_dir, this.extended_task + "_Atlas_mean.dscalar.nii"));
            end
            assert(~isemptytext(mg), stackstr())
            g = mg(1);
        end
        function g = get.thickness_dscalar_fqfn(this)
            g = fullfile(this.mninonlinear_dir, "fsaverage_LR32k", this.sub + ".thickness_MSMAll.32k_fs_LR.dscalar.nii");
            assert(isfile(g), stackstr())
        end
        function g = get.t1w_fqfn(this)
            g = fullfile(this.mninonlinear_dir, "T1w_restore.2.nii.gz");  % mm voxels
            assert(isfile(g), stackstr())
        end
        function g = get.wmparc_fqfn(this)
            g = fullfile(this.mninonlinear_dir, "ROIs", "wmparc.2.nii.gz");  % mm voxels
            assert(isfile(g), stackstr())
        end
    end
    
    methods
        function this = HCPAgingData(varargin)
            this = this@mlraut.CohortData(varargin{:});
        end

        function g = surf_gii_fqfn(this, hemis)
            arguments
                this mlraut.HCPAgingData
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
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
