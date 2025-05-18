classdef Test_CohortData < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 08-Feb-2024 18:43:56 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/test/+mlraut_unittest.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    
    properties
        testObj
    end
    
    methods (Test)
        function test_afun(this)
            import mlraut.*
            this.assumeEqual(1,1);
            this.verifyEqual(1,1);
            this.assertEqual(1,1);
        end
        function test_HCPYoungAdultData(this)
            ihcp = mlraut.AnalyticSignalHCP( ...
                subjects={'996782'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                hp_thresh=[], ...
                lp_thresh=0.1, ...
                source_physio="iFV-brightest", ...
                tags=stackstr(use_dashes=true));
            malloc(ihcp);
            c = mlraut.CohortData.create(ihcp);

            this.verifyEqual(        c.num_frames_to_trim, 4);
            this.verifyEqual(                        c.tr, 0.7200);
            this.verifyEqual(                 c.json_fqfn, "/Users/jjlee/Singularity/AnalyticSignalHCP/996782.json");
            this.verifyEqual(                   c.out_dir, "/Users/jjlee/Singularity/AnalyticSignalHCP");
            this.verifyEqual(                  c.root_dir, "/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200");
            this.verifyEqual(                c.stats_fqfn, "/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200/996782/MNINonLinear/Results/rfMRI_REST1_RL/rfMRI_REST1_RL_Atlas_stats.dscalar.nii");
            this.verifyEqual(        c.task_dtseries_fqfn, "/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200/996782/MNINonLinear/Results/rfMRI_REST1_RL/rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii");
            this.verifyEqual(           c.task_niigz_fqfn, "/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200/996782/MNINonLinear/Results/rfMRI_REST1_RL/rfMRI_REST1_RL_hp2000_clean.nii.gz");
            this.verifyEqual(       c.task_ref_niigz_fqfn, "/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200/996782/MNINonLinear/Results/rfMRI_REST1_RL/rfMRI_REST1_RL_SBRef.nii.gz");
            this.verifyEqual(     c.task_ref_dscalar_fqfn, "/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200/996782/MNINonLinear/Results/rfMRI_REST1_RL/rfMRI_REST1_RL_Atlas_hp2000_clean_vn.dscalar.nii");
            this.verifyEqual(    c.thickness_dscalar_fqfn, "/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200/996782/MNINonLinear/fsaverage_LR32k/996782.thickness_MSMAll.32k_fs_LR.dscalar.nii");
            this.verifyEqual(                  c.t1w_fqfn, "/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200/996782/MNINonLinear/T1w_restore.2.nii.gz");
            this.verifyEqual(               c.wmparc_fqfn, "/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200/996782/MNINonLinear/ROIs/wmparc.2.nii.gz");
            this.verifyEqual(                     c.is_7T, false);
            this.verifyEqual(          c.mninonlinear_dir, "/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200/996782/MNINonLinear");
            this.verifyEqual(                       c.sub, '996782');
            this.verifyEqual(                      c.task, 'rfMRI_REST1_RL');
            this.verifyEqual(                  c.task_dir, "/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200/996782/MNINonLinear/Results/rfMRI_REST1_RL");
        end
        function test_HCPAgingData(this)
            ihcp = mlraut.AnalyticSignalHCPAging( ...
                subjects={'HCA9992517_V1_MR'}, ...
                tasks={'fMRI_CONCAT_ALL'}, ...
                hp_thresh=[], ...
                lp_thresh=0.1, ...
                source_physio="iFV-brightest", ...
                tags=stackstr(use_dashes=true));
            malloc(ihcp);
            c = mlraut.CohortData.create(ihcp);

            this.verifyEqual(        c.num_frames_to_trim, 4);
            this.verifyEqual(                        c.tr, 0.8000);
            this.verifyEqual(             c.extended_task, 'fMRI_CONCAT_ALL');
            this.verifyEqual(         c.extended_task_dir, "/Volumes/PrecunealSSD2/HCPAging/rfMRIExtended/fmriresults01/HCA9992517_V1_MR/MNINonLinear/Results/fMRI_CONCAT_ALL");
            this.verifyEqual(                 c.json_fqfn, "/Users/jjlee/Singularity/AnalyticSignalHCPAging/HCA9992517_V1_MR/HCA9992517_V1_MR.json");
            this.verifyEqual(                   c.out_dir, "/Users/jjlee/Singularity/AnalyticSignalHCPAging");
            this.verifyEqual(                  c.root_dir, "/Volumes/PrecunealSSD2/HCPAging/HCPAgingRec/fmriresults01");
            this.verifyEqual(                c.stats_fqfn, "/Volumes/PrecunealSSD2/HCPAging/rfMRIExtended/fmriresults01/HCA9992517_V1_MR/MNINonLinear/Results/fMRI_CONCAT_ALL/fMRI_CONCAT_ALL_Atlas_MSMAll_mean.dscalar.nii");
            this.verifyEqual(        c.task_dtseries_fqfn, "/Volumes/PrecunealSSD2/HCPAging/rfMRIExtended/fmriresults01/HCA9992517_V1_MR/MNINonLinear/Results/fMRI_CONCAT_ALL/fMRI_CONCAT_ALL_Atlas_MSMAll_hp0_clean.dtseries.nii");
            this.verifyEqual(           c.task_niigz_fqfn, "/Volumes/PrecunealSSD2/HCPAging/rfMRIExtended/fmriresults01/HCA9992517_V1_MR/MNINonLinear/Results/fMRI_CONCAT_ALL/fMRI_CONCAT_ALL_hp0_clean.nii.gz");
            this.verifyEqual(       c.task_ref_niigz_fqfn, "/Volumes/PrecunealSSD2/HCPAging/rfMRIExtended/fmriresults01/HCA9992517_V1_MR/MNINonLinear/Results/fMRI_CONCAT_ALL/fMRI_CONCAT_ALL_SBRef.nii.gz");
            this.verifyEqual(     c.task_ref_dscalar_fqfn, "/Volumes/PrecunealSSD2/HCPAging/rfMRIExtended/fmriresults01/HCA9992517_V1_MR/MNINonLinear/Results/fMRI_CONCAT_ALL/fMRI_CONCAT_ALL_Atlas_MSMAll_hp0_clean_vn.dscalar.nii");
            this.verifyEqual(    c.thickness_dscalar_fqfn, "/Volumes/PrecunealSSD2/HCPAging/HCPAgingRec/fmriresults01/HCA9992517_V1_MR/MNINonLinear/fsaverage_LR32k/HCA9992517_V1_MR.thickness_MSMAll.32k_fs_LR.dscalar.nii");
            this.verifyEqual(                  c.t1w_fqfn, "/Volumes/PrecunealSSD2/HCPAging/HCPAgingRec/fmriresults01/HCA9992517_V1_MR/MNINonLinear/T1w_restore.2.nii.gz");
            this.verifyEqual(               c.wmparc_fqfn, "/Volumes/PrecunealSSD2/HCPAging/HCPAgingRec/fmriresults01/HCA9992517_V1_MR/MNINonLinear/ROIs/wmparc.2.nii.gz");
            this.verifyEqual(                     c.is_7T, false);
            this.verifyEqual(          c.mninonlinear_dir, "/Volumes/PrecunealSSD2/HCPAging/HCPAgingRec/fmriresults01/HCA9992517_V1_MR/MNINonLinear");
            this.verifyEqual(                       c.sub, 'HCA9992517_V1_MR');
            this.verifyEqual(                      c.task, 'fMRI_CONCAT_ALL');
            this.verifyEqual(                  c.task_dir, "/Volumes/PrecunealSSD2/HCPAging/HCPAgingRec/fmriresults01/HCA9992517_V1_MR/MNINonLinear/Results/fMRI_CONCAT_ALL");
        end
        function test_GBMCiftifyData2(this)
            ihcp = mlraut.AnalyticSignalGBM( ...
                subjects={'sub-I3CR1488'}, ...
                tasks={'ses-1_task-rest_run-01_desc-preproc', 'ses-1_task-rest_run-02_desc-preproc'}, ...
                do_resting=true, ...
                source_physio="iFV-brightest");
            malloc(ihcp);
            c = mlraut.CohortData.create(ihcp);

            this.verifyEqual(               c.num_frames_to_trim, 0);
            this.verifyEqual(                          c.CE_fqfn, "/Users/jjlee/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify/sub-I3CR1488/MNINonLinear/CE_on_T1w.nii.gz");
            this.verifyEqual(                      size(c.CE_ic), [91, 109, 91]);
            this.verifyEqual(    c.ciftify_subject_fmri_log_fqfn, "/Users/jjlee/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify/sub-I3CR1488/MNINonLinear/Results/ses-1_task-rest_run-01_desc-preproc/ciftify_subject_fmri.log");
            this.verifyEqual(                    c.datashare_dir, "/Users/jjlee/Singularity/AnalyticSignalGBM/GBM_datashare");
            this.verifyEqual(                       c.edema_fqfn, "/Users/jjlee/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify/sub-I3CR1488/MNINonLinear/edema_on_T1w.nii.gz");
            this.verifyEqual(                   size(c.edema_ic), [91, 109, 91]);
            this.verifyEqual(                        c.json_fqfn, "/Users/jjlee/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify/sub-I3CR1488/gbm.json");
            this.verifyEqual(                          c.out_dir, "/Users/jjlee/Singularity/AnalyticSignalGBM/analytic_signal/matlabout");
            this.verifyEqual(                 c.rsFC_PreProc_loc, "jjlee@linux1.neuroimage.wustl.edu:/data/nil-bluearc/shimony/bidhan/rsFC_PreProc");
            this.verifyEqual(                         c.root_dir, "/Users/jjlee/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify");
            this.verifyEqual(                       c.stats_fqfn, "/Users/jjlee/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify/sub-I3CR1488/MNINonLinear/Results/ses-1_task-rest_run-01_desc-preproc/ses-1_task-rest_run-01_desc-preproc_Atlas_s0.dtseries.nii");
            this.verifyEqual(               c.task_dtseries_fqfn, "/Users/jjlee/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify/sub-I3CR1488/MNINonLinear/Results/ses-1_task-rest_run-01_desc-preproc/ses-1_task-rest_run-01_desc-preproc_Atlas_s0.dtseries.nii");
            this.verifyEqual(                  c.task_niigz_fqfn, "/Users/jjlee/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify/sub-I3CR1488/MNINonLinear/Results/ses-1_task-rest_run-01_desc-preproc/ses-1_task-rest_run-01_desc-preproc.nii.gz");
            this.verifyEqual(              c.task_ref_niigz_fqfn, "/Users/jjlee/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify/sub-I3CR1488/MNINonLinear/Results/ses-1_task-rest_run-01_desc-preproc/ses-1_task-rest_run-01_desc-preproc.nii.gz");
            this.verifyEqual(            c.task_ref_dscalar_fqfn, "/Users/jjlee/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify/sub-I3CR1488/MNINonLinear/Results/ses-1_task-rest_run-01_desc-preproc/ses-1_task-rest_run-01_desc-preproc_Atlas_s0_avgt.dscalar.nii");
            this.verifyEqual(           c.thickness_dscalar_fqfn, "/Users/jjlee/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify/sub-I3CR1488/MNINonLinear/fsaverage_LR32k/sub-I3CR1488.thickness.32k_fs_LR.dscalar.nii");
            this.verifyEqual(                         c.t1w_fqfn, "/Users/jjlee/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify/sub-I3CR1488/MNINonLinear/T1w.nii.gz");
            this.verifyEqual(                               c.tr, 2.6100);
            this.verifyEqual(                      c.wmparc_fqfn, "/Users/jjlee/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify/sub-I3CR1488/MNINonLinear/wmparc.nii.gz");
            this.verifyEqual(                          c.WT_fqfn, "/Users/jjlee/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify/sub-I3CR1488/MNINonLinear/WT_on_T1w.nii.gz");
            this.verifyEqual(                      size(c.WT_ic), [91, 109, 91]);
            this.verifyEqual(                size(c.map_rt_i3cr), [281, 1]);
            this.verifyEqual(             size(c.table_excluded), [15, 117]);
            this.verifyEqual(                  size(c.table_gbm), [288, 116]);
            this.verifyEqual(              size(c.table_rt_i3cr), [142, 3]);
            this.verifyEqual(                            c.is_7T, false);
            this.verifyInstanceOf(                        c.json, "struct");
            this.verifyEqual(                 c.mninonlinear_dir, "/Users/jjlee/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify/sub-I3CR1488/MNINonLinear");
            this.verifyEqual(                              c.sub, 'sub-I3CR1488');
            this.verifyEqual(                             c.task, 'ses-1_task-rest_run-01_desc-preproc');
            this.verifyEqual(                         c.task_dir, "/Users/jjlee/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify/sub-I3CR1488/MNINonLinear/Results/ses-1_task-rest_run-01_desc-preproc");
        end
    end
    
    methods (TestClassSetup)
        function setupCohortData(this)
            import mlraut.*
            this.testObj_ = [];
        end
    end
    
    methods (TestMethodSetup)
        function setupCohortDataTest(this)
            this.testObj = this.testObj_;
            this.addTeardown(@this.cleanTestMethod)
        end
    end
    
    properties (Access = private)
        testObj_
    end
    
    methods (Access = private)
        function cleanTestMethod(this)
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
