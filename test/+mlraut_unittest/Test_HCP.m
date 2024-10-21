classdef Test_HCP < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 09-Feb-2024 20:43:35 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/test/+mlraut_unittest.
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
        function test_hcp(this)
            hcp = mlraut.HCP( ...
                subjects={'995174'}, tasks={'rfMRI_REST1_RL'});
            hcp.current_subject = '995174';
            hcp.current_task = 'rfMRI_REST1_RL';
            %disp(hcp)

            this.verifyEqual(mybasename(hcp.out_dir), "AnalyticSignalHCP")
            this.verifyEqual(mybasename(hcp.root_dir), "HCP_1200")
            this.verifyEqual(mybasename(hcp.task_dir), "rfMRI_REST1_RL")
            this.verifyEqual(mybasename(hcp.task_dtseries_fqfn, withext=true), "rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii")
            this.verifyEqual(mybasename(hcp.task_niigz_fqfn, withext=true), "rfMRI_REST1_RL_hp2000_clean.nii.gz")
            this.verifyEqual(mybasename(hcp.task_signal_reference_fqfn, withext=true), "rfMRI_REST1_RL_SBRef.nii.gz")
            this.verifyEqual(mybasename(hcp.t1w_fqfn, withext=true), "T1w_restore.2.nii.gz")
            this.verifyEqual(mybasename(hcp.wmparc_fqfn, withext=true), "wmparc.2.nii.gz")

            this.verifyEqual(size(hcp.task_dtseries), [1200, 91282])
            this.verifyEqual(hcp.task_niigz.filename, "rfMRI_REST1_RL_hp2000_clean.nii.gz")
            this.verifyEqual(hcp.task_signal_reference.filename, "rfMRI_REST1_RL_SBRef.nii.gz")

            %% DEPRECATED
            %this.verifyEqual(size(hcp.bold_fs_parcel), [])
            %this.verifyEqual(hcp.mask_fs_parcel, [])
        end
        function test_hcp_limited(this)

            %% empty ctor

            this.verifyTrue(isempty(this.testObj.subjects));
            this.verifyTrue(isempty(this.testObj.current_subject));
            this.verifyTrue(isempty(this.testObj.tasks));
            this.verifyTrue(isempty(this.testObj.current_task));

            this.verifyTrue(contains(this.testObj.out_dir, "AnalyticSignalHCP"));
            this.verifyTrue(contains(this.testObj.root_dir, "HCP_1200"));
            this.verifyTrue(contains(this.testObj.task_dir, fullfile("HCP_1200", "MNINonLinear", "Results")));
            this.verifyTrue(contains(this.testObj.waves_dir, fullfile("MATLAB-Drive", "arousal-waves-main")));

            return

            %% specify subject

            hcp = mlraut.HCP(subjects={'995174'});
            disp(hcp)

            %% specify subject & task

            hcp = mlraut.HCP(subjects={'995174'}, tasks={'rfMRI_REST1_RL'});
            disp(hcp)
        end
        function test_task_objects(this)
            hcp = mlraut.HCP(subjects={'995174'}, tasks={'rfMRI_REST1_RL'});

            %mysystem(sprintf("wb_view %s", hcp.task_dtseries_fqfn))
            mysystem(sprintf("fsleyes %s", hcp.task_niigz_fqfn))
            mysystem(sprintf("fsleyes %s %s %s", ...
                hcp.t1w_fqfn, ...
                hcp.wmparc_fqfn, ...
                hcp.task_signal_reference_fqfn))
        end
        function test_task_objects_7T(this)
            hcp = mlraut.HCP(subjects={'995174'}, tasks={'rfMRI_REST1_7T_PA'});

            %mysystem(sprintf("wb_view %s", hcp.task_dtseries_fqfn))
            mysystem(sprintf("fsleyes %s", hcp.task_niigz_fqfn)) % contains global signal?
            mysystem(sprintf("fsleyes %s %s %s", ...
                hcp.t1w_fqfn, ...
                hcp.wmparc_fqfn, ...
                hcp.task_signal_reference_fqfn))
        end
        function test_task_objects_Aging(this)
        end
        function test_task_objects_GBM(this)
        end
        function test_hcp_omit_late_frames(this)
        end

        %% migrate to tests of Cifti & Gifti

        function test_dlabel_nii(this)
            hcp = mlraut.HCP(subjects={'995174'}, tasks={'rfMRI_REST1_RL'});
            cii = mlraut.Cifti(hcp);
            disp(cii.aparc_a2009s_dlabel_nii())
        end
        function test_label_gii(this)
            hcp = mlraut.HCP(subjects={'995174'}, tasks={'rfMRI_REST1_RL'});
            gii = mlraut.Gifti(hcp);
            disp(gii.aparc_a2009s_label_gii())
        end
    end
    
    methods (TestClassSetup)
        function setupHCP(this)
            import mlraut.*
            this.testObj_ = HCP();
        end
    end
    
    methods (TestMethodSetup)
        function setupHCPTest(this)
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
