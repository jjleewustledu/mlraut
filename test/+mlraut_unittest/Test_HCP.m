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

        function test_ctor(this)
            hcp = this.testObj;  % HCP

            this.verifyEqual(hcp.max_frames, Inf);
            this.verifyEqual(hcp.current_subject, '995174')
            this.verifyEqual(hcp.current_task, 'rfMRI_REST1_RL')
            this.verifyEqual(hcp.Fs, 1.389, AbsTol=1e-3);
            this.verifyEqual(hcp.num_frames, 1196);
            this.verifyEqual(hcp.num_frames_ori, 1200);
            this.verifyEqual(hcp.num_frames_to_trim, 4);
            this.verifyEqual(hcp.num_nodes, 91282);
            this.verifyTrue(contains(hcp.out_dir, "AnalyticSignalHCP"));
            this.verifyTrue(contains(hcp.root_dir, "HCP_1200"));
            this.verifyTrue(contains(hcp.task_dir, fullfile("MNINonLinear", "Results", "rfMRI_REST1_RL")));
            this.verifyEqual(hcp.tr, 0.72);
            this.verifyEqual(mybasename(hcp.task_dtseries_fqfn, withext=true), "rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii")
            this.verifyEqual(mybasename(hcp.task_niigz_fqfn, withext=true), "rfMRI_REST1_RL_hp2000_clean.nii.gz")
            this.verifyEqual(mybasename(hcp.task_ref_niigz_fqfn, withext=true), "rfMRI_REST1_RL_SBRef.nii.gz")
            this.verifyEqual(mybasename(hcp.task_ref_dscalar_fqfn, withext=true), "rfMRI_REST1_RL_Atlas_hp2000_clean_vn.dscalar.nii")
            this.verifyEqual(mybasename(hcp.t1w_fqfn, withext=true), "T1w_restore.2.nii.gz")
            this.verifyTrue(contains(hcp.waves_dir, fullfile("MATLAB-Drive", "arousal-waves")));
            this.verifyEqual(mybasename(hcp.wmparc_fqfn, withext=true), "wmparc.2.nii.gz")

            %% template_cifti ~ Atlas_hp2000_clean_vn.dscalar, Cifti.average_times(task.dseries)
         
            cii = hcp.template_cifti;
            this.verifyTrue(isstruct(cii.diminfo{1}));
            this.verifyEqual(cii.diminfo{1}.type, 'dense');
            this.verifyTrue(isstruct(cii.diminfo{1}.vol));
            this.verifyTrue(iscell(cii.diminfo{1}.models));
            this.verifyEqual(cii.diminfo{1}.length, 91282);
            this.verifyEqual(cii.diminfo{2}.type, 'scalars');
            this.verifyEqual(cii.diminfo{2}.length, 1);
            this.verifyEqual(size(cii.cdata), [91282, 1]);
            this.verifyEqual(dipmax(cii.cdata), 772.802673339844, AbsTol=1e-3);

            %% template_niigz ~ wmparc

            this.verifyTrue(isfile(hcp.wmparc_fqfn))
            this.verifyInstanceOf(hcp.template_niigz, "mlfourd.ImagingContext2")
            this.verifyEqual(hcp.template_niigz.filename, "wmparc.2.nii.gz")
        end

        function test_malloc(this)
            hcp = this.testObj;

            this.verifyInstanceOf(hcp.bold_data, "mlraut.BOLDData");
            this.verifyEqual(hcp.bold_data.num_frames_ori, 1200);
            this.verifyEqual(hcp.bold_data.num_nodes, 91282);

            this.verifyInstanceOf(hcp.cifti, "mlraut.Cifti");
            this.verifyTrue(endsWith(hcp.out_dir, "AnalyticSignalHCP"));
            this.verifyEqual(size(hcp.template_cifti.cdata), [91282, 1]);
            this.verifyEqual(dipmax(hcp.template_cifti.cdata), 772.802673339844, AbsTol=1);

            this.verifyInstanceOf(hcp.cohort_data, "mlraut.HCPYoungAdultData");
            this.verifyEqual(hcp.cohort_data.tr, 0.72);
            this.verifyEqual(hcp.cohort_data.sub, '995174')
            this.verifyEqual(hcp.cohort_data.task, 'rfMRI_REST1_RL')

            this.verifyInstanceOf(hcp.twistors, "mlraut.Twistors");
        end

        function test_max_frames(this)
            hcp = this.testObj;
            hcp.max_frames = 320;
            this.verifyEqual(hcp.num_frames, 320);
        end

        function test_task_dtseries(this)
            hcp = this.testObj;
            dtseries = hcp.task_dtseries();
            this.verifyEqual(size(dtseries), [1196, 91282])
        end

        function test_task_niigz(this)
            hcp = this.testObj;
            ic = hcp.task_niigz();            
            this.verifyEqual(ic.filename, "rfMRI_REST1_RL_hp2000_clean.nii.gz")
            this.verifyEqual(ic.qfac, -1)
            this.verifyEqual(size(ic), [91, 109, 91, 1196])
            this.verifyInstanceOf(ic.imagingFormat.img, "single")
        end

        function test_task_mask_niigz(this)
            hcp = this.testObj;
            ic = hcp.task_mask_niigz();
            % hcp.task_mask_niigz.view_qc(ic)

            this.verifyEqual(ic.filename, "wmparc.2_binarized.nii.gz")
            this.verifyEqual(ic.qfac, -1)
            this.verifyEqual(size(ic), [91, 109, 91])
            this.verifyInstanceOf(ic.imagingFormat.img, "single")
            this.verifyEqual(dipmax(ic), 1)
            this.verifyEqual(dipsum(ic), 195405)
        end
        
        function test_task_ref_niigz(this)
            hcp = this.testObj;
            ic = hcp.task_ref_niigz;
            % ic.view_qc(hcp.task_mask_niigz);

            this.verifyEqual(ic.filename, "rfMRI_REST1_RL_SBRef.nii.gz")
            this.verifyEqual(ic.qfac, -1)
            this.verifyEqual(size(ic), [91, 109, 91])
            this.verifyInstanceOf(ic.imagingFormat.img, "single")
            this.verifyEqual(dipmax(ic), 40297.08203125, AbsTol=1e-3)
            this.verifyEqual(dipsum(ic), 2.272379999944073e+09, AbsTol=1)
        end

        function test_task_objects(this)
            return

            hcp = mlraut.HCP(subjects={'995174'}, tasks={'rfMRI_REST1_RL'});
            mysystem(sprintf("fsleyes %s", hcp.task_niigz_fqfn))
            mysystem(sprintf("fsleyes %s %s %s", ...
                hcp.t1w_fqfn, ...
                hcp.wmparc_fqfn, ...
                hcp.task_ref_niigz_fqfn))
        end
        function test_task_objects_7T(this)
            return

            hcp = mlraut.HCP(subjects={'995174'}, tasks={'rfMRI_REST1_7T_PA'});
            malloc(hcp);
            mysystem(sprintf("fsleyes %s", hcp.task_niigz_fqfn))
            mysystem(sprintf("fsleyes %s %s %s", ...
                hcp.t1w_fqfn, ...
                hcp.wmparc_fqfn, ...
                hcp.task_ref_niigz_fqfn))
        end

        function test_dlabel_nii(this)
            hcp = this.testObj;
            cii = mlraut.Cifti(hcp);
            dlabel = cii.aparc_a2009s_dlabel_nii();
            this.verifyEqual(size(dlabel.metadata), [1, 4])
            this.verifyEqual(size(dlabel.diminfo), [1, 2])
            this.verifyEqual(size(dlabel.cdata), [59412, 1])
        end

        function test_label_gii(this)
            hcp = this.testObj;
            gii = mlraut.Gifti(hcp);
            label = gii.aparc_a2009s_label_gii();
            this.verifyEqual(size(label.cdata), [32492, 1])
            this.verifyEqual(size(label.labels), [1, 1])
        end
    end
    
    methods (TestClassSetup)
        function setupHCP(this)
            import mlraut.*
            this.testObj_ = HCP(subjects="995174", tasks="rfMRI_REST1_RL");
        end
    end
    
    methods (TestMethodSetup)
        function setupHCPTest(this)
            this.testObj = copy(this.testObj_);
            malloc(this.testObj);
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
