classdef Test_AnalyticSignalHCPAging < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 12-Feb-2024 23:18:59 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/test/+mlraut_unittest.
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

        function test_task_mask_niigz(this)
            as = this.testObj;
            ic = as.task_mask_niigz;
            % as.task_ref_niigz.view_qc(ic);

            this.verifyEqual(ic.filename, "wmparc.2_binarized.nii.gz")
            this.verifyEqual(ic.qfac, -1)
            this.verifyEqual(size(ic), [91, 109, 91])
            this.verifyInstanceOf(ic.imagingFormat.img, "single")
            this.verifyEqual(dipmax(ic), 1)
            this.verifyEqual(dipsum(ic), 192280)
        end
        
        function test_task_ref_niigz(this)
            as = this.testObj;
            ic = as.task_ref_niigz;
            % ic.view_qc(as.task_mask_niigz);

            this.verifyEqual(ic.filename, "fMRI_CONCAT_ALL_SBRef.nii.gz")
            this.verifyEqual(ic.qfac, -1)
            this.verifyEqual(size(ic), [91, 109, 91])
            this.verifyInstanceOf(ic.imagingFormat.img, "single")
            this.verifyEqual(dipmax(ic), 38082.875, AbsTol=1e-3)
            this.verifyEqual(dipsum(ic), 2.305388571466583e+09, AbsTol=1)
        end

        function test_templates(this)
            as = this.testObj;

            % template_cifti ~ task_ref_dscalar_fqfn
            this.verifyTrue(isfile(as.task_ref_dscalar_fqfn));
            this.verifyTrue(isstruct(as.template_cifti.metadata));
            this.verifyTrue(iscell(as.template_cifti.diminfo));
            this.verifyEqual(as.template_cifti.diminfo{1}.models{1}.count, 29696);
            this.verifyEqual(as.template_cifti.diminfo{1}.models{1}.struct, 'CORTEX_LEFT');
            this.verifyEqual(as.template_cifti.diminfo{1}.models{1}.type, 'surf');
            this.verifyEqual(size(as.template_cifti.diminfo{1}.models{1}.vertlist), [1, 29696]);
            this.verifyEqual(as.template_cifti.diminfo{1}.models{2}.count, 29716);
            this.verifyEqual(as.template_cifti.diminfo{1}.models{2}.struct, 'CORTEX_RIGHT');
            this.verifyEqual(as.template_cifti.diminfo{1}.models{2}.type, 'surf');
            this.verifyEqual(size(as.template_cifti.diminfo{1}.models{2}.vertlist), [1, 29716]);
            this.verifyEqual(size(as.template_cifti.cdata), [91282, 1]);
            
            % template_niigz ~ wmparc
            this.verifyTrue(isfile(as.wmparc_fqfn))
            this.verifyInstanceOf(as.template_niigz, "mlfourd.ImagingContext2")
            this.verifyEqual(as.template_niigz.filename, "wmparc.2.nii.gz")
        end

        function test_memory_footprint(this)
            tic
            as = this.testObj;
            as.out_dir = "/Volumes/PrecunealSSD2/AnalyticSignalHCPAging";
            as.do_save = true;
            as.do_save_ciftis = true;
            as.do_save_subset = true;
            call(as);
            toc

            %% real memory max < 52 GB; Elapsed time is 141 seconds; 3.57 GB mat file
        end

        function test_memory_footprint_physio_suppl(this)
            %% HCP Aging has physio for individual scan sessions, not for CONCATALL,
            %  so exclude HRV, RV.

            tic
            as = this.testObj;
            as.source_physio = "iFV-brightest";
            as.source_physio_supplementary = [ ...
                "iFV-quantile", "sFV", "3rdV", "latV", "csf", "centrumsemiovale", "ctx"];
            as.out_dir = "/Volumes/PrecunealSSD2/AnalyticSignalHCPAging";
            as.do_save = true;
            as.do_save_ciftis = true;
            as.do_save_subset = true;
            call(as);
            toc

            %% real memory max < 59 GB; Elapsed time is 160 seconds; 3.57 GB mat file
        end

        function test_fultz_iFV(this)
            as = this.testObj;
            as.source_physio = "iFV";
            call_subject(as);
            
            tseries = ["bold", "-dbold/dt", "X", "Y", "Z"];
            for t = tseries
                as.plot_coherencyc(tseries=t);
            end
        end

        function test_exemplar(this)
            tic
            as = mlraut.AnalyticSignalHCPAging( ...
                subjects={'HCA9992517_V1_MR'}, ...
                tasks={'fMRI_CONCAT_ALL'}, ...
                do_global_signal_regression=false, ...
                do_save=true, ...
                do_save_ciftis=true, ...
                hp_thresh=[], ...
                lp_thresh=0.1, ...
                v_physio=50, ...
                plot_range=1:225, ...
                source_physio="iFV-brightest", ...
                tags=stackstr(use_dashes=true));
            call(as)
            toc
        end
    end
    
    methods (TestClassSetup)
        function setupAnalyticSignalHCPAging(this)
            import mlraut.*
            cd('/Volumes/PrecunealSSD2/HCPAging/HCPAgingRec/fmriresults01');
            this.testObj_ = AnalyticSignalHCPAging( ...
                subjects={'HCA9992517_V1_MR'}, ...
                tasks={'fMRI_CONCAT_ALL'}, ...
                do_global_signal_regression=true, ...
                hp_thresh=0.01, ...
                lp_thresh=0.1, ...
                v_physio=50, ...
                plot_range=1:225, ...
                source_physio="iFV", ...
                tags=stackstr(use_dashes=true));
        end
    end
    
    methods (TestMethodSetup)
        function setupAnalyticSignalHCPAgingTest(this)
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
