classdef Test_AnalyticSignalGBM < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 12-Feb-2024 23:18:30 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/test/+mlraut_unittest.
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
        function test_call_do_qc_I3CR0023(this)
            root_dir = '/Volumes/PrecunealSSD2/AnalyticSignalGBM/analytic_signal/dockerout/ciftify';
            cd(root_dir);
            mysystem('rsync -ra pascal.neuroimage.wustl.edu:~/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify/sub-I3CR0023 .')
            out_dir = '/Volumes/PrecunealSSD2/AnalyticSignalGBM/analytic_signal/matlabout/test-CE';
            ensuredir(out_dir);
            this.testObj_ = mlraut.AnalyticSignalGBM( ...
                subjects={'sub-I3CR0023'}, ...
                root_dir=root_dir, out_dir=out_dir, do_save=true, ...
                tasks={'ses-1_task-rest_run-01_desc-preproc', 'ses-1_task-rest_run-02_desc-preproc'}, ...
                source_physio='ROI', roi_fileprefix='CE_on_T1w');  
            call(this.testObj_, do_qc=true)
            disp(this.testObj_)
            cd(fullfile(out_dir, 'sub-I3CR0023'))
            saveFigures();

            % Elapsed time is ___ seconds.
        end
        function test_call_do_qc_I3CR0479(this)
            root_dir = '/Volumes/PrecunealSSD2/AnalyticSignalGBM/analytic_signal/dockerout/ciftify';
            cd(root_dir);
            mysystem('rsync -ra pascal.neuroimage.wustl.edu:~/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify/sub-I3CR0479 .')
            out_dir = '/Volumes/PrecunealSSD2/AnalyticSignalGBM/analytic_signal/matlabout/test-CE';
            ensuredir(out_dir);
            this.testObj_ = mlraut.AnalyticSignalGBM( ...
                subjects={'sub-I3CR0479'}, ...
                root_dir=root_dir, out_dir=out_dir, do_save=true, ...
                tasks={'ses-1_task-rest_run-01_desc-preproc', 'ses-1_task-rest_run-02_desc-preproc'}, ...
                source_physio='ROI', roi_fileprefix='CE_on_T1w');  
            call(this.testObj_, do_qc=true)
            disp(this.testObj_)
            cd(fullfile(out_dir, 'sub-I3CR0479'))
            saveFigures();

            % Elapsed time is ___ seconds.
        end
        function test_call_do_qc_RT126(this)
            cd('/Volumes/PrecunealSSD2/AnalyticSignalGBM/analytic_signal/dockerout_20230629111500/ciftify');
            this.testObj_ = mlraut.AnalyticSignalGBM(subjects={'sub-RT126'}, do_save=true);  
            call(this.testObj_, do_qc=true)
            disp(this.testObj_)

            % Elapsed time is ___ seconds.
        end
        function test_call_gbm_par(this)
            cd('/vgpool02/data2/jjlee/AnalyticSignalGBM/analytic_signal/dockerout/ciftify');
            out_dir = '/vgpool02/data2/jjlee/AnalyticSignalGBM/analytic_signal/matlabout';
            ensuredir(out_dir);
            asgbm = mlraut.AnalyticSignalGBMPar(subjects='sub-I3CR0479', root_dir=pwd, out_dir=out_dir, ...
                tasks={'ses-1_task-rest_run-01_desc-preproc', 'ses-1_task-rest_run-02_desc-preproc'});
            disp(asgbm)
            asgbm.call();
        end
        function test_call_I3CR0023(this)   
            %% left insula, no midline shift

            root_dir = '/Volumes/PrecunealSSD2/AnalyticSignalGBM/analytic_signal/dockerout/ciftify';
            cd(root_dir);
            mysystem('rsync -ra pascal.neuroimage.wustl.edu:~/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify/sub-I3CR0023 .')
            out_dir = '/Volumes/PrecunealSSD2/AnalyticSignalGBM/analytic_signal/matlabout/test-CE';
            ensuredir(out_dir);
            this.testObj_ = mlraut.AnalyticSignalGBM( ...
                subjects={'sub-I3CR0023'}, ...
                root_dir=root_dir, out_dir=out_dir, do_save=true, ...
                tasks={'ses-1_task-rest_run-01_desc-preproc', 'ses-1_task-rest_run-02_desc-preproc'}, ...
                source_physio='ROI', roi_fileprefix='CE_on_T1w');
            call(this.testObj_);
            disp(this.testObj_)
            mysystem('wb_view')

            % Elapsed time is ___ seconds.
        end
        function test_call_I3CR0023_physios(this)   
            %% left insula, no midline shift

            root_dir = '/Volumes/PrecunealSSD2/AnalyticSignalGBM/analytic_signal/dockerout/ciftify';
            cd(root_dir);

            out_dir = '/Volumes/PrecunealSSD2/AnalyticSignalGBM/analytic_signal/matlabout/physio_iFV';
            ensuredir(out_dir);

            this.testObj_ = mlraut.AnalyticSignalGBM( ...
                subjects={'sub-I3CR0023'}, ...
                tasks={'ses-1_task-rest_run-all_desc-preproc'}, ...
                do_save=true, ...
                do_save_ciftis=true, ...
                out_dir=out_dir, ...
                source_physio='iFV');
            call(this.testObj_);

            out_dir = '/Volumes/PrecunealSSD2/AnalyticSignalGBM/analytic_signal/matlabout/physio_CE';
            ensuredir(out_dir);

            ce = mlfourd.ImagingContext2( ...
                fullfile(root_dir, 'sub-I3CR0023', 'MNINonLinear', 'CE_on_T1w.nii.gz'));
            this.testObj_ = mlraut.AnalyticSignalGBM( ...
                subjects={'sub-I3CR0023'}, ...
                tasks={'ses-1_task-rest_run-all_desc-preproc'}, ...
                do_save=true, ...
                do_save_ciftis=true, ...
                out_dir=out_dir, ...
                roi=ce);
            call(this.testObj_);

            out_dir = '/Volumes/PrecunealSSD2/AnalyticSignalGBM/analytic_signal/matlabout/physio_WT';
            ensuredir(out_dir);

            wt = mlfourd.ImagingContext2( ...
                fullfile(root_dir, 'sub-I3CR0023', 'MNINonLinear', 'WT_on_T1w.nii.gz'));
            this.testObj_ = mlraut.AnalyticSignalGBM( ...
                subjects={'sub-I3CR0023'}, ...
                tasks={'ses-1_task-rest_run-all_desc-preproc'}, ...
                do_save=true, ...
                do_save_ciftis=true, ...
                out_dir=out_dir, ...
                roi=wt);
            call(this.testObj_);

            % Elapsed time is ___ seconds.
        end
        function test_call_I3CR0023_physio_DMN(this)  
        end
        function test_call_I3CR0479(this)   
            %% right frontal, midline shift
           
            root_dir = '/Volumes/PrecunealSSD2/AnalyticSignalGBM/analytic_signal/dockerout/ciftify';
            cd(root_dir);
            mysystem('rsync -ra pascal.neuroimage.wustl.edu:~/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify/sub-I3CR0479 .')
            out_dir = '/Volumes/PrecunealSSD2/AnalyticSignalGBM/analytic_signal/matlabout/test-CE';
            ensuredir(out_dir);
            this.testObj_ = mlraut.AnalyticSignalGBM( ...
                subjects={'sub-I3CR0479'}, ...
                root_dir=root_dir, out_dir=out_dir, do_save=true, ...
                tasks={'ses-1_task-rest_run-01_desc-preproc', 'ses-1_task-rest_run-02_desc-preproc'}, ...
                source_physio='ROI', roi_fileprefix='CE_on_T1w');
            call(this.testObj_);
            disp(this.testObj_)
            mysystem('wb_view')

            % Elapsed time is ___ seconds.
        end
        function test_call_I3CR0479_no_physio(this)   
            %% right frontal, midline shift
           
            root_dir = '~/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify';
            cd(root_dir);
            %mysystem('rsync -ra pascal.neuroimage.wustl.edu:~/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify/sub-I3CR0479 .')
            out_dir = '~/Singularity/AnalyticSignalGBM/analytic_signal/matlabout/test-CE';
            ensuredir(out_dir);
            this.testObj_ = mlraut.AnalyticSignalGBM( ...
                subjects={'sub-I3CR0479'}, ...
                root_dir=root_dir, out_dir=out_dir, do_save=true, ...
                tasks={'ses-1_task-rest_run-01_desc-preproc', 'ses-1_task-rest_run-02_desc-preproc'}, ...
                tag="no-physio", source_physio="none");
            call(this.testObj_);
            %disp(this.testObj_)
            %mysystem('wb_view')

            % Elapsed time is ___ seconds.
        end
        function test_call_RT126(this)
            cd('/Volumes/PrecunealSSD2/AnalyticSignalGBM/analytic_signal/dockerout_20230629111500/ciftify');
            this.testObj_ = mlraut.AnalyticSignalGBM(subjects={'sub-RT126'}, do_save=true);  
            call(this.testObj_)
            disp(this.testObj_)
            %mysystem('wb_view')

            % Elapsed time is ___ seconds
        end
        function test_metric_stats(this)
            metric_fun = @abs;
            metric_lbl = "abs";
            %metric_fun = @angle;
            %metric_lbl = "angle";

            out_dir = '/home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalGBM/analytic_signal/matlabout/lesionR-CE';
            %out_dir = '/home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignal';
            cd(out_dir);

            g = glob(fullfile(out_dir, 'sub-*', 'AnalyticSignalGBM_call_subject_lesionR-CE_as_proc-norm_xyzt-ROI_*.mat'));
            %g = glob(fullfile(out_dir, '*', 'AnalyticSignalPar_call_subject_as_proc-norm_xyzt-ROI_*.mat'));
            leng = length(g);
            metric_mu = zeros(1, 91282);
            metric_sig2 = zeros(1, 91282);

            for gidx = 1:leng
                try
                    ld = load(g{gidx});
                    metric_mu = metric_mu + mean(metric_fun(ld.as), 1)/leng; % mean
                catch ME
                    handwarning(ME)
                end
            end
            save(metric_lbl+"_mu.mat", "metric_mu");
            write_cifti_(metric_mu, convertStringsToChars(metric_lbl+"_mu.dscalar.nii"));
        
            for gidx = 1:leng
                try
                    ld = load(g{gidx});
                    metric_sig2 = metric_sig2 + (mean(metric_fun(ld.as),1) - metric_mu).^2/leng; % var
                catch ME
                    handwarning(ME)
                end
            end
            metric_sig = metric_sig2.^(0.5);
            save(metric_lbl+"_sig.mat", "metric_sig");
            write_cifti_(metric_sig, convertStringsToChars(metric_lbl+"_sig.dscalar.nii"));
        end
    end
    
    methods (TestClassSetup)
        function setupAnalyticSignalGBM(this)
        end
    end
    
    methods (TestMethodSetup)
        function setupAnalyticSignalGBMTest(this)
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
