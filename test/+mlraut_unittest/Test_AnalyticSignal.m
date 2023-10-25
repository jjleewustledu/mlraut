classdef Test_AnalyticSignal < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 04-Dec-2022 14:02:01 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/test/+mlraut_unittest.
    %  Developed on Matlab 9.13.0.2105380 (R2022b) Update 2 for MACI64.  Copyright 2022 John J. Lee.
    
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
        function test_call(this)
            cd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200');
            this.testObj_ = mlraut.AnalyticSignal(subjects={'995174'}, do_save=true, ...
                max_frames=1192, plot_range=1:572, ...
                tasks={'rfMRI_REST1_RL'});
            call(this.testObj_)
            disp(this.testObj_)
            %mysystem('wb_view')

            % Elapsed time is 327.790908 seconds.
        end
        function test_call_do_qc(this)
            cd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200');
            this.testObj_ = mlraut.AnalyticSignal(subjects={'995174'}, do_save=true, ...
                max_frames=1192, plot_range=1:572, ...
                tasks={'rfMRI_REST1_RL'});
            call(this.testObj_, do_qc=true)
            disp(this.testObj_)

            % Elapsed time is ___ seconds.
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
        function test_call_I3CR0479(this)   
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
        function test_sig(this)
        end
        function test_num_frames_to_trim(this)
        end
        function test_max_frames(this)
        end
    end
    
    methods (TestClassSetup)
        function setupAnalyticSignal(this)
        end
    end
    
    methods (TestMethodSetup)
        function setupAnalyticSignalTest(this)
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
