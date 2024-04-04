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
        function test_ctor(this)
            cd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200');
            as = mlraut.AnalyticSignal(subjects={'995174'}, tasks={'rfMRI_REST1_RL'});
            this.verifyEqual(as.num_nets, 9);
            this.verifyFalse(isemptytext(as.current_subject));
            this.verifyFalse(isemptytext(as.current_task));
            
            cifti_last = as.cifti_last;
            [a,b,c] = cifti_last.metadata.value;
            this.verifyTrue(contains(a, "Commit"));
            this.verifyTrue(contains(b, "wb_command -cifti-convert"));
            this.verifyTrue(strcmp(c, '/HCP/hcpdb/build_ssd/chpc/BUILD/HCP_Staging/DeDriftAndResample_1449032106_995174/995174/MNINonLinear/Results/rfMRI_REST1_RL/rfMRI_REST1_RL_hp2000.ica'));
            this.verifyTrue(isstruct(cifti_last.diminfo{1}));
            this.verifyEqual(cifti_last.diminfo{1}.type, 'dense');
            this.verifyTrue(isstruct(cifti_last.diminfo{1}.vol));
            this.verifyTrue(iscell(cifti_last.diminfo{1}.models));
            this.verifyEqual(cifti_last.diminfo{1}.length, 91282);
            this.verifyEqual(cifti_last.diminfo{2}.type, 'series');
            this.verifyEqual(cifti_last.diminfo{2}.length, 1200);
            this.verifyEqual(cifti_last.diminfo{2}.seriesStart, 0);
            this.verifyEqual(cifti_last.diminfo{2}.seriesStep, 0.7200);
            this.verifyEqual(cifti_last.diminfo{2}.seriesUnit, 'SECOND');
            this.verifyEqual(size(cifti_last.cdata), [91282, 1200]);

            this.verifyEqual(as.Fs, 1.38888888888889, RelTol=10*eps);
            this.verifyEqual(as.num_frames, 1192);
            this.verifyEqual(as.num_frames_ori, 1200);
            this.verifyEqual(as.num_frames_to_trim, 4);
            this.verifyEqual(as.num_nodes, 91282);
            this.verifyTrue(contains(as.out_dir, "AnalyticSignal"));
            this.verifyTrue(contains(as.root_dir, "HCP_1200"));
            this.verifyTrue(contains(as.task_dir, fullfile("Results", "rfMRI_REST1_RL")));
            this.verifyTrue(contains(as.task_dtseries_fqfn, "Atlas_MSMAll_hp2000_clean"));
            this.verifyTrue(contains(as.task_niigz_fqfn, "hp2000_clean"));
            this.verifyTrue(contains(as.task_signal_reference_fqfn, "SBRef"));
            this.verifyTrue(contains(as.t1w_fqfn, "T1w_restore.2.nii.gz"));
            this.verifyEqual(as.tr, 0.72);
            this.verifyTrue(contains(as.waves_dir, fullfile("MATLAB-Drive", "arousal-waves")));
            this.verifyTrue(contains(as.wmparc_fqfn, fullfile("ROIs", "wmparc.2.nii.gz")));
            this.verifyTrue(contains(as.workbench_dir, "workbench"));
            disp(as)
        end
        function test_ctor_7T(this)
            cd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200');
            as = mlraut.AnalyticSignal(subjects={'995174'}, tasks={'rfMRI_REST_7T_PA'});
            disp(as)
        end
        function test_analytic_bold(this)
            %% following hilbert transform, ensure:

            % - sensible amplitudes and angles

            % - sensible residual from global signal

            % - sensible normalization (centering & rescaling)

        end
        function test_iFV(this)
        end
        function test_physioROI(this)
        end
        function test_physioRV(this)
        end
        function test_physioHRV(this)
        end
        function test_bra_ket(this)
        end

        function test_call(this)
            cd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200');
            as = mlraut.AnalyticSignal( ...
                subjects={'995174'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                do_save=true, ...
                do_save_ciftis=true, ...
                tags=stackstr(use_dashes=true));
            call(as);
            disp(as)

            % Elapsed time is 327.790908 seconds.
        end
        function test_call_dphysio(this)
            cd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200');
            for p = ["iFV" "RV", "HRV"]
                as = mlraut.AnalyticSignal( ...
                    subjects={'995174'}, ...
                    tasks={'rfMRI_REST1_RL'}, ...
                    do_save=true, ...
                    do_save_ciftis=true, ...
                    force_band=false, ...
                    tags=stackstr(use_dashes=true), ...
                    source_physio=p(1));
                call(as);
                disp(as)
            end

            % Elapsed time is 327.790908 seconds.
        end
        function test_call_roi_from_wmparc_indices(this)
            cd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200');
            as = mlraut.AnalyticSignal( ...
                subjects={'995174'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                do_save=true, ...
                do_save_ciftis=true, ...
                force_band=false, ...
                tags=stackstr(use_dashes=true)+"-deepwhite", ...
                source_physio="ROI", ...
                roi=[5001, 5002]);
            call(as);
            disp(as)
            as = mlraut.AnalyticSignal( ...
                subjects={'995174'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                do_save=true, ...
                do_save_ciftis=true, ...
                force_band=false, ...
                tags=stackstr(use_dashes=true)+"-csf", ...
                source_physio="ROI", ...
                roi=[1, 4, 5, 24, 43, 44]);
            call(as);
            disp(as)

            % Elapsed time is 327.790908 seconds.
        end
        function test_call_no_physio(this)
            cd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200');
            as = mlraut.AnalyticSignal( ...
                subjects={'995174'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                do_save=true, ...
                do_save_ciftis=true, ...
                force_band=false, ...
                tags=stackstr(use_dashes=true), ...
                source_physio="none");
            call(as)
            disp(as)

            % Elapsed time is 327.790908 seconds.
        end
        function test_call_all_tasks(this)
            cd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200');
            as = mlraut.AnalyticSignal( ...
                subjects={'995174'}, ...
                tasks={}, ...
                do_save=true, ...
                do_save_ciftis=true, ...
                tags=stackstr(use_dashes=true), ...
                source_physio="iFV");
            call(as);
            disp(as)

            % Elapsed time is 327.790908 seconds.
        end
        function test_call_7T(this)
            cd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200');
            as = mlraut.AnalyticSignal( ...
                subjects={'995174'}, ...
                tasks={'rfMRI_REST1_7T_PA'}, ...
                do_save=true, ...
                do_save_ciftis=true, ...
                force_band=false, ...
                tags=stackstr(use_dashes=true), ...
                source_physio="iFV");
            call(as);
            disp(as)

            % Elapsed time is 327.790908 seconds.
        end
        function test_call_qc(this)
            cd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200');
            as = mlraut.AnalyticSignal( ...
                subjects={'995174'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                do_save=true, ...
                source_physio="iFV");
            call(as, do_qc=true)
            disp(as)

            % Elapsed time is ___ seconds.
        end
        function test_num_frames_to_trim(this)
        end
        function test_max_frames(this)
        end
        function test_average_network_signals(this)
            %% candidate supplemental figure, with concat runs

            cd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200');
            as = mlraut.AnalyticSignal( ...
                subjects={'995174'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                do_save=false, ...
                force_band=false, ...
                tags=stackstr(use_dashes=true));

            bold = as.task_dtseries();
            bold = as.build_global_signal_regressed(bold);

            yeos = as.average_network_signals(bold);
            yeo_names = mlraut.NetworkData.NETWORKS_YEO_NAMES;
            for anat = ["cbm", "ctx", "str", "thal"]
                figure
                hold on
                for y = 1:length(yeo_names)
                    plot(yeos.(anat(1))(:, y)); 
                end
                title(sprintf("yeos.%s", anat(1)));  
                legend(yeo_names)
                hold off

                figure;     
                tiledlayout(3,3);
                for y = 1:length(yeo_names)
                    nexttile
                    histogram(yeos.(anat(1))(:, y)); title(sprintf("yeos(%s)", yeo_names{y}), 100);
                    title(sprintf("yeos.%s(%s)", anat(1), yeo_names{y}));    
                end
            end
        end
        function test_task_dtseries(this)
            cd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200');
            as = mlraut.AnalyticSignal( ...
                subjects={'995174'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                do_save=false, ...
                force_band=false, ...
                tags=stackstr(use_dashes=true));

            bold = as.task_dtseries();
            gs = as.build_global_signal_for(bold);

            figure; plot(gs); title("gs");
            figure; plot(bold(:, 600)); title("bold(:, 600)");
            figure; plot(bold(:, 600) - gs); title("bold(:, 600) - gs");
            figure; histogram(bold); title("histogram(bold)");
            figure; histogram(gs); title("histogram(gs)");
            figure; histogram(bold - gs); title("histogram(bold - gs)")
        end
        function test_physio_HRV(this)
            cd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200');
            as = mlraut.AnalyticSignal( ...
                subjects={'995174'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                do_save=false, ...
                force_band=false, ...
                tags=stackstr(use_dashes=true), ...
                source_physio='HRV');

            bold = as.task_niigz();
            HRV = mlraut.PhysioHRV(as, bold);
            physio = HRV.call();
            gs = as.build_global_signal_for(physio);

            figure; plot(physio); title("physio");
            figure; plot(gs); title("gs");
            figure; histogram(physio); title("histogram(physio)");
            figure; histogram(gs); title("histogram(gs)");
        end
        function test_physio_iFV(this)
            cd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200');
            as = mlraut.AnalyticSignal( ...
                subjects={'995174'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                do_save=false, ...
                force_band=false, ...
                tags=stackstr(use_dashes=true), ...
                source_physio='iFV');

            bold = as.task_niigz();
            iFV = mlraut.IFourthVentricle(as, bold);
            iFV.wmparc.view_qc(iFV.ifv_mask);
            physio = iFV.call();   
            gs = as.build_global_signal_for(physio);

            figure; plot(physio); title("physio");
            figure; plot(gs); title("gs");
            figure; plot(physio - gs); title("physio - gs");
            figure; histogram(physio); title("histogram(physio)");
            figure; histogram(gs); title("histogram(gs)");
            figure; histogram(physio - gs); title("histogram(physio - gs)");
        end
        function test_physio_ROI(this)
            cd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200');
            as = mlraut.AnalyticSignal( ...
                subjects={'995174'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                do_save=false, ...
                force_band=false, ...
                tags=stackstr(use_dashes=true), ...close a
                source_physio='ROI');

            bold = as.task_niigz();
            Roi = mlraut.PhysioRoi(as, bold, from_wmparc_indices=[1 4 5 24 43 44]);  % csf
            Roi.wmparc.view();
            physio = Roi.call();
            gs = as.build_global_signal_for(physio);

            figure; plot(physio); title("physio");
            figure; plot(gs); title("gs");
            figure; plot(physio - gs); title("physio - gs");
            figure; histogram(physio); title("histogram(physio)");
            figure; histogram(gs); title("histogram(gs)");
            figure; histogram(physio - gs); title("histogram(physio - gs)");
        end
        function test_physio_RV(this)
            cd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200');
            as = mlraut.AnalyticSignal( ...
                subjects={'995174'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                do_save=false, ...
                force_band=false, ...
                tags=stackstr(use_dashes=true), ...
                source_physio='RV');

            bold = as.task_niigz();
            RV = mlraut.PhysioRV(as, bold);
            physio = RV.call();
            gs = as.build_global_signal_for(physio);

            figure; plot(physio); title("physio");
            figure; plot(gs); title("gs");
            figure; histogram(physio); title("histogram(physio)");
            figure; histogram(gs); title("histogram(gs)");
        end
        function test_task_signal_mask(this)
            cd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200');
            as = mlraut.AnalyticSignal( ...
                subjects={'995174'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                do_save=false, ...
                force_band=false, ...
                tags=stackstr(use_dashes=true), ...
                source_physio='RV');
            as.task_signal_reference.view( ...
                as.task_signal_mask);
        end
        function test_histograms(this)
            %figure; histogram(this.bold_signal_); title("bold_signal_", interpreter="none")
            %figure; histogram(this.physio_signal_); title("physio_signal_", interpreter="none")
            %figure; histogram(abs(this.bold_signal_)); title("abs(bold_signal_ in C)", interpreter="none")
            %figure; histogram(abs(this.physio_signal_)); title("abs(physio_signal_ in C)", interpreter="none")
            %figure; histogram(abs(this.analytic_signal_)); title("abs(analytic_signal_ in C)", interpreter="none")    
        end

        function test_ASPar_call(this)
            cd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200');
            as = mlraut.AnalyticSignalPar( ...
                subjects={'995174'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                do_save=true, ...
                do_save_ciftis=false, ...
                tags=stackstr(use_dashes=true));
            call(as);
            disp(as)

            % Elapsed time is 327.790908 seconds.
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
