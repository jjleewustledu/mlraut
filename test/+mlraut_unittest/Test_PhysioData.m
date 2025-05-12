classdef Test_PhysioData < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 09-Feb-2024 01:09:10 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/test/+mlraut_unittest.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    
    properties
        ashcp
        testObj
    end
    
    methods (Test)
        function test_afun(this)
            import mlraut.*
            this.assumeEqual(1,1);
            this.verifyEqual(1,1);
            this.assertEqual(1,1);
        end

        function test_task_niigz(this)
            %% tests that hcp.task_niigz is consistent with repeated calls

            this.verifyEqual(this.ashcp.num_frames_ori, 1200)
            this.verifyEqual(this.ashcp.num_frames_to_trim, 4)
            this.verifyEqual(this.ashcp.num_frames, 1196)
            this.verifyEqual(size(this.ashcp.task_niigz()), [91, 109, 91, 1196])
            this.verifyEqual(size(this.ashcp.task_niigz()), [91, 109, 91, 1196])
            this.verifyEqual(size(this.ashcp.task_niigz()), [91, 109, 91, 1196])
        end

        function test_iFV(this)
        end

        function test_gray(this)
        end

        function test_white(this)
        end

        function test_csf(this)
        end

        function test_thalamus(this)
            indices = [10, 49];
            physroi = mlraut.PhysioRoi( ...
                this.ashcp, from_wmparc_indices=indices);
            %physroi.view_qc();
            this.verifyEqual(crc32_adler(0, double(physroi.roi_mask)), uint32(1999844326))
            bold = call(physroi);  % toc ~ 19 s
            %plot(bold)
            this.verifyEqual(crc32_adler(0, bold), uint32(1375506269))
        end

        function test_iFV_vs_RV(this)
        end

        function test_iFV_vs_HRV(this)
        end

        function test_flipLR(this)
        end

        function test_centered(this)
            t = 1;
            hcp_ = this.ashcp;
            hcp_.current_task = hcp_.tasks{t};

            physio_0 = hcp_.task_physio();
            figure; 
            plot(physio_0);
            hold on

            physio_1 = ...
                hcp_.build_centered( ...
                hcp_.build_global_signal_regressed(physio_0));
            plot(physio_1);
            hold off

            legend(["task_physio()", "build_centered()"], ...
                Interpreter="none");
        end

        function test_rescaled(this)
            t = 1;
            hcp_ = this.ashcp;
            hcp_.current_task = hcp_.tasks{t};

            physio_0 = hcp_.task_physio();

            physio_1 = ...
                hcp_.build_centered( ...
                hcp_.build_global_signal_regressed(physio_0));
            figure; 
            plot(physio_1);
            hold on

            physio_2 = ...
                hcp_.build_rescaled(physio_1);
            plot(physio_2);
            hold off

            legend(["build_centered()", "build_rescaled()"], ...
                Interpreter="none");
        end

        function test_band_passed(this)
            t = 1;
            hcp_ = this.ashcp;
            hcp_.current_task = hcp_.tasks{t};

            physio_0 = hcp_.task_physio();

            physio_1 = ...
                hcp_.build_centered_and_rescaled( ...
                hcp_.build_global_signal_regressed(physio_0));
            figure; 
            plot(physio_1);
            hold on

            physio_2 = ...
                hcp_.build_band_passed(physio_1);
            plot(physio_2);
            hold off

            legend(["build_centered_and_rescaled()", "build_band_passed()"], ...
                Interpreter="none");
        end

        function test_analytic_signal(this)
            %% see also mlraut.AnalyticSignalHCP.call_subject_late_hilbert()

            t = 1;
            hcp_ = this.ashcp;
            hcp_.current_task = hcp_.tasks{t};

            % BOLD
            try
                bold_ = ...
                    hcp_.build_band_passed( ...
                    hcp_.build_centered_and_rescaled( ...
                    hcp_.build_global_signal_regressed(hcp_.task_dtseries())));
                bold_ = bold_(:, 91000);
            catch ME
                disp([hcp_.current_subject ' ' hcp_.current_task ' BOLD missing or defective:']);
                handerror(ME)
            end
            hcp_.plot3(z=hilbert(bold_));

            % physio
            try
                physio_ = ...
                    hcp_.build_band_passed( ...
                    hcp_.build_centered_and_rescaled( ...
                    hcp_.build_global_signal_regressed(hcp_.task_physio())));
            catch ME
                disp([hcp_.current_subject ' ' hcp_.current_task ' physio missing or defective:']);
                handerror(ME)
            end
            hcp_.plot3(z=hilbert(physio_));

            as_ = conj(hilbert(physio_)).*hilbert(bold_);
            hcp_.plot3(z=as_);

        end
    end
    
    methods (TestClassSetup)
        function setupPhysioData(this)
        end
    end
    
    methods (TestMethodSetup)
        function setupPhysioDataTest(this)
            cd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200');
            this.ashcp = mlraut.AnalyticSignalHCP( ...
                subjects={'995174'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                do_7T=false, ...
                do_resting=true, ...
                do_task=false, ...
                do_global_signal_regression=true, ...
                do_save=false, ...
                do_save_dynamic=false, ...
                do_save_ciftis=false, ...
                force_band=false, ...
                hp_thresh=0.01, ...
                lp_thresh=0.1, ...
                source_physio="iFV-brightest", ...
                tags=stackstr(use_dashes=true));
            this.ashcp.malloc();
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
