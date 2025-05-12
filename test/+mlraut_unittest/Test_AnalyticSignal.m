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
        function test_physio(this)

            this.verifyEqual(this.testObj.current_subject, this.testObj.subjects{1});

            size_bold = [1196, 91282];
            T = 1195 * this.testObj.tr;
            times = 0:this.testObj.tr:T;
            physios = "latV";
            %physios = [ ...
            %    "iFV-brightest", "iFV-quantile", "iFV", "sFV", "4thV", "3rdV", "latV"];
            for phys = physios
                this.testObj.source_physio = phys;
                [~,physio_vec,pROI] = this.testObj.task_physio(size_reference=size_bold);
                figure; plot(times, asrow(real(physio_vec)))
                pROI.view_qc();  % (this.testObj.t1w_fqfn);
            end
        end
        function test_build_band_passed(this)
            obj = this.testObj;
            size_bold = [1196, 91282];
            obj.filter_order = 8;
            obj.force_legacy_butter = false;

            try
                obj.malloc();
                obj.current_task = obj.tasks{1};
                [~,physio__] = obj.task_physio(size_reference=size_bold);
                figure; plot(physio__); 
            catch ME
                handexcept(ME)
            end
        end
        function test_ctor(this)
            as = mlraut.AnalyticSignalHCP(subjects={'995174'}, tasks={'rfMRI_REST1_RL'});
            this.verifyEqual(as.num_nets, 9);
            this.verifyFalse(isemptytext(as.current_subject));
            this.verifyFalse(isemptytext(as.current_task));
            
            template_cifti = as.template_cifti;
            [a,b,c] = template_cifti.metadata.value;
            this.verifyTrue(contains(a, "Commit"));
            this.verifyTrue(contains(b, "wb_command -cifti-convert"));
            this.verifyTrue(strcmp(c, '/HCP/hcpdb/build_ssd/chpc/BUILD/HCP_Staging/DeDriftAndResample_1449032106_995174/995174/MNINonLinear/Results/rfMRI_REST1_RL/rfMRI_REST1_RL_hp2000.ica'));
            this.verifyTrue(isstruct(template_cifti.diminfo{1}));
            this.verifyEqual(template_cifti.diminfo{1}.type, 'dense');
            this.verifyTrue(isstruct(template_cifti.diminfo{1}.vol));
            this.verifyTrue(iscell(template_cifti.diminfo{1}.models));
            this.verifyEqual(template_cifti.diminfo{1}.length, 91282);
            this.verifyEqual(template_cifti.diminfo{2}.type, 'series');
            this.verifyEqual(template_cifti.diminfo{2}.length, 1200);
            this.verifyEqual(template_cifti.diminfo{2}.seriesStart, 0);
            this.verifyEqual(template_cifti.diminfo{2}.seriesStep, 0.7200);
            this.verifyEqual(template_cifti.diminfo{2}.seriesUnit, 'SECOND');
            this.verifyEqual(size(template_cifti.cdata), [91282, 1200]);

            this.verifyEqual(as.Fs, 1.38888888888889, RelTol=10*eps);
            this.verifyEqual(as.num_frames, 1192);
            this.verifyEqual(as.num_frames_ori, 1200);
            this.verifyEqual(as.num_frames_to_trim, 4);
            this.verifyEqual(as.num_nodes, 91282);
            this.verifyTrue(contains(as.out_dir, "AnalyticSignalHCP"));
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
            as = mlraut.AnalyticSignalHCP( ...
                subjects={'995174'}, ...
                do_7T=true, ...
                tasks={'rfMRI_REST_7T_PA'}, ...
                do_save=true, ...
                do_save_ciftis=true);
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
            as = this.testObj;
            as.do_plot_emd = true;
            as.do_plot_networks = true;
            as.do_save = true;
            as.do_save_ciftis = true;
            as.do_save_ciftis_of_diffs = true;
            as.do_save_dynamic = true;

            call(as);
            disp(as)

            % Elapsed time is 327.790908 seconds.
        end
        function test_call_dphysio(this)
            for p = ["iFV-brightest" "RV", "HRV"]
                as = mlraut.AnalyticSignalHCP( ...
                    subjects={'995174'}, ...
                    tasks={'rfMRI_REST1_RL'}, ...
                    do_save=false, ...
                    do_save_ciftis=false, ...
                    force_band=false, ...
                    tags=stackstr(use_dashes=true), ...
                    source_physio=p(1));
                call(as);
                disp(as)
            end

            % Elapsed time is 327.790908 seconds.
        end
        function test_call_roi_from_wmparc_indices(this)
            as = mlraut.AnalyticSignalHCP( ...
                subjects={'995174'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                do_save=false, ...
                do_save_ciftis=false, ...
                force_band=false, ...
                tags=stackstr(use_dashes=true)+"-deepwhite", ...
                source_physio="ROI", ...
                roi=[5001, 5002]);
            call(as);
            disp(as)
            as = mlraut.AnalyticSignalHCP( ...
                subjects={'995174'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                do_save=false, ...
                do_save_ciftis=false, ...
                force_band=false, ...
                tags=stackstr(use_dashes=true)+"-csf", ...
                source_physio="ROI", ...
                roi=[1, 4, 5, 24, 43, 44]);
            call(as);
            disp(as)

            % Elapsed time is 327.790908 seconds.
        end
        function test_call_no_physio(this)
            as = mlraut.AnalyticSignalHCP( ...
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
        function test_call_task(this)
            as = mlraut.AnalyticSignalHCP( ...
                subjects={'995174'}, ...
                tasks={'tfMRI_LANGUAGE_LR'}, ...
                do_save=true, ...
                do_save_ciftis=true, ...
                do_save_dynamic=true, ...
                force_band=false, ...
                hp_thresh=0.005, ...
                lp_thresh=0.1, ...
                tags=stackstr(use_dashes=true), ...
                source_physio="iFV-brightest");
            call(as);
            disp(as)

            % Elapsed time is 327.790908 seconds.
        end
        function test_call_7T(this)
            as = mlraut.AnalyticSignalHCP( ...
                do_7T=true, ...
                subjects={'995174'}, ...
                tasks={'rfMRI_REST2_7T_AP'}, ...
                do_save=true, ...
                do_save_ciftis=true, ...
                do_save_dynamic=true, ...
                force_band=false, ...
                hp_thresh=0.005, ...
                lp_thresh=0.1, ...
                tags=stackstr(use_dashes=true), ...
                source_physio="iFV-brightest");
            call(as);
            disp(as)

            % Elapsed time is 327.790908 seconds.
        end

        function test_call_all(this)
            for do7t = [false, true]
                as = mlraut.AnalyticSignalHCP( ...
                    do_7T=do7t, ...
                    do_resting=true, ...
                    do_task=true, ...
                    subjects={'995174'}, ...
                    tasks={}, ...
                    do_save=true, ...
                    do_save_ciftis=true, ...
                    do_save_dynamic=true, ...
                    force_band=false, ...
                    hp_thresh=0.005, ...
                    lp_thresh=0.1, ...
                    tags=stackstr(use_dashes=true), ...
                    source_physio="iFV-brightest");
                call(as);
                disp(as.tasks)
            end

            % Elapsed time is 327.790908 seconds.
        end

        function test_call_qc(this)
            as = mlraut.AnalyticSignalHCP( ...
                subjects={'995174'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                do_save=true, ...
                source_physio="iFV-brightest");
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

            as = mlraut.AnalyticSignalHCP( ...
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

            % parameterize testing
            NETWORKS_YEO_NAMES = ...
                {'visual', 'somatomotor', 'dorsal attention', 'ventral attention', 'limbic', ...
                'frontoparietal', 'default mode', ...
                'task+', 'task-'};
            select_rsn = 'default mode';
            select_network_type = "cortical";
            selection = sprintf("%s, %s", select_network_type, select_rsn);

            as = this.testObj;
            as.tags_user=stackstr(use_dashes=true);
            as.source_physio = "none";

            bold_net_type = as.task_dtseries(network_type=select_network_type);
            bold_rsn = bold_net_type(:, contains(NETWORKS_YEO_NAMES, select_rsn));
            gs = as.build_global_signal_for(bold_rsn);

            figure; plot(gs); title("gs");
            figure; plot(bold_rsn); title("bold(:, " + selection + ")");
            figure; plot(bold_rsn - gs); title("bold(:, " + selection + ") - gs");

            figure; histogram(gs); title("histogram(gs)");
            figure; histogram(bold_rsn); title("histogram(bold(:, " + selection + "))");
            figure; histogram(bold_rsn - gs); title("histogram(bold(:, " + selection + ") - gs)");

            as.fit_power_law(x=gs,title="test_task_dtseries:  gs");
            as.fit_power_law(x=bold_rsn,title="test_task_dtseries:  bold(:, " + selection + ")");
            as.fit_power_law(x=(bold_rsn - gs),title="test_task_dtseries:  bold(:, " + selection + ") - gs");
        end
        function test_physio_HRV(this)
            as = this.testObj;
            as.tags_user=stackstr(use_dashes=true);
            as.source_physio = "HRV";

            bold = as.task_niigz();
            HRV = mlraut.PhysioHRV(as, bold);
            physio = HRV.call();
            gs = as.global_signal;

            figure; plot(physio); title("physio");
            figure; plot(gs); title("gs");

            figure; histogram(physio); title("histogram(physio)");
            figure; histogram(gs); title("histogram(gs)");

            as.fit_power_law(x=physio,title="test_physio_HRV:  physio");
            as.fit_power_law(x=gs,title="test_physio_HRV:  gs");
        end
        function test_physio_iFV(this)
            as = this.testObj;
            as.tags_user=stackstr(use_dashes=true);
            as.source_physio = "iFV-brightest";

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

            as.fit_power_law(x=physio,title="test_physio_iFV: physio");
            as.fit_power_law(x=gs,title="test_physio_iFV:  gs");
            as.fit_power_law(x=(physio - gs),title="test_physio_iFV:  physio - gs");
        end
        function test_physio_ROI(this)
            as = this.testObj;
            as.tags_user=stackstr(use_dashes=true);
            as.source_physio = "ROI";

            bold = as.task_niigz();
            Roi = mlraut.PhysioRoi(as, bold, from_wmparc_indices=[1 4 5 24 43 44]);  % csf
            Roi.wmparc.view_qc(Roi.roi_mask);
            physio = Roi.call();
            gs = as.global_signal;

            figure; plot(physio); title("physio");
            figure; plot(gs); title("gs");
            figure; plot(physio - gs); title("physio - gs");

            figure; histogram(physio); title("histogram(physio)");
            figure; histogram(gs); title("histogram(gs)");
            figure; histogram(physio - gs); title("histogram(physio - gs)");

            as.fit_power_law(x=physio,title="test_physio_ROI: physio");
            as.fit_power_law(x=gs,title="test_physio_ROI:  gs");
            as.fit_power_law(x=(physio - gs),title="test_physio_ROI:  physio - gs");
        end
        function test_physio_RV(this)
            as = this.testObj;
            as.tags_user=stackstr(use_dashes=true);
            as.source_physio = "RV";

            bold = as.task_niigz();
            RV = mlraut.PhysioRV(as, bold);
            physio = RV.call();
            gs = as.global_signal;

            figure; plot(physio); title("physio");
            figure; plot(gs); title("gs");

            figure; histogram(physio); title("histogram(physio)");
            figure; histogram(gs); title("histogram(gs)");

            as.fit_power_law(x=physio,title="test_physio_RV:  physio");
            as.fit_power_law(x=gs,title="test_physio_RV:  gs");
        end
        function test_task_signal_mask(this)
            as = mlraut.AnalyticSignalHCP( ...
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

        function test_fit_power_law(this)
            as = mlraut.AnalyticSignalHCP( ...
                subjects={'995174'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                do_save=false, ...
                force_band=true, ...
                tags=stackstr(use_dashes=true), ...
                source_physio='none');

            as.fit_power_law();
        end

        function test_power_spectrum_analysis(this)
            %% https://claude.ai/chat/7c0ba283-3bc7-4938-9dec-8acd7bd25e7a

            % Generate sample data (replace this with your actual time series)
            t = 0:0.72:(1200*0.72);
            x = randn(size(t));  % Random noise (you'd use your actual data here)

            % Compute the Fourier transform
            N = length(x);
            X = fft(x);

            % Compute the power spectrum
            P = abs(X).^2 / N;

            % Compute the corresponding frequencies
            fs = 1 / (t(2) - t(1));  % Sampling frequency
            f = (0:N-1)*(fs/N);      % Frequency range

            % Use only the first half of the spectrum (it's symmetric)
            P = P(1:floor(N/2)+1);
            f = f(1:floor(N/2)+1);

            % Plot the power spectrum on a log-log scale
            figure;
            loglog(f, P);
            xlabel('Frequency (Hz)');
            ylabel('Power');
            title('Power Spectrum');
            grid on;

            % Optional: Fit a power law
            % Select a range for fitting (adjust as needed)
            fit_range = f > 0.01 & f < 0.1;

            % Perform linear regression on log-log data
            p = polyfit(log10(f(fit_range)), log10(P(fit_range)), 1);

            % Add the fit line to the plot
            hold on;
            loglog(f(fit_range), 10.^(polyval(p, log10(f(fit_range)))), 'r--', 'LineWidth', 2);
            legend('Data', sprintf('Fit: slope = %.2f', p(1)));

            % Display the slope (which is the power law exponent)
            fprintf('Power law exponent: %.2f\n', p(1));
        end

        function test_complex_time_series_vis(this)
            %% https://claude.ai/chat/7c0ba283-3bc7-4938-9dec-8acd7bd25e7a

            % Generate a complex oscillatory time series
            t = linspace(0, 10, 1000);
            f1 = 1.0;
            f2 = 2.0;
            z = exp(1i * 2 * pi * f1 * t) + 0.5 * exp(1i * 2 * pi * f2 * t);

            % Create the 3D figure
            figure('Position', [100, 100, 1200, 500]);

            % 3D Line Plot
            subplot(1, 2, 1);
            civ = cividis;
            plot3(t, real(z), imag(z), LineWidth=2, Color=civ(1,:));
            xlabel('Time');
            ylabel('Real Part');
            zlabel('Imaginary Part');
            title('Complex Time Series: 3D View');
            grid on;

            % Add a 2D projection onto the complex plane
            subplot(1, 2, 2);
            scatter(real(z), imag(z), [], t, 'filled', 'o', MarkerFaceAlpha=0.8);
            xlabel('Real Part');
            ylabel('Imaginary Part');
            title('Complex Time Series: Complex Plane Projection');
            axis equal;
            colorbar;
            colormap('cividis');
            c = colorbar;
            c.Label.String = 'Time';

            % Adjust the layout
            fontsize(scale=1.5)
            sgtitle('Complex Oscillatory Time Series Visualization');
        end

        function test_late_hilbert_transform(this)

            pwd0 = pushd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200');

            as = mlraut.AnalyticSignalHCP( ...
                subjects={'995174'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                do_global_signal_regression=true, ...
                do_save=false, ...
                force_band=false, ...
                hp_thresh=0.005, ...
                lp_thresh=0.1, ...
                tags=stackstr(use_dashes=true), ...
                source_physio="iFV-brightest", ...

            as.current_subject = as.subjects{1};
            if ~contains(as.out_dir, as.current_subject)
                as.out_dir = fullfile(as.out_dir, as.current_subject);
                ensuredir(as.out_dir);
            end

            as.current_task = as.tasks{1};

            %as.plot3()

            %% parameters for aufbau

            select_network_type = "cortical";
            select_rsn = 'default mode';
            NETWORKS_YEO_NAMES = ...
                {'visual', 'somatomotor', 'dorsal attention', 'ventral attention', 'limbic', ...
                'frontoparietal', 'default mode', ...
                'task+', 'task-'};
            bool_select_rsn = contains(NETWORKS_YEO_NAMES, select_rsn);
            num_frames = 417;  % as.num_frames;

            tf = (num_frames - 1)*as.tr;
            t = 0:as.tr:tf;

            %% Late Hilbert transform

            ctx = as.average_network_signal(as.task_dtseries(), network_type="cortical");
            bold_gsr =  ...
                as.build_centered_and_rescaled( ...
                as.build_global_signal_regressed(ctx(:, bool_select_rsn)));

            bold_ = as.build_band_passed(bold_gsr);
            figure; plot(t, bold_(1:num_frames))
            ylabel("real(BOLD)")
            xlabel("time / s")

            bold_signal_ = ...
                hilbert(bold_);
            bold_signal__ = ...                
                bold_signal_(1:num_frames, :);
            as.plot3(num_frames=num_frames, t=t, z=bold_signal__, title="")
            as.fit_power_law(t=t, x=bold_signal__, title="power law:  BOLD cortical " + select_rsn)

            physio_gsr = ...
                as.build_centered_and_rescaled( ...
                as.build_global_signal_regressed(as.task_physio()));

            physio_ = as.build_band_passed(physio_gsr);
            figure; plot(t, physio_(1:num_frames))
            ylabel("real(arousal iFV)")
            xlabel("time / s")

            physio_signal_ = ...
                hilbert(physio_);
            physio_signal__ = ...
                physio_signal_(1:num_frames, :);
            as.plot3(num_frames=num_frames, t=t, z=physio_signal__, title="")
            as.fit_power_law(t=t, x=physio_signal__, title="power law:  Arousal iFV")

            analytic_signal_ = conj(physio_signal_).*bold_signal_;

            as.plot3(num_frames=num_frames, t=t, z=analytic_signal_, title="")
            as.fit_power_law(t=t, x=analytic_signal_, title="power law:  Analytic cortical " + select_rsn)

            popd(pwd0);
        end
        function test_mix_physio(this)
            p_0 = sin(0:0.1:2*pi);
            p_1 = 0.01*cos(0:0.2:4*pi);

            figure
            hold on
            for f = 0:.1:1
                this.testObj.frac_ext_physio = f;
                p = this.testObj.mix_physio(p_0, p_1);
                if f == 0.5
                    plot(p, LineWidth=3)
                else
                    plot(p)
                end
            end
        end
    end
    
    methods (TestClassSetup)
        function setupAnalyticSignal(this)
        end
    end
    
    methods (TestMethodSetup)
        function setupAnalyticSignalTest(this)
            cd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200');
            this.testObj = mlraut.AnalyticSignalHCP( ...
                subjects={'995174'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                do_7T=false, ...
                do_resting=true, ...
                do_task=false, ...
                do_save=false, ...
                do_save_dynamic=false, ...
                force_band=false, ...
                hp_thresh=0.01, ...
                lp_thresh=0.1, ...
                tags=stackstr(use_dashes=true));
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
