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
            as = this.testObj;  % AnalyticSignal

            this.verifyEqual(as.current_subject, '995174')
            this.verifyEqual(as.current_task, 'rfMRI_REST1_RL')
            this.verifyEqual(as.source_physio, "iFV-brightest")
            this.verifyEqual(as.v_physio, 50)
            this.verifyEqual(as.anatomy_list, {'ctx', 'str', 'thal', 'cbm'})
            this.verifyEqual(dipmax(as.global_signal), 10369.6112020772, AbsTol=1e-3)
            this.verifyEqual(as.hp_thresh, [])
            this.verifyEqual(as.lp_thresh, 0.1)
            this.verifyEqual(size(as.bold_signal), [1196, 91282])
            this.verifyEqual(size(as.physio_signal), [1196, 91282])
            this.verifyEqual(as.roi, [])
        end

        function test_global_signal(this)
            %% check magnitudes of {ctx, cbv, str, thal} vs gsr

            as = this.testObj;

            ctx = as.average_anat_signal(as.task_dtseries, network_type="ctx");
            crb = as.average_anat_signal(as.task_dtseries, network_type="crb");
            str = as.average_anat_signal(as.task_dtseries, network_type="str");
            thal = as.average_anat_signal(as.task_dtseries, network_type="thal");
            % plot([as.global_signal, ctx, crb, str, thal]);
            % legend(["global signal", "ctx", "crb", "str", "thal"]);

            this.verifyEqual(mean(as.global_signal, "all"), 10332.1418491366, AbsTol=1);
            this.verifyEqual(mean(ctx, "all"), single(1.0892601e+04), AbsTol=1);
            this.verifyEqual(mean(crb, "all"), single(8.5629580e+03), AbsTol=1);
            this.verifyEqual(mean(str, "all"), single(1.1747824e+04), AbsTol=1);
            this.verifyEqual(mean(thal, "all"), single(1.2746612e+04), AbsTol=1);
        end

        function test_build_global_signal_regressed(this)
            %% check magnitudes of task_niigz vs gsr

            as = this.testObj;

            ifc = as.task_niigz.imagingFormat;
            gsr = copy(ifc);
            gsr.fileprefix = ifc.fileprefix + "_gsr";
            gs = reshape(as.global_signal, [1, 1, 1, as.num_frames]);
            gsr.img = ifc.img - gs;
            % gsr.view_qc(ifc);

            % native BOLD values
            this.verifyEqual(dipmin(ifc.img), -1701.63549804688, AbsTol=1)
            this.verifyEqual(mean(ifc.img, "all"), single(2.5297590e+03), AbsTol=1)
            this.verifyEqual(dipmax(ifc.img), 26076.93359375, AbsTol=1)

            % gsr shifts native BOLD values down
            this.verifyEqual(dipmin(gsr.img), -12022.9560546875, AbsTol=1)
            this.verifyEqual(mean(gsr.img, "all"), single(-7.8023848e+03), AbsTol=1)
            this.verifyEqual(dipmax(gsr.img), 15745.298828125, AbsTol=1)
        end

        function test_build_band_passed(this)
            as = this.testObj;
            as.force_legacy_butter = false;
            gs = as.build_centered_and_rescaled(as.global_signal);
            gs1 = gs;

            orders = [2,4,8,16];
            for col_idx = 2:5
                as.filter_order = orders(col_idx-1);
                gs1(:,col_idx) = as.build_band_passed(gs, do_reset=true);
                this.verifyInstanceOf(as.digital_filter, "digitalFilter")
            end
            TT = array2timetable(gs1, Timestep=seconds(as.tr));
            signalAnalyzer(TT)
        end

        function test_mat_fqfn(this)
            this.verifyEqual(mybasename(this.testObj.mat_fqfn), ...
                "sub-995174_ses-rfMRI-REST1-RL_proc-v50-iFV-brightest-scaleiqr-Test-AnalyticSignal-setupAnalyticSignal")
        end

        function test_PhysioHRV(this)
            as = this.testObj;
            as.tags_user=stackstr(use_dashes=true);
            as.source_physio = "HRV";

            bold = as.task_niigz();
            HRV = mlraut.PhysioHRV(as, bold);
            physio = HRV.call();
            gs = as.build_centered_and_rescaled(as.global_signal);

            figure; plot(physio); title("physio");
            figure; plot(gs); title("gs");

            figure; histogram(physio); title("histogram(physio)");
            figure; histogram(gs); title("histogram(gs)");

            as.fit_power_law(x=physio,title="test_physio_HRV:  physio");
            as.fit_power_law(x=gs,title="test_physio_HRV:  gs");
        end

        function test_PhysioRV(this)
            as = this.testObj;
            as.tags_user=stackstr(use_dashes=true);
            as.source_physio = "RV";

            bold = as.task_niigz();
            RV = mlraut.PhysioRV(as, bold);
            physio = RV.call();
            gs = as.build_centered_and_rescaled(as.global_signal);

            figure; plot(physio); title("physio");
            figure; plot(gs); title("gs");

            figure; histogram(physio); title("histogram(physio)");
            figure; histogram(gs); title("histogram(gs)");

            as.fit_power_law(x=physio,title="test_physio_RV:  physio");
            as.fit_power_law(x=gs,title="test_physio_RV:  gs");
        end

        function test_IFourthVentricle(this)
            as = this.testObj;
            as.tags_user=stackstr(use_dashes=true);
            as.source_physio = "iFV-brightest";

            bold = as.task_niigz();
            iFV = mlraut.IFourthVentricle(as, bold);
            iFV.wmparc.view_qc(iFV.ifv_mask);
            physio = iFV.call();
            gs = as.build_centered_and_rescaled(as.global_signal);

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

        function test_task_physio_supplementary(this)
            as = this.testObj;
            as.tags_user=stackstr(use_dashes=true);
            as.source_physio = "iFV-brightest";
            as.source_physio_supplementary = [ ...
                "iFV-quantile", "sFV", "3rdV", "latV", "csf", "centrumsemiovale", "ctx"];
            physio = as.task_physio();
            cmap = as.task_physio_supplementary();
            keys = cmap.keys;

            figure;
            plot(physio(:, 1))
            hold on
            for cidx = 1:length(cmap)                
                mat = cmap.values(keys(cidx));
                plot(mat{1}(:, 1));
            end
            hold off
            legend([as.source_physio, cmap.keys])

            figure;
            pmtm(physio(:, 1))
            hold on
            for cidx = 1:length(cmap)                
                mat = cmap.values(keys(cidx));
                pmtm(mat{1}(:, 1));
            end
            hold off
            legend([as.source_physio, cmap.keys])
        end








        function test_power_spectrum_analysis(this)
            %% https://claude.ai/chat/7c0ba283-3bc7-4938-9dec-8acd7bd25e7a

            return

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

            return

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

            return

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
                source_physio="iFV-brightest");

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
            ylabel("real(arousal iFV-brightest)")
            xlabel("time / s")

            physio_signal_ = ...
                hilbert(physio_);
            physio_signal__ = ...
                physio_signal_(1:num_frames, :);
            as.plot3(num_frames=num_frames, t=t, z=physio_signal__, title="")
            as.fit_power_law(t=t, x=physio_signal__, title="power law:  Arousal iFV-brightest")

            analytic_signal_ = conj(physio_signal_).*bold_signal_;

            as.plot3(num_frames=num_frames, t=t, z=analytic_signal_, title="")
            as.fit_power_law(t=t, x=analytic_signal_, title="power law:  Analytic cortical " + select_rsn)

            popd(pwd0);
        end
        
        function test_mix_physio(this)
            return

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
            import mlraut.*
            cd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200');
            this.testObj_ = AnalyticSignal( ...
                subjects={'995174'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                hp_thresh=[], ...
                lp_thresh=0.1, ...
                source_physio="iFV-brightest", ...
                tags=stackstr(use_dashes=true));
        end
    end
    
    methods (TestMethodSetup)
        function setupAnalyticSignalTest(this)
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
