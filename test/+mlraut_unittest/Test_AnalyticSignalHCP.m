classdef Test_AnalyticSignalHCP < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 11-Nov-2024 23:22:32 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/test/+mlraut_unittest.
    %  Developed on Matlab 24.1.0.2689473 (R2024a) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    
    properties
        testObj
    end    

    methods (Static)
        function testObj = called(testObj)
            call_subject(testObj);
        end
    end

    methods (Test)
        function test_afun(this)
            import mlraut.*
            this.assumeEqual(1,1);
            this.verifyEqual(1,1);
            this.assertEqual(1,1);
        end

        function test_bin_by_physio_angle(this)
            out_dir = "/Volumes/PrecunealSSD2/AnalyticSignalHCP";
            mat = fullfile(out_dir, ...
                "sub-996782_ses-rfMRI-REST1-RL_proc-RV-gsr1-butter2-lp0p05-hp0p01-scaleiqr-Test-AnalyticSignalHCP-setupAnalyticSignalHCP.mat");            
            pwd0 = pushd(out_dir);

            ld = load(mat);
            ihcp = ld.this;
            phi = ihcp.physio_signal;
            psi = ihcp.bold_signal;
            c = mlraut.Cifti(ihcp);
            binned = c.bin_by_physio_angle(psi, phi);
            rbinned = single(real(binned));

            prev_cii = cifti_read("AnalyticSignalHCPPar_mean_bold_complexavg.dtseries.nii");
            this.verifyEqual(rbinned, prev_cii.cdata', RelTol=1e-6);

            popd(pwd0);
        end

        function test_call_20250529(this)
            tic
            out_dir = "/Volumes/PrecunealSSD2/AnalyticSignalHCP";
            as = mlraut.AnalyticSignalHCPPar( ...
                subjects={'996782'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                do_save=true, ...
                do_save_subset=true, ...
                do_save_ciftis=true, ...
                hp_thresh=0.01, ...
                lp_thresh=0.1, ...
                out_dir=out_dir, ...
                source_physio=[ ...
                "iFV-brightest", "iFV-quantile", "sFV", "3rdV", "latV", "centrumsemiovale", "ctx", "RV", "HRV"], ...
                tags=stackstr(use_dashes=true));
            call(as);
            toc

            [msg,id] = lastwarn();
            disp(msg)
            disp(id)

            %% real memory max < 34 GB; Elapsed time is 120 seconds; 1.69 GB mat file
        end

        function test_ctor(this)
            as = this.testObj;
            this.verifyEqual(as.num_nets, 9);
            this.verifyEqual(as.num_sub, 1);
            this.verifyEqual(as.num_tasks, 1);
        end

        function test_task_mask_niigz(this)
            as = this.testObj;
            ic = as.task_mask_niigz;
            % as.task_ref_niigz.view(ic);

            this.verifyEqual(ic.filename, "wmparc.2_binarized.nii.gz")
            this.verifyEqual(ic.qfac, -1)
            this.verifyEqual(size(ic), [91, 109, 91])
            this.verifyInstanceOf(ic.imagingFormat.img, "single")
            this.verifyEqual(dipmax(ic), 1)
            this.verifyEqual(dipsum(ic), 195505)
        end
        
        function test_task_ref_niigz(this)
            as = this.testObj;
            ic = as.task_ref_niigz;
            % ic.view_qc(as.task_mask_niigz);

            this.verifyEqual(ic.filename, "rfMRI_REST1_RL_SBRef.nii.gz")
            this.verifyEqual(ic.qfac, -1)
            this.verifyEqual(size(ic), [91, 109, 91])
            this.verifyInstanceOf(ic.imagingFormat.img, "single")
            this.verifyEqual(dipmax(ic), 49664.171875, AbsTol=1e-3)
            this.verifyEqual(dipsum(ic), 2.286530000098282e+09, AbsTol=1)
        end

        function test_templates(this)
            as = this.testObj;
            as.current_subject = '995174';
            as.current_task = 'rfMRI_REST1_RL';
            
            % template_niigz ~ wmparc
            this.verifyTrue(isfile(as.wmparc_fqfn))
            this.verifyInstanceOf(as.template_niigz, "mlfourd.ImagingContext2")
            this.verifyEqual(as.template_niigz.filename, "wmparc.2.nii.gz")

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
        end

        function test_average_network_signals(this)
            %% candidate supplemental figure, with concat runs

            as = this.testObj;
            call_subject(as);

            yeos = as.HCP_signals;
            yeo_names = mlraut.NetworkData.NETWORKS_YEO_NAMES;
            for anat = ["cbm", "ctx", "str", "thal"]
                figure
                hold on
                for y = 1:length(yeo_names)
                    plot(real(yeos.(anat).psi(:, y))); 
                end
                title(sprintf("yeos.%s", anat));  
                legend(yeo_names)
                hold off

                figure;     
                tiledlayout(3,3);
                for y = 1:length(yeo_names)
                    nexttile
                    histogram(real(yeos.(anat).psi(:, y))); 
                    title(sprintf("yeos(%s)", yeo_names{y}), 100);
                    title(sprintf("yeos.%s(%s)", anat, yeo_names{y}));    
                end
            end
        end

        function test_rescaled(this)
            hcp_ = mlraut.AnalyticSignalHCP( ...
                subjects={'996782'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                hp_thresh=0.01, ...
                lp_thresh=0.1, ...
                filter_order=8, ...
                source_physio="RV", ...
                tags=stackstr(use_dashes=true));
            malloc(hcp_);

            rescalings = ["none", "1-norm", "2-norm", "fro", "inf-norm", "iqr", "mean_ad", "median_ad", "std"];
            N = length(rescalings);
            colors = lines(N);  % or use other colormaps like hsv(N), parula(N), etc.
            ridx = 0;
            Nt = hcp_.num_frames;
            times = linspace(0, (Nt-1)*hcp_.tr, Nt);

            figure;
            hold on;
            for rescaling = rescalings

                ridx = ridx + 1;

                tic
                hcp_.rescaling = rescaling;
                bold_gsr_ = ...
                    hcp_.build_global_signal_regressed(hcp_.task_dtseries());
                bold_ = ...
                    hcp_.build_rescaled( ...
                    hcp_.build_band_passed( ...
                    hcp_.build_centered(bold_gsr_)));
                bold_ = mean(bold_, 2);
                bold_ = log10(abs(bold_));
                physio_ = hcp_.task_physio();
                physio_ = log10(abs(physio_));
                fprintf("rescaling %s: ", rescaling); toc

                % Plot physio (dashed line)
                plot(times, physio_, '--', 'Color', colors(ridx,:), 'LineWidth', 1.5);

                % Plot bold (solid line)
                plot(times, bold_, '-', 'Color', colors(ridx,:), 'LineWidth', 1.5);

            end
            hold off;

            % Create labels, legend entries
            ylabel("log_{10} physio, log_{10} bold")
            xlabel("time (s)")
            legend_entries = cell(1, 2*N);
            for i = 1:N
                legend_entries{2*i-1} = sprintf('%s (physio)', rescalings(i));
                legend_entries{2*i} = sprintf('%s (bold)', rescalings(i));
            end
            legend(legend_entries);
        end






        function test_fultz_iFV_brightest(this)
            as = mlraut.AnalyticSignalHCP( ...
                subjects={'996782'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                hp_thresh=0.01, ...
                lp_thresh=0.1, ...
                filter_order=8, ...
                source_physio="iFV-brightest", ...
                tags=stackstr(use_dashes=true));
            call_subject(as);
            
            tseries = "X";  % ["bold", "-dbold/dt", "X", "reY", "imY", "Z"];
            for t = tseries
                as.plot_coherencyc(tseries=t);
            end
        end

        function test_fultz_iFV(this)
            as = mlraut.AnalyticSignalHCP( ...
                subjects={'996782'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                hp_thresh=0.01, ...
                lp_thresh=0.1, ...
                filter_order=8, ...
                source_physio="iFV", ...
                tags=stackstr(use_dashes=true));
            call_subject(as);
            
            tseries = ["bold", "-dbold/dt", "X", "reY", "imY", "Z"];
            for t = tseries
                as.plot_coherencyc(tseries=t);
            end
        end

        function test_fultz_sFV(this)
            as = this.testObj;
            as.source_physio = "sFV";
            call_subject(as);
            
            tseries = ["bold", "-dbold/dt", "X", "reY", "imY", "Z"];
            for t = tseries
                as.plot_coherencyc(tseries=t);
            end
        end

        function test_fultz_latV(this)
            as = this.testObj;
            as.source_physio = "latV";
            call_subject(as);
            
            tseries = ["bold", "-dbold/dt", "X", "reY", "imY", "Z"];
            for t = tseries
                as.plot_coherencyc(tseries=t);
            end
        end

        function test_fultz_csf(this)
            as = this.testObj;
            as.source_physio = "csf";
            call_subject(as);
            
            tseries = ["bold", "-dbold/dt", "X", "reY", "imY", "Z"];
            for t = tseries
                as.plot_coherencyc(tseries=t);
            end
        end

        function test_fultz_HRV(this)
            as = mlraut.AnalyticSignalHCP( ...
                subjects={'996782'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                hp_thresh=0.01, ...
                lp_thresh=0.1, ...
                filter_order=8, ...
                source_physio="HRV", ...
                tags=stackstr(use_dashes=true));
            call_subject(as);
                        
            tseries = ["bold", "-dbold/dt", "X", "reY", "imY", "Z"];
            for t = tseries
                as.plot_coherencyc(tseries=t);
            end
        end

        function test_fultz_RV(this)
            as = mlraut.AnalyticSignalHCP( ...
                subjects={'996782'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                hp_thresh=0.01, ...
                lp_thresh=0.1, ...
                filter_order=8, ...
                source_physio="RV", ...
                tags=stackstr(use_dashes=true));
            call_subject(as);
            
            tseries = ["bold", "-dbold/dt", "X", "reY", "imY", "Z"];
            for t = tseries
                as.plot_coherencyc(tseries=t);
            end
        end

        function test_plot_coherency_par(this)
            as = mlraut.AnalyticSignalHCPPar( ...
                subjects={'996782'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                hp_thresh=0.01, ...
                lp_thresh=0.1, ...
                filter_order=8, ...
                source_physio="iFV-brightest", ...
                tags=stackstr(use_dashes=true));
            call_subject(as);

            psi_ld = load("/Volumes/PrecunealSSD2/AnalyticSignalHCP/FultzMulti_build_tseries_phi_psi_from_mat159_iFV-brightest_psi.mat");
            phi_ld = load("/Volumes/PrecunealSSD2/AnalyticSignalHCP/FultzMulti_build_tseries_phi_psi_from_mat159_iFV-brightest_phi.mat");
            psi = psi_ld.psi;
            phi = phi_ld.phi;
            % psi = [];
            % phi = [];
            
            tseries = "X";  % ["bold", "-dbold/dt", "X", "reY", "imY", "Z"];
            for t = tseries
                as.plot_coherencyc( ...
                    psi, ...
                    phi, ...
                    tseries=t);
            end
        end



        

        function test_call_wmparc(this)

            return

            % wmparc = 'precuneus';
            % wmparc = 'posteriorcingulate';
            % wmparc = 'hippocampus';
            % wmparc = 'entorhinal';
            % wmparc = 'medialorbitofrontal';
            % wmparc = 'insula';
            % wmparc = 'cuneus';
            % wmparc = 'thalamus';
            % wmparc = 'caudate';
            % wmparc = 'putamen';
            % wmparc = 'pallidum';
            % wmparc = 'cerebellum';
            % wmparc = 'ponsvermis';
            % wmparc = 'brainstem';
            % wmparc = 'brainstem+';
            % wmparc = 'csf';
            % wmparc = 'centrumsemiovale';
            % wmparc = 'corpuscallosum';

            wmparcs = { ...
                'cuneus' 'corpuscallosum' ...
                'precuneus' 'posteriorcingulate' 'hippocampus' 'entorhinal' 'medialorbitofrontal' ...
                'insula' ...
                'thalamus' 'caudate' 'putamen' 'pallidum' 'cerebellum' ...
                'ponsvermis' 'brainstem' 'brainstem+' 'csf' 'centrumsemiovale'};

            for w = wmparcs
                wmparc = w{1};

                as = this.testObj;
                as.do_save=false;
                as.do_save_dynamic=false;
                as.do_save_ciftis=true;
                as.source_physio=wmparc;
                as.out_dir = sprintf('/Volumes/PrecunealSSD2/AnalyticSignalHCP/physio_%s', wmparc);

                disp(as)
                call(as);
            end
        end
    end
    
    methods (TestClassSetup)
        function setupAnalyticSignalHCP(this)
            import mlraut.*
            % cd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200');
            this.testObj_ = AnalyticSignalHCP( ...
                subjects={'996782'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                hp_thresh=0.01, ...
                lp_thresh=0.05, ...
                filter_order=2, ...
                source_physio="RV", ...
                tags=stackstr(use_dashes=true));
        end
    end
    
    methods (TestMethodSetup)
        function setupAnalyticSignalHCPTest(this)
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
