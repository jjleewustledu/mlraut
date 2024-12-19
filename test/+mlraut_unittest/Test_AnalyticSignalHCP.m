classdef Test_AnalyticSignalHCP < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 11-Nov-2024 23:22:32 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/test/+mlraut_unittest.
    %  Developed on Matlab 24.1.0.2689473 (R2024a) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    
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
        
        function test_call_no_physio(this)
            as = this.testObj;
            as.do_save=true;
            as.do_save_dynamic=true;
            as.do_save_ciftis=true;
            as.do_plot_networks=true;
            as.source_physio="none";
            
            disp(as)            
            call(as);
        end

        function test_call_iFV(this)
            as = this.testObj;
            as.do_save=true;
            as.do_save_dynamic=true;
            as.do_save_ciftis=true;
            as.do_plot_networks=true;
            as.source_physio="iFV";
           
            disp(as)            
            call(as);

            % qc
            zeta = as.HCP_signals.ctx.psi(:,9) ./ as.HCP_signals.ctx.phi(:,9);
            as.plot3(z=zeta)  % re(psi) vaguely resemble ECG :-)
            % as.plot3(z=mean(as.bold_signal, 2))
            as.plot3(z=as.HCP_signals.ctx.psi(:,9), symbol="\psi")  % ctx, task-
            as.plot3(z=as.HCP_signals.ctx.phi(:,9), symbol="\phi")  % ctx, task-
            figure; imagesc(angle(as.physio_signal));
        end

        function test_call_HRV(this)
            as = this.testObj;
            as.do_save=true;
            as.do_save_dynamic=true;
            as.do_save_ciftis=true;
            as.do_plot_networks=true;
            as.source_physio="HRV";
            
            disp(as)            
            call(as);
        end

        function test_call_RV(this)
            as = this.testObj;
            as.do_save=true;
            as.do_save_dynamic=true;
            as.do_save_ciftis=true;
            as.do_plot_networks=true;
            as.source_physio="RV";
            
            disp(as)            
            call(as);
        end

        function test_call_wmparc(this)

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
                as.do_save=true;
                as.do_save_dynamic=true;
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
            this.testObj_ = [];
        end
    end
    
    methods (TestMethodSetup)
        function setupAnalyticSignalHCPTest(this)
            cd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200');
            this.testObj = mlraut.AnalyticSignalHCP( ...
                subjects={'995174'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                do_7T=false, ...
                do_resting=true, ...
                do_task=false, ...
                do_save=false, ...
                do_save_dynamic=false, ...
                do_save_ciftis=false, ...
                final_normalization="none", ...
                force_band=false, ...
                hp_thresh=0.01, ...
                lp_thresh=0.1, ...
                v_physio=50, ...
                plot_range=1:250, ...
                global_signal_regression=true, ...
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
