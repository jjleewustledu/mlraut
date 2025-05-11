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

        function test_templates(this)
            as = this.testObj;

            % template_cifti ~ thickness
            this.verifyTrue(isfile(as.thickness_dscalar_fqfn));
            this.verifyTrue(isstruct(as.template_cifti.metadata));
            this.verifyTrue(iscell(as.template_cifti.diminfo));
            this.verifyEqual(as.template_cifti.diminfo{1}.models{1}.count, 149141);
            this.verifyEqual(as.template_cifti.diminfo{1}.models{1}.struct, 'CORTEX_LEFT');
            this.verifyEqual(as.template_cifti.diminfo{1}.models{1}.type, 'surf');
            this.verifyEqual(size(as.template_cifti.diminfo{1}.models{1}.vertlist), [1, 149141]);
            this.verifyEqual(as.template_cifti.diminfo{1}.models{2}.count, 149120);
            this.verifyEqual(as.template_cifti.diminfo{1}.models{2}.struct, 'CORTEX_RIGHT');
            this.verifyEqual(as.template_cifti.diminfo{1}.models{2}.type, 'surf');
            this.verifyEqual(size(as.template_cifti.diminfo{1}.models{2}.vertlist), [1, 149120]);
            this.verifyEqual(as.template_cifti.diminfo{2}.maps.name, 'HCA9992517_V1_MR_Thickness')
            this.verifyEqual(size(as.template_cifti.cdata), [298261, 1]);
            
            % template_niigz ~ wmparc
            this.verifyTrue(isfile(as.wmparc_fqfn))
            this.verifyInstanceOf(as.template_niigz, "mlfourd.ImagingContext2")
            this.verifyEqual(as.template_niigz.filename, "wmparc.2.nii.gz")
        end

        function test_memory_footprint(this)
            tic
            as = this.testObj;
            as.do_save=false;
            as.do_save_dynamic=false;
            as.do_save_ciftis=false;
            as.do_plot_networks=false;
            as.source_physio="iFV-brightest";
            as.out_dir = '/Volumes/PrecunealSSD2/AnalyticSignalHCPAging';
            
            disp(as)            
            call(as);
            toc
        end

        function test_call_iFV(this)
            as = this.testObj;
            as.do_save=false;
            as.do_save_ciftis=false;
            as.source_physio="iFV-brightest";
            % as.force_legacy_butter = true;  % legacy butter has less contrast for X(rsn7)
            as.out_dir = '/Volumes/PrecunealSSD2/AnalyticSignalHCPAging';
            
            tic
            disp(as)
            call(as)
            toc

            % qc
            % pwd0 = pushd(as.out_dir);
            % zeta = as.HCP_signals.ctx.psi(:,9) ./ as.HCP_signals.ctx.phi(:,9);
            % as.fit_power_law(x=zeta, title="\zeta = \psi / \phi");
            % as.plot3(z=zeta, symbol="\zeta")  % re(psi) vaguely resemble ECG :-)
            % as.plot3(z=as.HCP_signals.ctx.psi(:,9), symbol="\psi")  % ctx, task-
            % as.plot3(z=as.HCP_signals.ctx.phi(:,9), symbol="\phi")  % ctx, task-
            % % figure; imagesc(angle(as.physio_signal));
            % as.plotting.saveFigures("qc");
            % popd(pwd0);
        end

        function test_call_physio(this)

            for phys = {'HRV' 'RV'}
                as = this.testObj;
                as.do_save=true;
                as.do_save_ciftis=true;
                as.source_physio=phys{1};
                as.out_dir = '/Volumes/PrecunealSSD2/AnalyticSignalHCPAging';

                disp(as)
                call(as)
            end

            % Elapsed time is 273 seconds on twistor.  Peak memory is 64 GB.  Files saved < 8 GB.
            % Elapsed time is ~791 seconds on vglab2.
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
            % wmparc = 'brainstem';
            % wmparc = 'brainstem+';
            % wmparc = 'csf';
            % wmparc = 'centrumsemiovale';
            % wmparc = 'corpuscallosum';

            wmparcs = { ...
                'precuneus' 'posteriorcingulate' 'hippocampus' 'entorhinal' 'medialorbitofrontal' ...
                'insula' ...
                'cuneus' ...
                'thalamus' 'caudate' 'putamen' 'pallidum' 'cerebellum' ...
                'brainstem' 'brainstem+' 'csf' 'centrumsemiovale' 'corpuscallosum'};

            for w = wmparcs
                wmp = w{1};

                as = this.testObj;
                as.do_save=true;
                as.do_save_dynamic=true;
                as.do_save_ciftis=true;
                as.source_physio=wmp;
                as.out_dir = sprintf('/Volumes/PrecunealSSD2/AnalyticSignalHCPAging/physio_%s', wmp);

                disp(as)
                call(as);
            end
        end
    end
    
    methods (TestClassSetup)
        function setupAnalyticSignalHCPAging(this)
        end
    end
    
    methods (TestMethodSetup)
        function setupAnalyticSignalHCPAgingTest(this)
            cd('/Volumes/PrecunealSSD2/HCPAging/HCPAgingRec/fmriresults01');
            this.testObj = mlraut.AnalyticSignalHCPAging( ...
                subjects={'HCA9992517_V1_MR'}, ...
                tasks={'fMRI_CONCAT_ALL'}, ...
                do_save=false, ...
                do_save_dynamic=false, ...
                do_save_ciftis=false, ...
                hp_thresh=0.01, ...
                lp_thresh=0.1, ...
                v_physio=50, ...
                plot_range=1:225, ...
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
