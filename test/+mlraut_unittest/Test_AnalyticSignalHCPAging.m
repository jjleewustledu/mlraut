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

        function test_load(this)
            fqfn_mat = fullfile( ...
                "/Volumes/PrecunealSSD2/AnalyticSignalHCPAging/HCA9992517_V1_MR", ...
                "sub-HCA9992517_V1_MR_ses-fMRI-CONCAT-ALL_proc-iFV-brightest-gsr1-butter8-lp0p1-hpnone-scaleiqr-subset-Test-AnalyticSignalHCPAging-setupAnalyticSignalHCPAging.mat");
            ld = load(fqfn_mat);
            that_subset = ld.this_subset;
            that = mlraut.AnalyticSignalHCPAging.load(fqfn_mat);

            %% adapted to correspond to AnalyticSignalHCP.load()

            this.verifyEqual(that_subset.num_nets, that.num_nets);
            this.verifyEqual(that_subset.num_sub, that.num_sub);
            this.verifyEqual(that_subset.num_tasks, that.num_tasks);
            this.verifyEqual(that_subset.comparator, that.comparator);
            this.verifyEqual(that_subset.HCP_signals, that.HCP_signals);
            this.verifyEqual(that_subset.do_global_signal_regression, that.do_global_signal_regression);
            this.verifyEqual(that_subset.do_plot_emd, that.do_plot_emd);
            this.verifyEqual(that_subset.do_plot_global_physio, that.do_plot_global_physio);
            this.verifyEqual(that_subset.do_plot_networks, that.do_plot_networks);
            this.verifyEqual(that_subset.do_plot_radar, that.do_plot_radar);
            this.verifyEqual(that_subset.do_plot_wavelets, that.do_plot_wavelets);
            this.verifyEqual(that_subset.do_save, that.do_save);
            this.verifyEqual(that_subset.do_save_ciftis, that.do_save_ciftis);
            % this.verifyEqual(that_subset.do_save_ciftis_alt, that.do_save_ciftis_alt);  % bug from 4/17/2025: missing in save_subset()
            this.verifyEqual(that_subset.do_save_ciftis_of_diffs, that.do_save_ciftis_of_diffs);
            this.verifyEqual(that_subset.do_save_dynamic, that.do_save_dynamic);
            this.verifyEqual(that_subset.do_save_regularized, that.do_save_regularized);
            this.verifyEqual(that_subset.do_save_subset, that.do_save_subset);
            this.verifyEqual(that_subset.filter_order, that.filter_order);
            this.verifyEqual(that_subset.force_band, that.force_band);
            this.verifyEqual(that_subset.force_legacy_butter, that.force_legacy_butter);
            this.verifyEqual(that_subset.frac_ext_physio, that.frac_ext_physio);
            this.verifyEqual(that_subset.scale_of_rescaling, that.scale_of_rescaling, RelTol=1e-5);
            this.verifyEqual(that_subset.source_physio, that.source_physio);
            this.verifyEqual(that_subset.source_physio_supplementary, that.source_physio_supplementary);
            this.verifyEqual(that_subset.v_physio, that.v_physio);
            this.verifyEqual(that_subset.anatomy_list, that.anatomy_list);
            this.verifyEqual(that_subset.digital_filter.Coefficients, that.digital_filter.Coefficients);
            this.verifyEqual(that_subset.global_signal, that.global_signal);
            % this.verifyEqual(that_subset.hp_thresh, that.hp_thresh);  % bug from 4/17/2025: hp_thresh <- lp_thresh 
            % this.verifyEqual(that_subset.lp_thresh, that.lp_thresh);
            this.verifyEqual(that_subset.rescaling, that.rescaling);
            this.verifyEqual(that_subset.rsn_list, that.rsn_list);
            this.verifyEqual(that_subset.tags_user, that.tags_user);
            this.verifyEqual(that_subset.bold_signal, that.bold_signal);
            this.verifyEqual(that_subset.physio_angle, that.physio_angle);
            this.verifyEqual(that_subset.physio_signal, that.physio_signal);
            this.verifyEqual(that_subset.physio_supplementary, that.physio_supplementary);
            this.verifyEqual(that_subset.roi, that.roi);
            this.verifyEqual(that_subset.v_physio_is_inf, that.v_physio_is_inf);
            this.verifyEqual(that_subset.do_7T, that.do_7T);
            this.verifyEqual(that_subset.do_resting, that.do_resting);
            this.verifyEqual(that_subset.do_task, that.do_task);
            this.verifyEqual(that_subset.max_frames, that.max_frames);
            this.verifyEqual(that_subset.current_subject, that.current_subject);
            this.verifyEqual(that_subset.current_task, that.current_task);
            this.verifyEqual(that_subset.subjects, that.subjects);
            this.verifyEqual(that_subset.tasks, that.tasks);
            this.verifyEqual(that_subset.extended_task_dir, that.extended_task_dir);
            this.verifyEqual(that_subset.Fs, that.Fs);
            this.verifyEqual(that_subset.num_frames, that.num_frames);
            this.verifyEqual(that_subset.num_frames_ori, that.num_frames_ori);
            this.verifyEqual(that_subset.num_frames_to_trim, that.num_frames_to_trim);
            this.verifyEqual(that_subset.num_nodes, that.num_nodes);
            this.verifyEqual(that_subset.out_dir, that.out_dir);
            this.verifyEqual(that_subset.root_dir, that.root_dir);
            this.verifyEqual(that_subset.stats_fqfn, that.stats_fqfn);
            this.verifyEqual(that_subset.task_dir, that.task_dir);
            this.verifyEqual(that_subset.task_dtseries_fqfn, that.task_dtseries_fqfn);
            this.verifyEqual(that_subset.task_niigz_fqfn, that.task_niigz_fqfn);
            this.verifyEqual(that_subset.task_ref_niigz_fqfn, that.task_ref_niigz_fqfn);
            this.verifyEqual(that_subset.task_ref_dscalar_fqfn, that.task_ref_dscalar_fqfn);
            this.verifyEqual(that_subset.thickness_dscalar_fqfn, that.thickness_dscalar_fqfn);
            this.verifyEqual(that_subset.t1w_fqfn, that.t1w_fqfn);
            this.verifyEqual(that_subset.tr, that.tr);
            this.verifyEqual(that_subset.waves_dir, that.waves_dir);
            this.verifyEqual(that_subset.wmparc_fqfn, that.wmparc_fqfn);
            this.verifyEqual(that_subset.workbench_dir, that.workbench_dir);

            % this.verifyEqual(that_subset.class, class(that));  % TODO
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

        function test_call(this)
            tic
            as = this.testObj;
            as.out_dir = "/Volumes/PrecunealSSD2/AnalyticSignalHCPAging";
            as.do_save = false;
            call(as);
            toc

            %% real memory max < 53 GB; Elapsed time is 133 seconds
        end

        function test_memory_footprint(this)
            tic
            as = this.testObj;
            as.out_dir = "/Volumes/PrecunealSSD2/AnalyticSignalHCPAging";
            as.do_save = true;
            as.do_save_subset = true;
            call(as);
            toc

            %% real memory max < 53 GB; Elapsed time is 133 seconds; 3.83 GB mat file
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
            as.do_save_subset = true;
            call(as);
            toc

            %% real memory max < 46 GB; Elapsed time is 337 seconds; 3.83 GB mat file
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

        function test_exemplar_20250415(this)
            out_dir = '/Volumes/PrecunealSSD2/AnalyticSignalHCPAging/exemplar_20250415';
            ensuredir(out_dir);

            tic
            as = mlraut.AnalyticSignalHCPAging( ...
                subjects={'HCA9992517_V1_MR'}, ...
                tasks={'fMRI_CONCAT_ALL'}, ...
                source_physio="iFV-brightest", ...
                do_global_signal_regression=true, ...
                hp_thresh=[], ...
                lp_thresh=0.1, ...
                filter_order=8, ...
                out_dir=out_dir, ...
                do_save=true, ...
                do_save_ciftis=true, ...
                do_save_ciftis_alt=true, ...
                do_save_dynamic=true, ...
                do_save_regularized=true, ...
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
                hp_thresh=[], ...
                lp_thresh=0.1, ...
                source_physio="iFV-brightest", ...
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
