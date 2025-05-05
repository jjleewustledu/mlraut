classdef Test_Physio < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 12-Sep-2022 20:02:14 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/test/+mlraut_unittest.
    %  Developed on Matlab 9.12.0.2039608 (R2022a) Update 5 for MACI64.  Copyright 2022 John J. Lee.
    
    properties
        rsfMRI_scene
        testObj
    end
    
    methods (Test)
        function test_ctor(this)
            disp(this.testObj)
        end
        function test_afun(this)
            obj = this.testObj;
            task = obj.tasks{1};
            sub = obj.subjects{1};
            BOLD = obj.task_dtseries(sub, task); 
            BOLD = BOLD(5:end-5,:);

            len = length(obj.dmn_parcels);
            prec_signals = single(nan(obj.num_frames,len,obj.num_sub,obj.num_tasks));
            for pidx = 1:len
                msk = obj.mask_fs(obj.subjects{1}, obj.dmn_parcels{pidx});
                prec_signals(:,pidx,1,1) = mean(BOLD(:, msk), 2, 'omitnan');
                figure; plot(prec_signals(:,pidx)); title(obj.dmn_parcels{pidx});
                obj.write_cifti(msk, sprintf('msk_%s', obj.dmn_parcels{pidx}));
            end     
            figure; plot(obj.plvs(:,1,1)); title('plvs');
        end
        function test_call(this)
            obj = this.testObj;
            call(obj);
            mlbash(strcat('wb_view ', this.rsfMRI_scene))
        end
        function test_call_iFV_multi(this)
            pwd1 = pushd(fullfile(getenv('HCP_HOME'), '995174'));
            lp_thr = [0.05 0.1 0.2 ];
            hp_thr = [0.01 0.005 ];
            %lp_thr = 0.5;
            %hp_thr = 1/512;
            parfor idx = 1:length(lp_thr)
                fold = sprintf('arousal-waves-%g-%g', ...
                    round(hp_thr(idx),2,'significant'), round(lp_thr(idx),2,'significant'));
                fold = strrep(fold, '.', 'p');
                if ~isfolder(fold)
                    mkdir(fold);
                    pwd0 = pushd(fold);
                    fprintf('mlraut.Physio.call():  working in %s\n', pwd)
                    obj = mlraut.Physio('hp_thresh', hp_thr(idx), 'lp_thresh', lp_thr(idx), ...
                        'source_physio', '', ...
                        'gs_subtract', true);
                    call(obj)
                    popd(pwd0);
                end
            end
            mlbash(strcat('wb_view ', this.rsfMRI_scene))
            popd(pwd1)
        end
        function test_mask_fs(this)
            sub = '100307';
            parc = 'L_S_parieto_occipital';
            m = this.testObj.mask_fs(sub, parc);
            this.testObj.write_cifti(single(m), sprintf('msk_%s', this.testObj.dmn_parcels{1}));
        end
        function test_PhysioTool(this)
            root_dir_ = this.testObj.root_dir;
            sub_ = this.testObj.subjects{1};
            task_ = this.testObj.tasks{2};
            bold = fullfile(root_dir_, sub_, 'MNINonLinear', 'Results', task_, ...
                sprintf('%s_hp2000_clean.nii.gz', task_));
            wmparc = fullfile(root_dir_, sub_, 'MNINonLinear', 'ROIs', ...
                'wmparc.2.nii.gz');
            fd = fullfile(root_dir_, sub_, 'MNINonLinear', 'Results', task_, ...
                'Movement_RelativeRMS.txt');
            physio = fullfile(root_dir_, sub_, 'MNINonLinear', 'Results', task_, ...
                sprintf('%s_Physio_log.txt', task_));
            this.assertTrue(isfile(bold));
            this.assertTrue(isfile(wmparc));

            tool = mlfourd.PhysioTool(bold, wmparc, ...
                frame_displacement=fd, physio=physio);
            disp(tool)
            
            figure; plot(1:1200, tool.resp_belt); title('resp_belt');
            figure; plot(1:1200, tool.rv); title('rv');
            figure; plot(1:1200, tool.pulse_oximeter); title('pulse_oximeter');
            figure; plot(1:1200, tool.peak_amplitude); title('peak_amplitude');
            figure; plot(1:1200, tool.heart_rate); title('heart_rate');
            figure; plot(1:1200, tool.hrv); title('hrv');
        end
    end
    
    methods (TestClassSetup)
        function setupPhysio(this)
            import mlraut.*
            this.testObj_ = Physio();
            assert(isfolder(getenv('HCP_HOME')))
            this.rsfMRI_scene = fullfile(getenv('HCP_HOME'), '995174', '995174_angles_m2pi.scene');
        end
    end
    
    methods (TestMethodSetup)
        function setupPhysioTest(this)
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
