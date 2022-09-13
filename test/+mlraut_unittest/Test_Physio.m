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
            mlbash(['wb_view ' this.rsfMRI_scene])
        end
        function test_mask_fs(this)
            sub = '100307';
            parc = 'L_S_parieto_occipital';
            m = this.testObj.mask_fs(sub, parc);
            this.testObj.write_cifti(single(m), sprintf('msk_%s', this.testObj.dmn_parcels{1}));
        end
    end
    
    methods (TestClassSetup)
        function setupPhysio(this)
            import mlraut.*
            this.testObj_ = Physio();
            this.rsfMRI_scene = fullfile(getenv('SINGULARITY_HOME'), '..', ...
                'HCP', 'HCP_WB_Tutorial_1.5_Pr_kN3mg', 'WB_1.5_SCENES.32k_fs_LR.scene');
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
