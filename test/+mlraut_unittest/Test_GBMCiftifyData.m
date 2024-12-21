classdef Test_GBMCiftifyData < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 20-Dec-2024 16:11:50 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/test/+mlraut_unittest.
    %  Developed on Matlab 24.2.0.2806996 (R2024b) Update 3 for MACA64.  Copyright 2024 John J. Lee.
    
    properties
        testObj
        asGBM
    end
    
    methods (Test)
        function test_afun(this)
            import mlraut.*
            this.assumeEqual(1,1);
            this.verifyEqual(1,1);
            this.assertEqual(1,1);
        end

        function test_copy(this)
            disp(this.testObj)
        end

        function test_build_gbm_json(this)

            fqfn = this.testObj.json_fqfn;
            deleteExisting(fqfn);
            this.testObj.build_gbm_json();

            % Step 1: Read the file content
            jsonText = fileread(fqfn);

            % Step 2: Decode the JSON content into MATLAB data
            jsonData = jsondecode(jsonText);

            % Step 3: Verify format (optional)
            if ~isstruct(jsonData) && ~iscell(jsonData) && ~isnumeric(jsonData)
                error('mlraut:UnittestError', '%s:  Invalid JSON format.', stackstr());
            end

            % Step 4: Print for verification
            disp('gbm.json:');
            disp(jsonData);  % For human-readable output
        end

        function test_build_CE_WT_edema(this)
        end

        function test_table_gbm(this)            
            excluded = string(this.testObj.table_excluded.I3CRID);
            t = this.testObj.table_gbm;
            this.verifyFalse(sum(contains(t.I3CRID, excluded)) > 0)
        end
    end
    
    methods (TestClassSetup)
        function setupGBMCiftifyData(this)
            SUB = {'sub-I3CR1156'};  % OS ~ 18 days, 74 yo
            v_physio_iFV = 50;
            root_dir = '/Volumes/PrecunealSSD2/AnalyticSignalGBM/analytic_signal/dockerout/ciftify';
            cd(root_dir);
            out_dir = '/Volumes/PrecunealSSD2/AnalyticSignalGBM/analytic_signal/matlabout';
            ensuredir(out_dir);

            this.asGBM = mlraut.AnalyticSignalGBM( ...
                subjects=SUB, ...
                tasks={'ses-1_task-rest_run-01_desc-preproc'}, ...
                do_resting=true, ...
                do_save=false, ...
                do_save_dynamic=false, ...
                do_save_ciftis=false, ...
                out_dir=out_dir, ...
                v_physio=v_physio_iFV, ...
                plot_range=1:69, ...
                do_plot_networks=true, ...
                source_physio='iFV');
            this.testObj_ = mlraut.GBMCiftifyData(this.asGBM);
        end
    end
    
    methods (TestMethodSetup)
        function setupGBMCiftifyDataTest(this)
            this.testObj = copy(this.testObj_);
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
