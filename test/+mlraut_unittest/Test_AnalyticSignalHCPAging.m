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
        function test_ctor(this)            
            cd('/Volumes/PrecunealSSD2/HCPAging/HCPAgingRec/fmriresults01');
            as = mlraut.AnalyticSignalHCPAging(subjects={'HCA9992517_V1_MR'}, tasks={'fMRI_CONCAT_ALL'});
            disp(as)
        end
        function test_call(this)
            tic
            % cd('/Volumes/PrecunealSSD2/HCPAging/HCPAgingRec/fmriresults01');
            as = mlraut.AnalyticSignalHCPAging( ...
                subjects={'HCA9992517_V1_MR'}, ...
                tasks={'fMRI_CONCAT_ALL'}, ...
                do_save=true, ...
                do_save_ciftis=true, ...
                tags=stackstr(use_dashes=true));
            call(as)
            toc
            disp(as)

            % Elapsed time is 273 seconds on twistor.  Peak memory is 64 GB.  Files saved < 8 GB.
            % Elapsed time is ~791 seconds on vglab2.
        end
    end
    
    methods (TestClassSetup)
        function setupAnalyticSignalHCPAging(this)
        end
    end
    
    methods (TestMethodSetup)
        function setupAnalyticSignalHCPAgingTest(this)
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
