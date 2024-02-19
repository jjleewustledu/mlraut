classdef Test_CohortData < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 08-Feb-2024 18:43:56 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/test/+mlraut_unittest.
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
        function test_HCPYoungAdult(this)
        end
        function test_HCPAging(this)
        end
        function test_GBMCiftify(this)
        end
    end
    
    methods (TestClassSetup)
        function setupCohortData(this)
            import mlraut.*
            this.testObj_ = CohortData();
        end
    end
    
    methods (TestMethodSetup)
        function setupCohortDataTest(this)
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
