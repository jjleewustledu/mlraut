classdef Test_Ciftify < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 20-Jun-2023 14:03:05 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/test/+mlraut_unittest.
    %  Developed on Matlab 9.14.0.2286388 (R2023a) Update 3 for MACI64.  Copyright 2023 John J. Lee.
    
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
    end
    
    methods (TestClassSetup)
        function setupCiftify(this)
            import mlraut.*
            this.testObj_ = Ciftify();
        end
    end
    
    methods (TestMethodSetup)
        function setupCiftifyTest(this)
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
