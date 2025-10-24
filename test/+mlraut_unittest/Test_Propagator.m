classdef Test_Propagator < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 19-Jul-2025 18:16:42 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/test/+mlraut_unittest.
    %  Developed on Matlab 25.1.0.2943329 (R2025a) for MACA64.  Copyright 2025 John J. Lee.
    
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
        function setupPropagator(this)
            import mlraut.*
            this.testObj_ = Propagator();
        end
    end
    
    methods (TestMethodSetup)
        function setupPropagatorTest(this)
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
