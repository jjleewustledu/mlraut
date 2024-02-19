classdef Test_BOLDData < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 08-Feb-2024 15:20:10 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/test/+mlraut_unittest.
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
        function test_global_signal(this)
            %% check conformance with Petersen's expectations

        end
        function test_iFV(this)
            %% check conformance with results of Bandettini

            this.verifyTrue(size(bold, 2) > 1);
            plot(sum(bold, 2));
        end
        function test_PhysioRoi(this)
            %% overlap with iFV

            %% verify key structures
        end
        function test_GBM(this)
            %% visualize GBM in BOLD
        end
        function test_signal_lateralization(this)
            %% verify flipping, avoidance of double flips
        end
        function test_signal_normalizations(this)
            %% show appropriateness of chosen norm
        end
        function test_signal_histograms(this)
        end
        function test_signal_powerlaw(this)
        end
    end
    
    methods (TestClassSetup)
        function setupBOLDData(this)
            import mlraut.*
            this.testObj_ = BOLDData();
        end
    end
    
    methods (TestMethodSetup)
        function setupBOLDDataTest(this)
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
