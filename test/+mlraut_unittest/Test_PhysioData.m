classdef Test_PhysioData < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 09-Feb-2024 01:09:10 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/test/+mlraut_unittest.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    
    properties
        hcp
        testObj
    end
    
    methods (Test)
        function test_afun(this)
            import mlraut.*
            this.assumeEqual(1,1);
            this.verifyEqual(1,1);
            this.assertEqual(1,1);
        end
        function test_iFV(this)
        end
        function test_gray(this)
        end
        function test_white(this)
        end
        function test_csf(this)
        end
        function test_thalamus(this)
            indices = [10, 49];
            physroi = mlraut.PhysioRoi( ...
                this.hcp, from_wmparc_indices=indices);
            %physroi.view_qc();
            this.verifyEqual(crc32_adler(0, double(physroi.roi_mask)), uint32(1999844326))
            bold = call(physroi);  % toc ~ 19 s
            %plot(bold)
            this.verifyEqual(crc32_adler(0, bold), uint32(1375506269))
        end

        function test_iFV_vs_RV(this)
        end
        function test_iFV_vs_HRV(this)
        end

        function test_flipLR(this)
        end
    end
    
    methods (TestClassSetup)
        function setupPhysioData(this)
            import mlraut.*
            this.hcp = mlraut.HCP( ...
                subjects={'995174'}, tasks={'rfMRI_REST1_7T_PA'});
            %this.testObj_ = PhysioData();
        end
    end
    
    methods (TestMethodSetup)
        function setupPhysioDataTest(this)
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
