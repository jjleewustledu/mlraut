classdef Test_CiftiGsp < matlab.unittest.TestCase
	%% TEST_CIFTIGSP 

	%  Usage:  >> results = run(mlraut_unittest.Test_CiftiGsp)
 	%          >> result  = run(mlraut_unittest.Test_CiftiGsp, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 26-Mar-2021 14:26:02 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlraut/test/+mlraut_unittest.
 	%% It was developed on Matlab 9.9.0.1592791 (R2020b) Update 5 for MACI64.  Copyright 2021 John Joowon Lee.
 	
	properties
 		registry
 		testObj
 	end

	methods (Test)
		function test_afun(this)
 			import mlraut.*;
 			ciftis = this.testObj.template();
            this.assertTrue(~isempty(ciftis))
            disp(ciftis{1})
        end
        function test_writeClusters(this)
            ld = load('MSC_clusters.mat');
            this.testObj.writeClusters(ld.clusters, 'MSC_clusters')
        end
	end

 	methods (TestClassSetup)
		function setupCiftiGsp(this)
 			import mlraut.*;
 			this.testObj_ = CiftiGsp('outdir', pwd);
 		end
	end

 	methods (TestMethodSetup)
		function setupCiftiGspTest(this)
 			this.testObj = this.testObj_;
 			this.addTeardown(@this.cleanTestMethod);
 		end
	end

	properties (Access = private)
 		testObj_
 	end

	methods (Access = private)
		function cleanTestMethod(this)
 		end
	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

