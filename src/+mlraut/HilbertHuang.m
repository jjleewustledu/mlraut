classdef HilbertHuang < handle
    %% line1
    %  line2
    %  
    %  Created 04-Sep-2022 21:25:59 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 9.12.0.2039608 (R2022a) Update 5 for MACI64.  Copyright 2022 John J. Lee.
    
    methods
        function this = call(this)
            
            ld = load('mlraut_Physio_call.mat');

            for rsn = 1:7
    
                plvs_xt_1_1 = ld.this.plvs_xt(:,:,1,1);
                plvs_rsn = plvs_xt_1_1(:, rsn == ld.this.networks_HCP);
                plvs_rsn = mean(plvs_rsn, 2);
    
                [imf_rsn,~,info_rsn] = emd(real(plvs_rsn), 'MaxNumIMF', 8, 'Display', 1);
                disp(info_rsn)
    
                figure; hht(imf_rsn(:,1), 1/0.72)
                figure; hht(imf_rsn(:,2), 1/0.72)
                figure; hht(imf_rsn(:,3), 1/0.72)
                figure; hht(imf_rsn(:,4), 1/0.72, 'FrequencyLimits', [0 0.05])
                figure; hht(imf_rsn(:,5), 1/0.72, 'FrequencyLimits', [0 0.05])
                figure; hht(imf_rsn(:,6), 1/0.72, 'FrequencyLimits', [0 0.05])
                figure; hht(imf_rsn(:,7), 1/0.72, 'FrequencyLimits', [0 0.05])
                figure; hht(imf_rsn, 1/0.72)

            end

        end
        function this = HilbertHuang(varargin)
            %% HILBERTHUANG 
            %  Args:
            %      arg1 (its_class): Description of arg1.
            
            ip = inputParser;
            addParameter(ip, "arg1", [], @(x) true)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
