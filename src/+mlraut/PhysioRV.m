classdef PhysioRV < handle & mlraut.PhysioData
    %% line1
    %  line2
    %  
    %  Created 09-Feb-2024 00:52:06 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    
    
    methods
        function physio = call(this)
            %% resp. belt samples <- *_Physio_log.txt, data_(:,2)
            %  Returns:
            %      physio (numeric): from <task>_Physio_log.txt

            data_ = this.physio_log();
            time_vec_bold = this.ihcp_.tr*(1:size(this.bold_,4))';
            time_vec_phys = (0:size(data_, 1)-1)'/this.physFs;
            physio = zeros(size(time_vec_bold));       
            data_(:,2) = zscore(data_(:,2)); % of pulse ox tracing

            for idx = 1:length(physio)                
                % For RV, get 6 sec windows
                [~,phys_start] = min(abs(time_vec_phys - (time_vec_bold(idx) - 3)));
                [~,phys_end] = min(abs(time_vec_phys - (time_vec_bold(idx) + 3)));
                physio(idx) = this.options.statistic(data_(phys_start:phys_end,2)); % of resp. belt tracing
            end
        end
        function this = PhysioRV(varargin)
            this = this@mlraut.PhysioData(varargin{:});
            if isempty(this.options)
                this.options.statistic = @std;  % Ryan's measure in Sci. Adv. 2021
            end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
