classdef PhysioHRV < handle & mlraut.PhysioData
    %% line1
    %  line2
    %  
    %  Created 09-Feb-2024 00:52:13 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    
    methods
        function physio = call(this)
            %% pulse ox samples <- *_Physio_log.txt, data_(:,3)
            %  Returns:
            %      physio (numeric): from <task>_Physio_log.txt

            data_ = this.physio_log();
            time_vec_bold = this.ihcp_.tr*(1:size(this.bold_,4))';
            time_vec_phys = (0:length(data_)-1)'/this.physFs;
            physio = zeros(size(time_vec_bold));            
            data_(:,3) = zscore(data_(:,3)); % of pulse ox tracing

            for idx = 1:length(physio)
                
                % For HRV, get 6 sec windows
                [~,phys_start] = min(abs(time_vec_phys - (time_vec_bold(idx) - 3)));
                [~,phys_end] = min(abs(time_vec_phys - (time_vec_bold(idx) + 3)));
               
                % For HRV, use native to avoid conflicts with chronux
                path_0 = path;
                rmpath(fullfile(getenv("HOME"), "MATLAB-Drive", "chronux_2_12", "spectral_analysis", "continuous"))
                [pks,locs] = findpeaks(data_(phys_start:phys_end,3), ...
                    'MinPeakDistance', round(this.physFs/(180/60)));
                             %,'minpeakwidth',400/(1/(200/60))); 
                             % max heart rate = 180 bpm; at 400 Hz, minimum of 100 samples apart
                locs = locs(pks > prctile(data_(phys_start:phys_end,3), 60));
                physio(idx) = mean(diff(locs), 'omitnan') / this.physFs;
                path(path_0)
            end
        end

        function this = PhysioHRV(varargin)
            this = this@mlraut.PhysioData(varargin{:});
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
