classdef PhysioData < handle & mlsystem.IHandle
    %% line1
    %  line2
    %  
    %  Created 08-Feb-2024 19:45:15 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    

    properties
        min_physN = 860  % min physio samples to accept
        physFs = 400  % Physio sampling rate, Hz
    end

    properties (Dependent)
        is_7T
        sub
        task
        wmparc
    end

    methods %% GET
        function g = get.is_7T(this)
            g = contains(this.task, '7T');
        end
        function g = get.sub(this)
            g = this.ihcp_.current_subject;
        end
        function g = get.task(this)
            g = this.ihcp_.current_task;
        end
        function g = get.wmparc(this)
            % if ~isempty(this.wmparc_)
            %     g = this.wmparc_;
            % end

            this.wmparc_ = mlfourd.ImagingContext2(this.ihcp_.wmparc_fqfn);
            g = this.wmparc_;
        end
    end

    methods
        function data = physio_log(this)
            fqfn = fullfile( ...
                this.ihcp_.task_dir, strcat(this.task, '_Physio_log.txt'));
            data = importdata(fqfn);
            assert(length(data)/this.physFs >= this.min_physN, stackstr(2))
        end

        function this = PhysioData(ihcp, bold)
            arguments
                ihcp mlraut.HCP {mustBeNonempty}
                bold {mustBeValidBold}
            end

            this.bold_ = bold;
            this.ihcp_ = ihcp;
        end
    end

    %% PROTECTED

    properties (Access = protected)
        bold_
        ihcp_
        wmparc_    
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end



function mustBeValidBold(b)
    assert(isnumeric(b) || isa(b, "mlfourd.ImagingContext2"))
end
