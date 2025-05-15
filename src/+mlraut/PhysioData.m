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
        v_physio_is_inf
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
        function g = get.v_physio_is_inf(this)
            g = this.ihcp_.v_physio_is_inf;
        end
        function g = get.wmparc(this)
            if ~isempty(this.wmparc_)
                g = this.wmparc_;
            end

            this.wmparc_ = mlfourd.ImagingContext2(this.ihcp_.wmparc_fqfn);
            g = this.wmparc_;
        end
    end

    methods
        function physio = build_physio_from_ROI(this)
            %% Constructs physio vector|matrix for implementations of PhysioData 
            %  (PhysioROI, IFourthVentricle, LateralVentricles, ...)
            %
            %  Returns:
            %      physio double ~ Nt x 1 if this.v_physio_is_inf else Nt x Ngo

            physio_vec_ = this.call();
            physio_vec_gsr_ = ...
                this.ihcp_.build_global_signal_regressed(physio_vec_);
            physio_vec = ...
                this.ihcp_.build_rescaled( ...
                this.ihcp_.build_band_passed( ...
                this.ihcp_.build_centered(physio_vec_gsr_)));

            if this.v_physio_is_inf
                physio = physio_vec;
            else
                size_dtseries = size(this.ihcp_.task_dtseries());
                proi_pos = this.ihcp_.twistors.center_of_mass_position(this.roi_mask);
                physio = this.ihcp_.twistors.propagate_physio( ...
                    physio_vec, ...
                    size_bold_signal=size_dtseries, ...
                    physio_pos=proi_pos);
                assert(all(isfinite(physio), "all"), "likely that Twistors.propagate_physio is faulty")
            end
        end

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

    methods (Access = protected)
        function that = copyElement(this)
            that = copyElement@matlab.mixin.Copyable(this);
            if ~isempty(this.bold_) && ishandle(this.bold_)
                that.bold_ = copy(this.bold_); end
            if ~isempty(this.wmparc_) && ishandle(this.wmparc_)
                that.wmparc_ = copy(this.wmparc_); end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end



function mustBeValidBold(b)
    assert(isnumeric(b) || isa(b, "mlfourd.ImagingContext2"))
end
