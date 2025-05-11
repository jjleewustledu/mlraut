classdef LateralVentricles < handle & mlraut.PhysioData
    %% Supports Gonzalez-Castillo, J., Fernandez, I. S., Handwerker, D. A. & Bandettini, P. A. 
    %  Ultra-slow fMRI fluctuations in the fourth ventricle as a marker of drowsiness. 
    %  NeuroImage 259, 119424 (2022).
    %  
    %  Created 08-May-2025 14:24:03 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 24.2.0.2923080 (R2024b) Update 6 for MACA64.  Copyright 2025 John J. Lee.
    
    properties (Dependent)
        ifv_mask
        roi_mask
    end

    methods %% GET/SET
        function g = get.ifv_mask(this)

            if ~isempty(this.ifv_mask_)
                g = this.ifv_mask_;
                return
            end

            this.ifv_mask_ = ...
                this.wmparc.numeq(4) | ...
                this.wmparc.numeq(5) | ...
                this.wmparc.numeq(43) | ...
                this.wmparc.numeq(44);  % [0 1]; FS indices 4, 5, 43, 44 are lateral ventricles
            g = this.ifv_mask_;
        end

        function g = get.roi_mask(this)
            g = this.ifv_mask;
        end
    end

    methods
        function bold = call(this)
            bold = call(this.physio_roi_);
        end
        function view_qc(this)
            this.ifv_mask.view_qc(this.ihcp_.task_signal_reference)
        end

        function this = LateralVentricles(ihcp, bold, opts)
            %% ROIs for arousal signals
            %  arguments:
            %      ihcp mlraut.HCP {mustBeNonempty}
            %      bold mlfourd.ImagingContext2
            %      opts.flipLR logical = false
            %  returns:
            %      this
                 
            arguments
                ihcp mlraut.HCP {mustBeNonempty}
                bold mlfourd.ImagingContext2
                opts.flipLR logical = false
            end
            this = this@mlraut.PhysioData(ihcp, bold);
            this.physio_roi_ = mlraut.PhysioRoi(ihcp, bold, ...
                flipLR=opts.flipLR, ...
                from_imaging_context=this.ifv_mask);
        end
    end

    %% PRIVATE

    properties (Access = private)
        ifv_mask_
        physio_roi_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
