classdef PhysioRoi < handle & mlraut.PhysioData
    %% Extends the concepts underlying mlraut.IFourthVentricle to arbitrary ROIs,
    %  especially pathophysiological ROIs such as samples from within tumors.
    %  Internally manages whether to perform flipLR.
    %  call(this) to obtain physio signal as col vector for time-series.
    %  
    %  Created 10-Aug-2023 18:11:35 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 9.14.0.2306882 (R2023a) Update 4 for MACI64.  Copyright 2023 John J. Lee.
    
    properties (Dependent)
        roi_mask
    end

    methods %% GET/SET
        function g = get.roi_mask(this)
            g = this.roi_mask_;
        end
    end
    
    methods
        function ic = build_roi_mask(this, opts)
            %% Performs exactly one flipLR as needed.
            %  Args:
            %      this mlraut.PhysioRoi
            %      opts.from_imaging_context = [] % arbitary mask
            %      opts.from_wmparc_indices double = [] % from wmparc

            arguments
                this mlraut.PhysioRoi
                opts.from_imaging_context = []
                opts.from_wmparc_indices double = []
            end

            if ~isempty(opts.from_imaging_context)
                ic = mlfourd.ImagingContext2(opts.from_imaging_context);
                SBmsk = this.ihcp_.task_signal_mask();  % minimize extra-parenchymal ROI
                ic = ic .* SBmsk;
                ic = ic.binarized();
                if this.flipLR_
                    ic = flip(ic, 1);
                end
                return
            end
            if ~isempty(opts.from_wmparc_indices)
                wmp = mlfourd.ImagingContext2(this.ihcp_.wmparc_fqfn);
                ic = zeros(wmp);                
                for idx = opts.from_wmparc_indices
                    ic = ic + wmp.numeq(idx);
                end
                ic = ic.binarized();
                ic.fileprefix = wmp.fileprefix + "-idx" + strrep(num2str(opts.from_wmparc_indices), "  ", ",");
                if this.flipLR_
                    ic = flip(ic, 1);
                end
                return
            end

            % return trivial unit mask
            wmp = mlfourd.ImagingContext2(this.ihcp_.wmparc_fqfn);
            ic = ones(wmp);
        end

        function bold = call(this)
            %% Performs exactly one flipLR as needed.

            fMRI = this.bold_;
            if this.flipLR_
                assert(isa(fMRI, "mlfourd.ImagingContext2"))
                fMRI = flip(fMRI, 1);
            end
            ic = fMRI.volumeAveraged(this.roi_mask);
            bold = ascol(ic.nifti.img);
        end

        function view_qc(this)
            %% Performs exactly one flipLR as needed.
            %  QC by comparison to task_signal_reference.

            task_signal_reference = this.ihcp_.task_signal_reference;
            if this.flipLR_
                assert(isa(task_signal_reference, "mlfourd.ImagingContext2"))
                task_signal_reference = flip(task_signal_reference, 1);
            end
            this.roi_mask.view_qc(task_signal_reference)
        end

        function this = PhysioRoi(ihcp, bold, opts)
            %% PHYSIOROI 
            %  Args:
            %      ihcp mlraut.HCP : client possessing HCP information, esp. filesystem information.
            %      opts.flipLR logical = false : flipping may be needed to match task_dtseries, task_niigz
            %      opts.from_imaging_context = [] : arbitary mask
            %      opts.from_wmparc_indices double = [] : from wmparc
            
            arguments
                ihcp mlraut.HCP {mustBeNonempty}
                bold {mustBeValidBold}
                opts.flipLR logical = false
                opts.from_imaging_context = []
                opts.from_wmparc_indices double = []
            end
            this = this@mlraut.PhysioData(ihcp, bold);
            this.flipLR_ = opts.flipLR;
            this.roi_mask_ = this.build_roi_mask( ...
                from_imaging_context=opts.from_imaging_context, ...
                from_wmparc_indices=opts.from_wmparc_indices);
        end
    end

    %% PROTECTED

    properties (Access = private)
        flipLR_
        roi_mask_
    end

    methods (Access = protected)
        function that = copyElement(this)
            that = copyElement@matlab.mixin.Copyable(this);
            if ~isempty(this.roi_mask_) && ishandle(this.roi_mask_)
                that.roi_mask_ = copy(this.roi_mask_); end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end



function mustBeValidBold(b)
    assert(isnumeric(b) || isa(b, "mlfourd.ImagingContext2"))
end
