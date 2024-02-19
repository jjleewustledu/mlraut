classdef IFourthVentricle < handle & mlraut.PhysioData
    %% Supports Gonzalez-Castillo, J., Fernandez, I. S., Handwerker, D. A. & Bandettini, P. A. 
    %  Ultra-slow fMRI fluctuations in the fourth ventricle as a marker of drowsiness. 
    %  NeuroImage 259, 119424 (2022).
    %  
    %  Created 05-Oct-2022 14:21:11 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 9.13.0.2049777 (R2022b) for MACI64.  Copyright 2022 John J. Lee.
    
    properties (Dependent)
        ifv_mask
    end

    methods %% GET/SET
        function g = get.ifv_mask(this)
            % if ~isempty(this.ifv_mask_)
            %     g = this.ifv_mask_;
            % end

            ic = this.wmparc.numeq(15); % [0 1]; FS index 15 is 4th ventricle
            idx_girth = this.find_idx_girth(ic);
            ifc = ic.imagingFormat;
            ifc.img(:,:,idx_girth+2:end) = 0;
            if this.is_7T
                ifc.fileprefix = 'ifv.1.60';  %% mm of voxels
            else
                ifc.fileprefix = 'ifv.2';  %% mm of voxels
            end
            this.ifv_mask_ = mlfourd.ImagingContext2(ifc);
            g = this.ifv_mask_;
        end
    end

    methods
        function bold = call(this)
            bold = call(this.physio_roi_);
        end
        function view_qc(this)
            this.ifv_mask.view_qc(this.ihcp_.task_signal_reference)
        end

        function this = IFourthVentricle(ihcp, bold, opts)
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

    methods (Access = private)
        function idx = find_idx_girth(~, ic)
            ic = max(ic, [], 1); 
            img_yz = squeeze(ic.imagingFormat.img);
            img_z = sum(img_yz, 1);
            [~,idx_] = max(img_z);  % idx_ of widest part of mip of 4th ventricle

            [~,idxL] = max(img_z > 0);  % idxL of img_z, a histogram
            [~,idxR_] = max(flip(img_z) > 0);
            idxR = length(img_z) - idxR_;  % idxR of img_z, a histogram
            idx_median = ceil(1 + idxL + (idxR - idxL)/2);  % idx_median of 4th ventricle along z

            idx = max(idx_, idx_median);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
