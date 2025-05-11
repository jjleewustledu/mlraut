classdef IFourthVentricle < handle & mlraut.PhysioData
    %% Supports Gonzalez-Castillo, J., Fernandez, I. S., Handwerker, D. A. & Bandettini, P. A. 
    %  Ultra-slow fMRI fluctuations in the fourth ventricle as a marker of drowsiness. 
    %  NeuroImage 259, 119424 (2022).
    %  
    %  Created 05-Oct-2022 14:21:11 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 9.13.0.2049777 (R2022b) for MACI64.  Copyright 2022 John J. Lee.
    
    properties
        prefer_brightest  % if there exists file MNINonLinear/brightest.nii.gz, manually created, preferentially use it.
    end

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

            ic = this.wmparc.numeq(15); % [0 1]; FS index 15 is 4th ventricle
            ifc = ic.imagingFormat;
            switch char(this.selected_voxels_)
                case 'brightest'
                    newifc = this.build_brightest_voxels(ifc);
                    newifc.fileprefix = strrep(newifc.fileprefix, 'ifv', 'ifvbrightest');
                case 'inferior-quantile'
                    brightest = this.build_brightest_voxels(ifc);
                    brightest_msk = logical(brightest.img);

                    idx_girth = this.find_idx_girth(ic);
                    inferior_msk = logical(ifc.img);
                    inferior_msk(:,:,idx_girth+1:end) = false;  % remove superior half

                    newifc = copy(ifc);
                    newifc.img = single(inferior_msk & ~brightest_msk);
                    newifc.fileprefix = strrep(newifc.fileprefix, 'ifv', 'ifviq');
                case 'inferior-half'
                    idx_girth = this.find_idx_girth(ic);
                    newifc = copy(ifc);
                    newifc.img(:,:,idx_girth+1:end) = 0;  % remove superior half
                case 'superior-half'
                    idx_girth = this.find_idx_girth(ic);
                    newifc = copy(ifc);
                    newifc.img(:,:,1:idx_girth) = 0;  % remove inferior half
                    newifc.fileprefix = strrep(newifc.fileprefix, 'ifv', 'sfv');
                case 'all'
                    newifc = copy(ifc);
                    newifc.fileprefix = strrep(newifc.fileprefix, 'ifv', 'fv');
                otherwise
                    error("mlraut:ValueError", stackstr())
            end
            if this.is_7T
                newifc.fileprefix = 'ifv.1.60';  %% mm of voxels
            else
                newifc.fileprefix = 'ifv.2';  %% mm of voxels
            end
            this.ifv_mask_ = mlfourd.ImagingContext2(newifc);
            g = this.ifv_mask_;
        end
        function g = get.roi_mask(this)
            g = this.ifv_mask;
        end
    end

    methods
        function newmask_ifc = build_brightest_voxels(this, mask_ifc)
            %% select best voxels in mask corresponding to inferior-most 4 slices of mask_ifc;
            %  rationale follows Fultz, et al.  https://www.science.org/doi/10.1126/science.aax5440.

            b_fqfn = fullfile(this.ihcp_.cohort_data.mninonlinear_dir, "brightest.nii.gz");
            if this.prefer_brightest && isfile(b_fqfn)
                b1_fqfn = strrep(b_fqfn, "brightest", "brightest_nifti1");
                mlbash(sprintf("fslchfiletype NIFTI1_GZ %s %s", b_fqfn, b1_fqfn));
                newmask_ifc = mlfourd.ImagingFormatContext2(b1_fqfn);
                if dipmax(newmask_ifc.img) > 1
                    newmask_ifc.img = single(newmask_ifc.img > 0);  % precautionary
                end
                assert(dipmax(newmask_ifc.img) == 1)
                return
            end
            
            % Fultz' 4 slices spanned 10 mm.
            Nslices = ceil(10 / mask_ifc.mmppix(3));
            idx_inferior = this.find_idx_inferior(mlfourd.ImagingContext2(mask_ifc));
            indices = (0:(Nslices-1)) + idx_inferior;
            
            newmask = zeros(size(mask_ifc.img));
            newmask(:,:,indices) = mask_ifc.img(:,:,indices);

            newmask_ifc = copy(mask_ifc);
            newmask_ifc.img = newmask;
            newmask_ifc.fileprefix = mask_ifc.fileprefix + "_" + stackstr();
        end

        function bold = call(this)
            bold = call(this.physio_roi_);
        end

        function view_qc(this, ref)
            arguments
                this mlraut.IFourthVentricle
                ref = this.ihcp_.task_signal_reference
            end

            this.ifv_mask.view_qc(ref)
        end

        function this = IFourthVentricle(ihcp, bold, opts)
            %% ROIs for arousal signals
            %  arguments:
            %      ihcp mlraut.HCP {mustBeNonempty}
            %      bold mlfourd.ImagingContext2
            %      opts.prefer_brightest logical = true
            %      opts.best_voxels {mustBeTextScalar} = "brightest"  % "inferior-quantile", "inferior-half", "superior-half", "all"
            %      opts.flipLR logical = false
            %  returns:
            %      this
                 
            arguments
                ihcp mlraut.HCP {mustBeNonempty}
                bold mlfourd.ImagingContext2
                opts.prefer_brightest logical = true
                opts.selected_voxels {mustBeTextScalar} = "brightest"  % "inferior-quantile", "inferior-half", "superior-half", "all"
                opts.flipLR logical = false
            end
            this = this@mlraut.PhysioData(ihcp, bold);

            this.prefer_brightest = opts.prefer_brightest;
            this.selected_voxels_ = opts.selected_voxels;
            this.physio_roi_ = mlraut.PhysioRoi(ihcp, bold, ...
                flipLR=opts.flipLR, ...
                from_imaging_context=this.ifv_mask);
        end
    end

    %% PRIVATE

    properties (Access = private)
        ifv_mask_
        physio_roi_    
        selected_voxels_
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

        function idx = find_idx_inferior(~, ic)
            ic = max(max(ic, [], 1), [], 2);  % mip projected to z
            img_z = squeeze(ic.imagingFormat.img);
            idx = find(img_z, 1, 'first');
            assert(~isempty(idx))
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
