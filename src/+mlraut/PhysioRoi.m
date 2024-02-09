classdef PhysioRoi < handle & mlraut.PhysioData
    %% Extends the concepts underlying mlraut.IFourthVentricle to arbitrary ROIs,
    %  especially pathophysiological ROIs such as samples from within tumors.
    %  
    %  Created 10-Aug-2023 18:11:35 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 9.14.0.2306882 (R2023a) Update 4 for MACI64.  Copyright 2023 John J. Lee.
    
    properties (Dependent)
        is_7T
        roi_mask
        subject
    end

    methods %% GET/SET
        function g = get.roi_mask(this)
            if ~isempty(this.roi_mask_)
                g = this.roi_mask_;
            end

            fqfn = fullfile(this.ihcp_.root_dir, this.subject, 'MNINonLinear', 'ROIs', [this.roi_fileprefix_, '.nii.gz']);
            if ~isfile(fqfn)
                fqfn = fullfile(this.ihcp_.root_dir, this.subject, 'MNINonLinear',  [this.roi_fileprefix_, '.nii.gz']);
            end
            this.roi_mask_ = mlfourd.ImagingContext2(fqfn);
            SBRef = this.ihcp_.task_signal_reference();
            mask2 = SBRef.blurred(7).thresh(100).binarized();
            if this.flipLR_
                this.roi_mask_ = flip(this.roi_mask_, 1);
                mask2 = flip(mask2, 1);
            end
            this.roi_mask_ = this.roi_mask_ .* mask2;
            g = this.roi_mask_;
        end
        function g = get.is_7T(this)
            g = contains(this.task_, '7T');
        end
        function g = get.subject(this)
            g = this.subject_;
        end
    end
    
    methods
        function bold = call(this)
            fMRI = this.task_niigz();
            ic = fMRI.volumeAveraged(this.roi_mask);
            bold = ascol(ic.nifti.img);
        end
        function ic = task_niigz(this)
            ic = this.ihcp_.task_niigz();
            if this.flipLR_
                ic = flip(ic, 1);
            end
        end
        function view_qc(this)
            this.roi_mask.view_qc(flip(this.ihcp_.task_signal_reference, 1))
        end

        function this = PhysioRoi(ihcp, subject, task, fileprefix, opts)
            %% PHYSIOROI 
            %  Args:
            %      ihcp mlraut.HCP : client possessing HCP information, esp. filesystem information.
            %      subject {mustBeTextScalar} : e.g., 995174.
            %      task {mustBeTextScalar} : e.g., rfMRI_REST1_LR.
            %      fileprefix }mustBeTextScalar} : e.g., WT_on_T1w.
            %      opts.flipLR logical = false : flipping may be needed to match task_dtseries, task_niigz
            
            arguments
                ihcp mlraut.HCP {mustBeNonempty}
                subject {mustBeTextScalar}
                task {mustBeTextScalar}
                fileprefix {mustBeTextScalar}
                opts.flipLR logical = false
            end
            this = this@mlraut.PhysioData(ihcp);
            this.subject_ = subject;
            this.task_ = task;
            this.roi_fileprefix_ = fileprefix;
            this.flipLR_ = opts.flipLR;
        end
    end

    %% PRIVATE

    properties (Access = private)
        flipLR_
        roi_fileprefix_
        roi_mask_
        subject_
        task_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
