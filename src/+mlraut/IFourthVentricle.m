classdef IFourthVentricle < handle
    %% line1
    %  line2
    %  
    %  Created 05-Oct-2022 14:21:11 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 9.13.0.2049777 (R2022b) for MACI64.  Copyright 2022 John J. Lee.
    
    properties (Dependent)
        aparc_a2009s
        fMRI
        ifv_mask
        is_7T
        subject
    end

    methods

        %% GET/SET

        function g = get.fMRI(this)
            if ~isempty(this.fMRI_)
                g = this.fMRI_;
                return
            end

            fqfn = fullfile(this.physio_.root_dir, this.subject, 'MNINonLinear', 'Results', this.task_, ...
                sprintf('%s_hp2000_clean.nii.gz', this.task_));
            this.fMRI_ = mlfourd.ImagingContext2(fqfn);
            g = this.fMRI_;
        end
        function g = get.aparc_a2009s(this)
            if ~isempty(this.aparc_a2009s_)
                g = this.aparc_a2009s_;
            end

            if this.is_7T
                fqfn = fullfile(this.physio_.root_dir, this.subject, 'MNINonLinear', 'ROIs', 'wmparc.1.60.nii.gz');
            else
                fqfn = fullfile(this.physio_.root_dir, this.subject, 'MNINonLinear', 'ROIs', 'wmparc.2.nii.gz');
            end
            this.aparc_a2009s_ = mlfourd.ImagingContext2(fqfn);
            g = this.aparc_a2009s_;
        end
        function g = get.ifv_mask(this)
            ic = this.aparc_a2009s.numeq(15); % [0 1]; 15 is 4th ventricle
            ifc = ic.nifti;
            ifc.img(:,:,23:end) = 0;
            if this.is_7T
                ifc.fileprefix = 'ifv.1.60';
            else
                ifc.fileprefix = 'ifv.2';
            end
            g = mlfourd.ImagingContext2(ifc);
        end
        function g = get.is_7T(this)
            g = contains(this.task_, '7T');
        end
        function g = get.subject(this)
            g = this.subject_;
        end

        %%

        function bold = call(this)
            ic = this.fMRI.volumeAveraged(this.ifv_mask);
            bold = ascol(ic.nifti.img);
        end

        function this = IFourthVentricle(physio, subject, task)
            %% IFOURTHVENTRICLE 
            %  Args:
            %      physio mlraut.Physio : client possessing filesystem information.
            %      task mustBeTextScalar : e.g., rfMRI_REST1_LR
            
            arguments
                physio mlraut.Physio
                subject {mustBeTextScalar}
                task {mustBeTextScalar}
            end
            this.physio_ = physio;
            this.subject_ = subject;
            this.task_ = task;
        end
    end

    %% PRIVATE

    properties (Access = private)
        aparc_a2009s_
        fMRI_
        physio_
        subject_
        task_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
