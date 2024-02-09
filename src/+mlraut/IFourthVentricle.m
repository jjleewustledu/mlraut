classdef IFourthVentricle < handle & mlraut.PhysioData
    %% Supports Gonzalez-Castillo, J., Fernandez, I. S., Handwerker, D. A. & Bandettini, P. A. 
    %  Ultra-slow fMRI fluctuations in the fourth ventricle as a marker of drowsiness. 
    %  NeuroImage 259, 119424 (2022).
    %  
    %  Created 05-Oct-2022 14:21:11 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 9.13.0.2049777 (R2022b) for MACI64.  Copyright 2022 John J. Lee.
    
    properties (Dependent)
        aparc_a2009s
        ifv_mask
        is_7T
        subject
        wmparc
    end

    methods %% GET/SET
        function g = get.aparc_a2009s(this)
            g = this.wmparc;
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
        function g = get.wmparc(this)
            if ~isempty(this.aparc_a2009s_)
                g = this.aparc_a2009s_;
            end

            if this.is_7T
                fqfn = fullfile(this.ihcp_.root_dir, this.subject, 'MNINonLinear', 'ROIs', 'wmparc.1.60.nii.gz');
            else
                fqfn = fullfile(this.ihcp_.root_dir, this.subject, 'MNINonLinear', 'ROIs', 'wmparc.2.nii.gz');
                if ~isfile(fqfn)
                    fqfn = fullfile(this.ihcp_.root_dir, this.subject, 'MNINonLinear', 'wmparc.nii.gz');
                end
            end
            this.aparc_a2009s_ = mlfourd.ImagingContext2(fqfn);
            g = this.aparc_a2009s_;
        end
    end

    methods
        function bold = call(this)
            fMRI = this.task_niigz();
            ic = fMRI.volumeAveraged(this.ifv_mask);
            bold = ascol(ic.nifti.img);
        end
        function ic = task_niigz(this)
            ic = this.ihcp_.task_niigz();
        end
        function view_qc(this)
            this.ifv_mask.view_qc(this.ihcp_.task_signal_reference)
        end

        function this = IFourthVentricle(ihcp, subject, task)
            %% IFOURTHVENTRICLE 
            %  Args:
            %      ihcp mlraut.HCP : client possessing HCP information, esp. filesystem information.
            %      subject {mustBeTextScalar} : e.g., 995174.
            %      task mustBeTextScalar : e.g., rfMRI_REST1_LR.
            
            arguments
                ihcp mlraut.HCP {mustBeNonempty}
                subject {mustBeTextScalar}
                task {mustBeTextScalar}
            end
            this = this@mlraut.PhysioData(ihcp);
            this.subject_ = subject;
            this.task_ = task;
        end
    end

    %% PRIVATE

    properties (Access = private)
        aparc_a2009s_
        task_niigz_
        subject_
        task_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
