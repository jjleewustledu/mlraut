classdef BOLDData < handle & mlsystem.IHandle
    %% line1
    %  line2
    %  
    %  Created 08-Feb-2024 15:20:10 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    
    properties (Dependent)
    end

    methods
        function this = BOLDData(ihcp)
            arguments
                ihcp mlraut.HCP
            end

            this.ihcp_ = ihcp;
        end

        function mat = task_dtseries(this, varargin)
            if ~isempty(this.task_dtseries_)
                mat = this.task_dtseries_;
                return
            end

            try
                cifti = cifti_read(this.ihcp_.task_dtseries_fqfn);
                this.ihcp_.cifti_last = cifti;
                mat = cifti.cdata';
                this.ihcp_.num_frames_ori = size(mat, 1);
                mat = this.trim_frames(mat);
            catch ME
                disp("%s: %s not found.", stackstr(), fqfn)
                disp(ME)
                mat = [];
            end
        end
        function ic = task_niigz(this)
            if ~isempty(this.task_niigz_)
                ic = this.task_niigz_;
                return
            end

            ifc = mlfourd.ImagingFormatContext2(this.ihcp_.task_niigz_fqfn);  % HCP Young Adult
            s = struct("timesMid", ascol(0:size(ifc.img, 4)));
            ifc.addJsonMetadata(s);
            ic = mlfourd.ImagingContext2(ifc);
            this.task_niigz_ = copy(ic);
        end
        function ic = task_signal_reference(this)
            if ~isempty(this.task_signal_reference_)
                ic = copy(this.task_signal_reference_);
                return
            end

            ic = mlfourd.ImagingContext2(this.ihcp_.task_signal_reference_fqfn);
            if 4 == ndims(ic)
                ic = ic.timeAveraged();
                ic.save();
            end
            this.task_signal_reference_ = copy(ic);
        end
    end

    %% PROTECTED

    properties (Access = protected)
        ihcp_
        task_dtseries_
        task_niigz_
        task_signal_reference_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
