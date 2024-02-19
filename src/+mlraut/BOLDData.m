classdef BOLDData < handle & mlsystem.IHandle
    %% line1
    %  line2
    %  
    %  Created 08-Feb-2024 15:20:10 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    
    properties (Dependent)
        is_7T
        num_frames_ori % expected:  1200 for HCP, 160 for RT GBM
        num_nodes
        sub
        task
    end

    methods %% GET
        function g = get.is_7T(this)
            g = contains(this.task, '7T');
        end
        function g = get.num_frames_ori(this)
            if ~isempty(this.num_frames_ori_)
                g = this.num_frames_ori_;
                return
            end

            this.num_frames_ori_ = size(this.task_dtseries, 1);
            g = this.num_frames_ori_;
        end
        function g = get.num_nodes(this)
            if this.is_7T
                g = 170494;  % HCP standard 1.6mm "grayordinates" 
            else
                g = 91282;  % HCP standard 2mm "grayordinates" 
            end
        end
        function g = get.sub(this)
            g = this.ihcp_.current_subject;
        end
        function g = get.task(this)
            g = this.ihcp_.current_task;
        end
    end

    methods
        function this = BOLDData(ihcp)
            arguments
                ihcp mlraut.HCP
            end

            this.ihcp_ = ihcp;
        end

        function mat = task_dtseries(this, varargin)
            % if ~isempty(this.task_dtseries_)
            %     mat = this.task_dtseries_;
            %     return
            % end

            try
                cifti = cifti_read(this.ihcp_.task_dtseries_fqfn);
                this.ihcp_.cifti_last = cifti;
                mat = cifti.cdata';
                this.num_frames_ori_ = size(mat, 1);
            catch ME
                disp("%s: error while attempting to cifti_read %s.", stackstr(), this.ihcp_.task_dtseries_fqfn)
                handexcept(ME)
            end
        end
        function ic = task_niigz(this)
            % if ~isempty(this.task_niigz_)
            %     ic = this.task_niigz_;
            %     return
            % end

            ifc = mlfourd.ImagingFormatContext2(this.ihcp_.task_niigz_fqfn);  % HCP Young Adult
            tr = this.ihcp_.tr;
            T = tr*(size(ifc.img, 4) - 1);
            s = struct("timesMid", ascol(0:tr:T));
            ifc.addJsonMetadata(s);
            ic = mlfourd.ImagingContext2(ifc);
            this.task_niigz_ = copy(ic);
        end
        function ic = task_signal_reference(this)
            % if ~isempty(this.task_signal_reference_)
            %     ic = copy(this.task_signal_reference_);
            %     return
            % end

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
        num_frames_ori_
        task_dtseries_
        task_niigz_
        task_signal_reference_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
