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

            this.task_dtseries;  % obtains this.num_frames_ori_ from cifti_read()
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
            if ~isempty(this.task_dtseries_)
                mat = this.task_dtseries_;
                return
            end   % caching requires BOLDData to be reset by HPC.malloc()

            try
                cifti = cifti_read(this.ihcp_.task_dtseries_fqfn);
                mat = cifti.cdata';
                this.task_dtseries_ = mat;
                this.num_frames_ori_ = size(mat, 1);
            catch ME
                disp("%s: error while attempting to cifti_read %s.", stackstr(), this.ihcp_.task_dtseries_fqfn)
                handexcept(ME)
            end
        end
        function ic = task_niigz(this)
            if ~isempty(this.task_niigz_)
                ic = this.task_niigz_;
                return
            end   % caching requires BOLDData to be reset by HPC.malloc()

            ifc = mlfourd.ImagingFormatContext2(this.ihcp_.task_niigz_fqfn);  % HCP Young Adult
            tr = this.ihcp_.tr;
            T = tr*(size(ifc.img, 4) - 1);
            s = struct("timesMid", ascol(0:tr:T));
            ifc.addJsonMetadata(s);
            ic = mlfourd.ImagingContext2(ifc);
            this.task_niigz_ = copy(ic);
        end
        function ic = task_signal_reference(this)
            if ~isempty(this.task_signal_reference_)
                ic = copy(this.task_signal_reference_);
                return
            end   % caching requires BOLDData to be reset by HPC.malloc()

            fqfn_avgt = strrep(this.ihcp_.task_signal_reference_fqfn, ".nii.gz", "_avgt.nii.gz");
            if isfile(fqfn_avgt)
                ic = mlfourd.ImagingContext2(fqfn_avgt);
                return
            end

            ic = mlfourd.ImagingContext2(this.ihcp_.task_signal_reference_fqfn);
            if 4 == ndims(ic)
                % very slow
                fprintf("%s:  time-averaging %s\n", stackstr(), ic.fqfn)
                fprintf("\tplease wait...\n")
                ic = ic.timeAveraged();
                ic.save();
                fprintf("\tcomplete!\n")
            end
            this.task_signal_reference_ = copy(ic);
        end
    end

    methods (Static)
        function ifc = concat_nii(fn1, fn2, fn, opts)
            arguments
                fn1 {mustBeFile}
                fn2 {mustBeFile}
                fn {mustBeTextScalar} = ""
                opts.do_save logical = true
            end
            assert(endsWith(fn1, ".nii.gz"))
            assert(endsWith(fn2, ".nii.gz"))
            if isemptytext(fn)
                fn = regexprep(fn1, 'run-\d+', 'run-all');
            end

            ifc1 = mlfourd.ImagingFormatContext2(fn1);            
            ifc2 = mlfourd.ImagingFormatContext2(fn2);
            ifc = copy(ifc1);
            ifc.img = cat(4, ifc1.img, ifc2.img);
            ifc.fileprefix = mybasename(fn);

            if opts.do_save
                ifc.save();
            end            
        end

        function c = concat_dtseries(fn1, fn2, fn, opts)
            arguments
                fn1 {mustBeFile}
                fn2 {mustBeFile}
                fn {mustBeTextScalar} = ""
                opts.do_save logical = true
            end
            assert(endsWith(fn1, ".dtseries.nii"))
            assert(endsWith(fn2, ".dtseries.nii"))
            if isemptytext(fn)
                fn = regexprep(fn1, 'run-\d+', 'run-all');
            end

            c1 = cifti_read(fn1);
            c2 = cifti_read(fn2);
            seriesStep1 = c1.diminfo{2}.seriesStep;
            seriesStep2 = c2.diminfo{2}.seriesStep;
            assert(seriesStep1 == seriesStep2)
            seriesUnit1 = c1.diminfo{2}.seriesUnit;
            seriesUnit2 = c2.diminfo{2}.seriesUnit;
            assert(seriesUnit1 == seriesUnit2)
            length1 = c1.diminfo{2}.length;
            length2 = c2.diminfo{2}.length;
            
            c = c1;
            c.cdata = [c1.cdata, c2.cdata];
            c.diminfo{2} = cifti_diminfo_make_series( ...
                length1 + length2, 0, seriesStep1, seriesUnit1);

            if opts.do_save
                cifti_write(c, convertStringsToChars(fn));
            end
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

    methods (Access = protected)
        function that = copyElement(this)
            that = copyElement@matlab.mixin.Copyable(this);
            if ~isempty(this.task_dtseries_)
                that.task_dtseries_ = copy(this.task_dtseries_); end
            if ~isempty(this.task_niigz_)
                that.task_niigz_ = copy(this.task_niigz_); end
            if ~isempty(this.task_signal_reference_)
                that.task_signal_reference_ = copy(this.task_signal_reference_); end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
