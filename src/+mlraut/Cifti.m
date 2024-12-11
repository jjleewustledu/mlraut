classdef Cifti < handle & mlsystem.IHandle
    %% line1
    %  line2
    %  
    %  Created 10-Feb-2024 13:47:24 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    

    properties (Dependent)
        is_7T
        out_dir
        template_dscalar
        template_niigz
    end

    methods %% GET
        function g = get.is_7T(this)
            g = contains(this.ihcp_.task, '7T');
        end
        function g = get.out_dir(this)
            g = this.ihcp_.out_dir;
        end
        function g = get.template_dscalar(this)
            g = this.ihcp_.task_dtseries_fqfn;
        end
        function g = get.template_niigz(this)
            g = this.ihcp_.wmparc_fqfn;
        end
    end

    methods
        function cii = aparc_a2009s_dlabel_nii(this, sub)
            %  e.g.:
            %  cd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200/100307/MNINonLinear/fsaverage_LR32k')
            %  cii = cifti_read('100307.aparc.a2009s.32k_fs_LR.dlabel.nii')

            arguments
                this mlraut.Cifti
                sub {mustBeTextScalar} = this.ihcp_.current_subject
            end
            pth = fullfile(this.ihcp_.root_dir, sub, 'MNINonLinear', 'fsaverage_LR32k');
            fn = sprintf('%s.aparc.a2009s.32k_fs_LR.dlabel.nii', sub);
            fn = fullfile(pth, fn);
            cii = cifti_read(fn);
        end
        function c_ = write_cifti(this, c1_data, fn)
            cifti_last = cifti_read(this.template_dscalar);
            
            sz = size(c1_data);
            if sz(2) > sz(1)
                c1_data = c1_data'; % grey-ordinates x series
                sz = sort(sz, 'descend');
            end
            if sz(2) > 1 && ~contains(fn, '.dtseries')
                [pth,fp,ext] = myfileparts(fn);
                if isempty(pth) || "" == pth
                    pth = this.out_dir;
                end
                fn = strcat(fullfile(pth, fp), '.dtseries', ext);
            end
            if sz(2) == 1 && ~contains(fn, '.dscalar')
                [pth,fp,ext] = myfileparts(fn);
                if isempty(pth) || "" == pth
                    pth = this.out_dir;
                end
                fn = strcat(fullfile(pth, fp), '.dscalar', ext);
            end
            if ~contains(fn, '.nii')
                fn = strcat(fn, '.nii');
            end
            c_ = cifti_last;
            c_.cdata = c1_data;
            if contains(fn, '.dtseries')
                c_.diminfo{2} = cifti_diminfo_make_series(sz(2), 0, this.ihcp_.tr, 'SECOND');
            else
                c_.diminfo{2} = cifti_diminfo_make_scalars(1);
            end
            cifti_write(c_, convertStringsToChars(fn));
        end 
        function c_ = write_cifti_previous(this, c1_data, fn)
            sz = size(c1_data);
            if sz(2) > sz(1)
                c1_data = c1_data'; % grey-ordinates x series
                sz = sort(sz, 'descend');
            end
            if sz(2) > 1 && ~contains(fn, '.dtseries')
                [pth,fp,ext] = myfileparts(fn);
                if isempty(pth) || "" == pth
                    pth = this.out_dir;
                end
                fn = strcat(fullfile(pth, fp), '.dtseries', ext);
            end
            if sz(2) == 1 && ~contains(fn, '.dscalar')
                [pth,fp,ext] = myfileparts(fn);
                if isempty(pth) || "" == pth
                    pth = this.out_dir;
                end
                fn = strcat(fullfile(pth, fp), '.dscalar', ext);
            end
            if ~contains(fn, '.nii')
                fn = strcat(fn, '.nii');
            end
            c_ = this.cifti_last;
            c_.cdata = c1_data;
            if contains(fn, '.dtseries')
                c_.diminfo{2} = cifti_diminfo_make_series(sz(2), 0, this.ihcp_.tr, 'SECOND');
            else
                c_.diminfo{2} = cifti_diminfo_make_scalars(1);
            end
            cifti_write(c_, convertStringsToChars(fn));
        end
        function [cdata,cdata1] = write_ciftis(this, cdata, fp, opts)
            arguments
                this mlraut.Cifti
                cdata {mustBeNumericOrLogical}
                fp {mustBeTextScalar}
                opts.do_save_dynamic logical = false
                opts.averaging_method function_handle = @median
            end

            try
                cdata = this.ihcp_.build_final_normalization(cdata);
                if opts.do_save_dynamic
                    this.write_cifti(cdata, fp); % mlraut.HCP
                end
                cdata1 = opts.averaging_method(cdata, 1);
                fp1 = strcat(fp, '_avgt');
                this.write_cifti(cdata1, fp1); % mlraut.HCP
            catch ME
                handwarning(ME)
            end
        end
        function ic = write_nii(this, img, fp)
            arguments
                this mlraut.Cifti
                img {mustBeNumericOrLogical}
                fp {mustBeTextScalar}
            end

            try
                ifc = mlfourd.ImagingFormatContext2(this.template_niigz);
                ifc.img = img;
                [pth,fp] = myfileparts(fp);
                if isempty(pth) || "" == pth
                    pth = this.out_dir;
                end
                ifc.filepath = pth;
                ifc.fileprefix = fp;
                tr = this.ihcp_.tr;
                Nt = length(img);
                ifc.json_metadata.timesMid = 0:tr:tr*(Nt - 1);
                ifc.save();

                ic = mlfourd.ImagingContext2(ifc);
            catch ME
                handwarning(ME)
            end
        end

        function this = Cifti(ihcp)
            arguments
                ihcp mlraut.HCP
            end

            this.ihcp_ = ihcp;
        end
    end

    methods (Static)
    end

    %% PROTECTED

    properties (Access = protected)
        ihcp_
    end

    %% HIDDEN

    methods (Hidden)
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
