classdef Cifti < handle & mlsystem.IHandle
    %% line1
    %  line2
    %  
    %  Created 10-Feb-2024 13:47:24 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    

    properties (Dependent)
        is_7T
        out_dir
        template_cifti  % struct for CIFTI
    end

    methods %% GET
        function g = get.is_7T(this)
            g = contains(this.ihcp_.task, '7T');
        end

        function g = get.out_dir(this)
            g = this.ihcp_.out_dir;
        end

        function g = get.template_cifti(this)
            if ~isempty(this.template_cifti_) && isstruct(this.template_cifti_)
                g = this.template_cifti_;
                return
            end

            g = cifti_read(this.ihcp_.task_ref_dscalar_fqfn);
            this.template_cifti_ = g;
        end

        function     set.template_cifti(this, s)
            assert(isstruct(s))
            this.template_cifti_ = s;
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
        
        function fqfn = average_times(~, fqfn0)
            cii = cifti_read(fqfn0);
            cii.cdata = mean(cii.cdata, 2);  % cifti ~ Ngo x Nt
            cii.diminfo{2} = cifti_diminfo_make_scalars(1);
            fqfn = strrep(fqfn0, ".dscalar.nii", "_avgt.dscalar.nii");
            cifti_write(cii, convertStringsToChars(fqfn));
        end
        
        function cii = write_cifti(this, c1_data, fn)
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

            fn = strrep(fn, 'sub-sub', 'sub');  % clean-up legacy nomenclature
            fn = strrep(fn, 'ses-ses', 'ses');

            cii = this.template_cifti;
            cii.cdata = c1_data;
            if contains(fn, '.dtseries')
                cii.diminfo{2} = cifti_diminfo_make_series(sz(2), 0, this.ihcp_.tr, 'SECOND');
            else
                cii.diminfo{2} = cifti_diminfo_make_scalars(1);
            end
            cifti_write(cii, convertStringsToChars(fn));
        end 
        
        function [cdata,cdata1] = write_ciftis(this, cdata, fp, opts)
            arguments
                this mlraut.Cifti
                cdata {mustBeNumericOrLogical}
                fp {mustBeTextScalar}
                opts.partitions logical = []
                opts.do_save_dynamic logical = false
                opts.averaging_method function_handle = @mean
            end

            try
                if opts.do_save_dynamic
                    this.write_cifti(cdata, fp); % mlraut.HCP
                end

                if ~isempty(opts.partitions)
                    cdata1 = cdata(opts.partitions, :);  % select times of interest
                    cdata1 = opts.averaging_method(cdata1, 1);  % average over t
                    fp1 = strcat(fp, '_avgt');
                    this.write_cifti(cdata1, fp1); 
                else
                    cdata1 = opts.averaging_method(cdata, 1);  % average over t
                    fp1 = strcat(fp, '_avgt');   
                    this.write_cifti(cdata1, fp1);  
                end
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
        template_cifti_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
