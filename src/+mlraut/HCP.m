classdef HCP < handle
    %% Supports Ryan Raut's use of the Human Connectome Project.  See also:
    %  https://www.science.org/doi/10.1126/sciadv.abf2709
    %  physio_phase_mapping.m
    %  
    %  Created 29-Nov-2022 00:11:57 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 9.13.0.2105380 (R2022b) Update 2 for MACI64.  Copyright 2022 John J. Lee.
    
    properties (Abstract)
        HCP_signals
        subjects
        tasks
    end

    properties (Constant)
        RSN_NAMES = ...
            {'visual', 'somatomotor', 'dorsal attention', 'ventral attention', 'limbic', ...
             'frontoparietal', 'default mode', ...
             'task+', 'task-'}
    end

    properties
        num_frames_to_trim = 60; % Ryan used 4
        num_nodes = 91282; % HCP standard 2mm "grayordinates"
        tr = 0.72; % sampling interval (s)
    end

    properties (Dependent)
        cifti_last % configures cifti historically
        Fs % BOLD sampling rate (Hz)
        mask_cbm_HCP
        mask_ctx_HCP
        mask_str_HCP
        mask_thal_HCP
        masks_HCP
        networks_HCP
        num_frames
        num_nets
        out_dir
        root_dir % HCP data directory
        waves_dir
        workbench_dir
    end

    methods

        %% GET

        function g = get.cifti_last(this)
            if ~isempty(this.cifti_last_)
                g = this.cifti_last_;
                return
            end
            this.cifti_last_ = cifti_read(this.task_dtseries(this.subjects{1},this.tasks{1}));
            g = this.cifti_last_;
        end
        function     set.cifti_last(this, s)
            assert(isstruct(s));
            this.cifti_last_ = s;
        end
        function g = get.Fs(this)
            g = 1/this.tr;
        end
        function g = get.mask_cbm_HCP(this)
            if ~isempty(this.mask_cbm_HCP_)
                g = this.mask_cbm_HCP_;
                return
            end
            ld = load(fullfile(this.waves_dir, 'supporting_files', 'mask_cbm_HCP.mat'));
            this.mask_cbm_HCP_ = ld.mask_cbm;
            g = this.mask_cbm_HCP_;
        end
        function g = get.mask_ctx_HCP(this)
            if ~isempty(this.mask_ctx_HCP_)
                g = this.mask_ctx_HCP_;
                return
            end
            ld = load(fullfile(this.waves_dir, 'supporting_files', 'mask_ctx_HCP.mat'));
            this.mask_ctx_HCP_ = ld.mask_ctx;
            g = this.mask_ctx_HCP_;
        end
        function g = get.mask_str_HCP(this)
            if ~isempty(this.mask_str_HCP_)
                g = this.mask_str_HCP_;
                return
            end
            ld = load(fullfile(this.waves_dir, 'supporting_files', 'mask_str_HCP.mat'));
            this.mask_str_HCP_ = ld.mask_str;
            g = this.mask_str_HCP_;
        end
        function g = get.mask_thal_HCP(this)
            if ~isempty(this.mask_thal_HCP_)
                g = this.mask_thal_HCP_;
                return
            end
            ld = load(fullfile(this.waves_dir, 'supporting_files', 'mask_thal_HCP.mat'));
            this.mask_thal_HCP_ = ld.mask_thal;
            g = this.mask_thal_HCP_;
        end        
        function g = get.masks_HCP(this)
            if ~isempty(this.masks_HCP_)
                g = this.masks_HCP_;
                return
            end
            ld = load(fullfile(this.waves_dir, 'supporting_files', 'mask_cbm_HCP.mat'));
            this.masks_HCP_.cbm = ld.mask_cbm;
            ld = load(fullfile(this.waves_dir, 'supporting_files', 'mask_ctx_HCP.mat'));
            this.masks_HCP_.ctx = ld.mask_ctx;
            ld = load(fullfile(this.waves_dir, 'supporting_files', 'mask_str_HCP.mat'));
            this.masks_HCP_.str = ld.mask_str;
            ld = load(fullfile(this.waves_dir, 'supporting_files', 'mask_thal_HCP.mat'));
            this.masks_HCP_.thal = ld.mask_thal;
            g = this.masks_HCP_;
        end
        function g = get.networks_HCP(this)
            if ~isempty(this.networks_HCP_)
                g = this.networks_HCP_;
                return
            end
            ld = load(fullfile(this.waves_dir, 'supporting_files', 'networks_HCP.mat'));
            this.networks_HCP_ = ld.assns2;
            g = this.networks_HCP_;
        end
        function g = get.num_nets(this)
            g = length(this.RSN_NAMES);
        end
        function g = get.num_frames(this)
            g = 1200 - 2*this.num_frames_to_trim - 1;
        end
        function g = get.out_dir(this)
            if isempty(this.out_dir_)
                g = pwd;
                return
            end
            g = this.out_dir_;
        end
        function     set.out_dir(this, s)
            assert(isfolder(s))
            this.out_dir_ = s;
        end
        function g = get.root_dir(~)
            if contains(hostname, 'precuneal')
                g = '/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200';
                assert(isfolder(g));
                return
            end
            if contains(hostname, 'pascal')
                g = '/mnt/CHPC_hcpdb/packages/unzip/HCP_1200';
                assert(isfolder(g));
                return
            end
            if contains(hostname, 'cluster')
                g = '/ceph/hcpdb/packages/unzip/HCP_1200';
                assert(isfolder(g));
                return
            end
            error('mlraut:NotImplementedError', 'HCP.get.subjects');
        end
        function g = get.waves_dir(~)
            g = fullfile(getenv('HOME'), 'MATLAB-Drive', 'arousal-waves-main', '');
        end
        function g = get.workbench_dir(~)
            if contains(hostname, 'precuneal')
                g = '/Applications/workbench/bin_macosx64';
                return
            end
            if strcmp(computer, 'GLNXA64')
                [~,r] = system('cat /etc/os-release | grep PRETTY_NAME');
                if contains(r, 'Rocky') || contains(r, 'CentOS')
                    g = '/data/nil-bluearc/raichle/jjlee/.local/workbench-rh_linux64-v1.5.0/workbench/exe_rh_linux64';
                end
                if contains(r, 'Ubuntu')
                    g = '/data/nil-bluearc/raichle/jjlee/.local/workbench-linux64-v1.5.0/workbench/bin_linux64';
                end
                return
            end
            error('mlraut:NotImplementedError', 'HCP.get.workbench_dir');
        end

        %%

        function fn = aparc_a2009s_label_gii(this, sub, hemi)
            %  e.g.:
            %  cd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200/100307/MNINonLinear/fsaverage_LR32k')
            %  g = gifti('100307.L.aparc.a2009s.32k_fs_LR.label.gii')

            if contains(hemi, 'L', IgnoreCase=true)
                hemi = 'L';
            end
            if contains(hemi, 'R', IgnoreCase=true)
                hemi = 'R';
            end

            pth = fullfile(this.root_dir, sub, 'MNINonLinear', 'fsaverage_LR32k');
            fn = sprintf('%s.%s.aparc.a2009s.32k_fs_LR.label.gii', sub, hemi);
            fn = fullfile(pth, fn);
        end
        function fn = aparc_a2009s_dlabel_nii(this, sub)
            %  e.g.:
            %  cd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200/100307/MNINonLinear/fsaverage_LR32k')
            %  g = gifti('100307.aparc.a2009s.32k_fs_LR.dlabel.gii')

            pth = fullfile(this.root_dir, sub, 'MNINonLinear', 'fsaverage_LR32k');
            fn = sprintf('%s.aparc.a2009s.32k_fs_LR.dlabel.nii', sub);
            fn = fullfile(pth, fn);
        end
        function bold = bold_fs_parcel(this, sub, task, parc)
            bold = this.task_dtseries(sub, task);
            mask = this.mask_fs_parcel(sub, parc);
            bold = mean(bold(:, mask), 2, 'omitnan');
        end
        function dd = data_dir(this, subj, task)
            assert(istext(subj));
            assert(istext(task));
            dd = fullfile(this.root_dir, subj, 'MNINonLinear', 'Results', task);
        end
        function m = mask_fs_parcel(this, sub, parc)
            %  sub (text):  e.g., 100307
            %  parc (char):  e.g., 'L_S_parieto_occipital'

            hemi = parc(1);
            g = gifti(this.aparc_a2009s_label_gii(sub, hemi));
            tf = matches(g.labels.name, parc);
            key = g.labels.key(tf');
            m = false(this.num_nodes, 1);
            switch hemi
                case 'L'                   
                    m(1:32492,1) = g.cdata == key;
                case 'R'                  
                    m(32493:64984) = g.cdata == key;
                otherwise 
                    error('mlraut:ValueError', 'HCP.mask_fs_parcel.hemi->%s', hemi);
            end
        end
        function bold = task_dtseries(this, subj, task)
            %  Args:
            %      subj (text)
            %      task (text)
            %  Returns:
            %      BOLD (numeric):  time x grayordinate

            assert(istext(subj));
            assert(istext(task));
            fqfn = fullfile( ...
                this.data_dir(subj, task), strcat(task, '_Atlas_MSMAll_hp2000_clean.dtseries.nii'));
            if ~isfile(fqfn)
                fqfn = fullfile( ...
                    this.data_dir(subj, task), strcat(task, '_Atlas_MSMAll.dtseries.nii'));
            end

            try
                cifti = cifti_read(fqfn);
                this.cifti_last_ = cifti;
                bold = cifti.cdata';
                    
                nt = this.num_frames_to_trim + 1;
                bold = bold(nt:end-nt,:);
            catch ME
                disp([subj ' ' task ' missing:']);
                disp(ME)
                bold = [];
                return
            end
        end
        function write_cifti(this, c1_data, fn)
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
                c_.diminfo{2} = cifti_diminfo_make_series(sz(2));
            else
                c_.diminfo{2} = cifti_diminfo_make_scalars(sz(2));
            end
            cifti_write(c_, fn);
        end
    end

    %% PRIVATE

    properties (Access = private)
        cifti_last_
        mask_ctx_HCP_
        mask_str_HCP_
        mask_thal_HCP_
        mask_cbm_HCP_
        masks_HCP_
        networks_HCP_
        out_dir_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
