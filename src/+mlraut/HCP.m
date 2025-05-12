classdef HCP < handle & mlsystem.IHandle
    %% Supports Ryan Raut's use of the Human Connectome Project.  See also:
    %  https://www.science.org/doi/10.1126/sciadv.abf2709
    %  physio_phase_mapping.m
    %  
    %  Created 29-Nov-2022 00:11:57 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 9.13.0.2105380 (R2022b) Update 2 for MACI64.  Copyright 2022 John J. Lee.
    

    properties 
        do_7T
        do_resting
        do_task       
        max_frames  % max(num_frames) to impose
    end

    properties (Dependent)
        current_subject  % defers to subjects{1} as needed
        current_task  % defers to tasks{1} as needed
        subjects
        tasks

        bold_data
        cohort_data
        twistors

        extended_task_dir  % supports HCPAging/rfMRIExtended/fmriresults01/HCA*
        Fs  % BOLD sampling rate (Hz)
        num_frames
        num_frames_ori  % set by BOLDData.task_dtseries()
        num_frames_to_trim  % set by CohortData.num_frames_to_trim
        num_nodes  % set by BOLDData
        out_dir
        root_dir  % HCP data directory
        stats_fqfn
        task_dir  % e.g., subject/MNINonlinear/Results/rfMRI_REST1_RL
        task_dtseries_fqfn
        task_niigz_fqfn
        task_signal_reference_fqfn
        template_cifti
        template_niigz
        thickness_dscalar_fqfn
        t1w_fqfn
        tr  % sampling interval (s), 0.72 for HCP Y.A., 0.8 for HCP Aging, 2.71 for RT GBM
        waves_dir
        wmparc_fqfn
        workbench_dir
    end

    methods %% GET, SET 
        function g = get.current_subject(this)
            if ~isempty(this.current_subject_)
                g = this.current_subject_;
                g = strrep(g, filesep, '');  % kludge for legacy .mat
                return
            end
            if ~isempty(this.subjects)
                g = this.subjects{1};
                g = strrep(g, filesep, '');  % kludge for legacy .mat
                return
            end
            g = [];
        end
        function     set.current_subject(this, s)
            arguments
                this mlraut.HCP
                s {mustBeTextScalar}
            end
            this.current_subject_ = s;
        end
        function g = get.current_task(this)
            if ~isempty(this.current_task_)
                g = this.current_task_;
                g = strrep(g, filesep, '');  % kludge for legacy .mat
                return
            end
            if ~isempty(this.tasks)
                g = this.tasks{1};
                g = strrep(g, filesep, '');  % kludge for legacy .mat
                return
            end
            g = [];
        end
        function     set.current_task(this, s)
            arguments
                this mlraut.HCP
                s {mustBeTextScalar}
            end
            this.current_task_ = s;
        end
        function g = get.subjects(this)
            g = this.subjects_;
        end
        function     set.subjects(this, s)
            arguments
                this mlraut.HCP
                s {mustBeText}
            end
            this.subjects_ = s;
        end 
        function g = get.tasks(this)
            if ~isempty(this.tasks_)
                g = this.tasks_;
                return
            end

            tasks_dir = fullfile(this.root_dir, this.current_subject, 'MNINonLinear', 'Results');
            if ~isfolder(tasks_dir)
                g = this.tasks_;
                return
            end
            this.tasks_ = {};
            if this.do_resting
                this.tasks_ = [ ...
                    this.tasks_; ...
                    cellfun(@basename, glob(fullfile(tasks_dir, 'rfMRI*')), UniformOutput=false)];
                if isempty(this.tasks_)
                    this.tasks_ = [ ...
                        this.tasks_; ...
                        cellfun(@basename, glob(fullfile(tasks_dir, '*rest*')), UniformOutput=false)];
                end
            end
            if this.do_task
                this.tasks_ = [ ...
                    this.tasks_; ...
                    cellfun(@basename, glob(fullfile(tasks_dir, 'tfMRI*')), UniformOutput=false)];
            end
            if this.do_7T
                this.tasks_ = this.tasks_(contains(this.tasks_, '7T'));
            else
                this.tasks_ = this.tasks_(~contains(this.tasks_, '7T'));
            end
            g = this.tasks_;
        end
        function     set.tasks(this, s)
            assert(istext(s))
            this.tasks_ = s;
        end

        function g = get.bold_data(this)
            g = this.bold_data_;
        end
        function g = get.cohort_data(this)
            g = this.cohort_data_;
        end
        function g = get.twistors(this)
            g = this.twistors_;
        end

        function g = get.extended_task_dir(this)
            g = this.cohort_data_.extended_task_dir;
        end
        function g = get.Fs(this)
            g = 1/this.tr;
        end
        function g = get.num_frames(this)
            trimmed = this.num_frames_ori - this.num_frames_to_trim;
            g = min(trimmed, this.max_frames);
        end
        function g = get.num_frames_ori(this)
            g = this.bold_data_.num_frames_ori;
        end
        function g = get.num_frames_to_trim(this)
            g = this.cohort_data_.num_frames_to_trim;
        end
        function g = get.num_nodes(this)
            g = this.bold_data_.num_nodes;
        end
        function g = get.out_dir(this)
            g = this.cohort_data_.out_dir;
        end
        function     set.out_dir(this, s)
            this.cohort_data_.out_dir = s;
        end
        function g = get.root_dir(this)
            g = this.cohort_data_.root_dir;
        end
        function g = get.stats_fqfn(this)
            g = this.cohort_data_.stats_fqfn;
        end
        function g = get.task_dir(this)
            g = this.cohort_data_.task_dir;
        end
        function g = get.task_dtseries_fqfn(this)
            g = this.cohort_data_.task_dtseries_fqfn;
        end
        function g = get.task_niigz_fqfn(this)
            g = this.cohort_data_.task_niigz_fqfn;
        end
        function g = get.task_signal_reference_fqfn(this)
            g = this.cohort_data_.task_signal_reference_fqfn;
        end
        function g = get.template_cifti(this)
            if ~isempty(this.template_cifti_) && isstruct(this.template_cifti_)
                g = this.template_cifti_;
                return
            end
            g = cifti_read(this.thickness_dscalar_fqfn);
            this.template_cifti_ = g;
        end
        function     set.template_cifti(this, s)
            assert(isstruct(s))
            this.template_cifti_ = s;
        end
        function g = get.thickness_dscalar_fqfn(this)
            g = this.cohort_data_.thickness_dscalar_fqfn;
        end
        function g = get.template_niigz(this)
            if ~isempty(this.template_niigz_) && isa(this.template_niigz_, 'mlfourd.ImagingContext2')
                g = this.template_niigz_;
                return
            end
            g = mlfourd.ImagingContext2(this.wmparc_fqfn);
            this.template_niigz_ = g;
        end
        function     set.template_niigz(this, s)
            s = mlfourd.ImagingContext2(s);
            this.template_niigz_ = s;
        end
        function g = get.t1w_fqfn(this)
            g = this.cohort_data_.t1w_fqfn;
        end
        function g = get.tr(this)
            g = this.cohort_data_.tr;
        end
        function g = get.waves_dir(~)
            g = fullfile(getenv('HOME'), 'MATLAB-Drive', 'arousal-waves-main', '');
        end
        function g = get.wmparc_fqfn(this)
            g = this.cohort_data_.wmparc_fqfn;
        end
        function g = get.workbench_dir(~)
            if contains(computer, 'MAC')
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
    end

    methods
        function this = malloc(this)
            %% reset for new tasks or new subjects

            this.bold_data_ = mlraut.BOLDData(this);
            this.cohort_data_ = mlraut.CohortData.create(this);
            this.twistors_ = mlraut.Twistors(this);

            this.task_dtseries_ = [];
            this.task_niigz_ = [];
            this.task_signal_mask_ = [];
            this.task_signal_reference_ = [];
            this.template_cifti_ = [];
            this.template_niigz_ = [];
        end

        function mat = task_dtseries(this, varargin)
            %% supports interface

            if ~isempty(this.task_dtseries_)
                mat = this.task_dtseries_;
                return
            end

            mat = this.task_dtseries_simple(varargin{:});
            this.task_dtseries_ = mat;
        end

        function mat = task_dtseries_simple(this, sub, task, opts)
            %  Args:
            %      this mlraut.HCP
            %      sub {mustBeTextScalar} = this.current_subject
            %      task {mustBeTextScalar} = this.current_task
            %      opts.network_type {mustBeText} = ""
            %  Returns:
            %      mat (numeric):  time x grayordinate from BOLDData

            arguments
                this mlraut.HCP
                sub {mustBeTextScalar} = this.current_subject
                task {mustBeTextScalar} = this.current_task
                opts.network_type {mustBeText} = ""
            end
            this.current_subject = sub;
            this.current_task = task;

            mat = this.bold_data_.task_dtseries();
            mat = this.trim_frames(mat);

            if ~isemptytext(opts.network_type)
                mat = this.average_network_signal(mat, network_type=opts.network_type);
            end
        end

        function ic = task_niigz(this)
            if ~isempty(this.task_niigz_)
                ic = this.task_niigz_;
                return
            end

            ic = this.bold_data_.task_niigz();
            ic = this.trim_frames(ic);
            ic.ensureSingle();
            this.task_niigz_ = ic;
        end

        function ic = task_signal_mask(this)
            if ~isempty(this.task_signal_mask_)
                ic = this.task_signal_mask_;
                return
            end

            ic = this.task_signal_reference();
            ic = ic.binarized();
            
            ic1 = mlfourd.ImagingContext2(this.wmparc_fqfn);
            ic1 = ic1.binarized();
            % ic1 = ic1.blurred(6);
            % ic1 = ic1.thresh(0.1);
            % ic1 = ic1.binarized();
            if ~isempty(getenv("DEBUG"))
                ic1.view_qc(ic);
            end

            % mask should not have greater coverage than binarized wmparc
            if dipsum(ic) > dipsum(ic1)
                ic = ic1;
            end
            ic.ensureSingle();
            this.task_signal_mask_ = ic;
        end

        function ic = task_signal_reference(this)
            if ~isempty(this.task_signal_reference_)
                ic = this.task_signal_reference_;
                return
            end

            ic = this.bold_data_.task_signal_reference();
            ic.ensureSingle();
            this.task_signal_reference_ = ic;
        end

        function tseries = trim_frames(this, tseries)
            idx_start = 1 + this.num_frames_to_trim;
            if isnumeric(tseries)
                tseries = tseries(idx_start:end,:);
                tseries = tseries(1:this.num_frames,:);
                return
            end
            if isa(tseries, "mlfourd.ImagingContext2")
                img = double(tseries);
                img = img(:,:,:,idx_start:end);
                img = img(:,:,:,1:this.num_frames);
                tseries.selectImagingTool(img=img);
                j = tseries.json_metadata;
                j.timesMid = j.timesMid(idx_start:end);
                tseries.addJsonMetadata(j);
                return
            end
            error("mlraut:TypeError", stackstr())
        end

        function this = HCP(opts)
            %%
            %  Args:
            %      opts.max_frames double = Inf
            %      opts.subjects cell {mustBeText} = {}
            %      opts.tasks cell {mustBeText} = {}
            %      opts.network {mustBeText} = ""

            arguments
                opts.max_frames double = Inf
                opts.subjects {mustBeText} = {}
                opts.tasks {mustBeText} = {}
            end
            opts.subjects = strrep(opts.subjects, filesep, "");
            opts.subjects = convertStringsToChars(opts.subjects);
            opts.subjects = ensureCell(opts.subjects);
            opts.tasks = strrep(opts.tasks, filesep, "");
            opts.tasks = convertStringsToChars(opts.tasks);
            opts.tasks = ensureCell(opts.tasks);
            
            this.do_7T = false;
            this.do_resting = true;
            this.do_task = false;

            this.max_frames = opts.max_frames;
            this.subjects_ = opts.subjects;
            this.tasks_ = opts.tasks;

            if ischar(this.subjects_)
                this.subjects_ = {this.subjects_};
            end
            if isstring(this.subjects_)
                this.subjects_ = convertStringsToChars(this.subjects_);
            end
        end
    end

    %% PROTECTED

    properties (Access = protected)
        current_subject_  % defers to subjects{1} as needed
        current_task_  % defers to tasks{1} as needed

        bold_data_
        cohort_data_
        subjects_
        tasks_
        task_dtseries_
        task_niigz_  % cached here and in mlraut.BOLDData
        task_signal_mask_
        task_signal_reference_
        template_cifti_
        template_niigz_
        twistors_
    end

    methods (Access = protected)
        function that = copyElement(this)
            that = copyElement@matlab.mixin.Copyable(this);
            if ~isempty(this.bold_data_)
                that.bold_data_ = copy(this.bold_data_); end
            if ~isempty(this.cohort_data_)
                that.cohort_data_ = copy(this.cohort_data_); end
    %% HIDDEN & DEPRECATED

    methods (Hidden)
        function bold = bold_fs_parcel(this, sub, task, parc)
            arguments
                this mlraut.HCP
                sub {mustBeTextScalar} = this.current_subject
                task {mustBeTextScalar} = this.current_task
                parc char = 'L_G_precuneus'
            end
            bold = this.task_dtseries(sub, task);
            mask = this.mask_fs_parcel(sub, parc);
            bold = mean(bold(:, mask), 2, 'omitnan');
        end

        function m = mask_fs_parcel(this, sub, parc)
            % this mlraut.HCP
            % sub {mustBeTextScalar} = this.current_subject
            % parc double = 15 % fourth ventricle
            % hemi {mustBeTextScalar} = 'LR'

            arguments
                this mlraut.HCP
                sub {mustBeTextScalar} = this.current_subject
                parc char = 'L_G_precuneus'
            end

            hemi = parc(1);
            g = gifti(convertStringsToChars(this.aparc_a2009s_label_gii(sub, hemi)));
            tf = matches(g.labels.name, parc);
            key = g.labels.key(tf');
            m = false(this.num_nodes, 1);
            switch hemi
                case 'L'                   
                    m(1:32492) = g.cdata == key;
                case 'R'                  
                    m(32493:64984) = g.cdata == key;
                otherwise 
                    error('mlraut:ValueError', stackstr());
            end
            if ~isempty(this.task_niigz_)
                that.task_niigz_ = copy(this.task_niigz_); end
            if ~isempty(this.task_signal_mask_)
                that.task_signal_mask_ = copy(this.task_signal_mask_); end
            if ~isempty(this.task_signal_reference_)
                that.task_signal_reference_ = copy(this.task_signal_reference_); end
            if ~isempty(this.twistors_)
                that.twistors_ = copy(this.twistors_); end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
