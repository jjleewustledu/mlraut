classdef HCP < handle & mlsystem.IHandle
    %% Supports Ryan Raut's use of the Human Connectome Project.  See also:
    %  https://www.science.org/doi/10.1126/sciadv.abf2709
    %  physio_phase_mapping.m
    %  
    %  Created 29-Nov-2022 00:11:57 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 9.13.0.2105380 (R2022b) Update 2 for MACI64.  Copyright 2022 John J. Lee.
    

    properties        
        max_frames = Inf  % max(num_frames) to enforce, used by omit_late_frames()
    end

    properties (Dependent)
        current_subject  % defers to subjects{1} as needed
        current_task  % defers to tasks{1} as needed
        subjects
        tasks

        bold_data
        cohort_data
        cifti_last  % configures cifti historically
        Fs  % BOLD sampling rate (Hz)
        num_frames
        num_frames_ori  % set by BOLDData.task_dtseries()
        num_frames_to_trim  % used by HCP.task_dtseries, HCP.trim_frames, AnalyticSignal.physio_*(); Ryan used 4
        num_nodes  % set by BOLDData
        out_dir
        root_dir  % HCP data directory
        task_dir  % e.g., subject/MNINonlinear/Results/rfMRI_REST1_RL
        task_dtseries_fqfn
        task_niigz_fqfn
        task_signal_reference_fqfn
        t1w_fqfn
        tr  % sampling interval (s), 0.72 for HCP Y.A., 0.8 for HCP Aging, 2.71 for RT GBM
        waves_dir
        wmparc_fqfn
        workbench_dir
    end

    methods %% GET, SET
        function g = get.bold_data(this)
            g = this.bold_data_;
        end
        function g = get.cohort_data(this)
            g = this.cohort_data_;
        end
        function g = get.current_subject(this)
            if ~isempty(this.current_subject_)
                g = this.current_subject_;
                return
            end
            if ~isempty(this.subjects)
                g = this.subjects{1};
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
                return
            end
            if ~isempty(this.tasks)
                g = this.tasks{1};
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
            this.tasks_ = cellfun(@basename, glob(fullfile(tasks_dir, '*')), UniformOutput=false);
            %this.tasks_ = this.tasks_(~contains(this.tasks_, '7T'));
            this.tasks_ = this.tasks_(~contains(this.tasks_, 'tfMRI'));
            g = this.tasks_;
        end

        function g = get.cifti_last(this)
            if ~isempty(this.cifti_last_)
                g = this.cifti_last_;
                return
            end

            this.cifti_last_ = cifti_read(this.task_dtseries_fqfn);
            g = this.cifti_last_;
        end
        function     set.cifti_last(this, s)
            assert(isstruct(s));
            this.cifti_last_ = s;
        end
        function g = get.Fs(this)
            g = 1/this.tr;
        end
        function g = get.num_frames(this)
            trimmed = this.num_frames_ori - 2*this.num_frames_to_trim;
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
        function b = omit_late_frames(this, b)
            %% Keep frames 1:this.max_frames, following use of trim_frames() to remove this.num_frames_to_trim
            %  from start and end of frames, for purposes of omitting brain/cognitive responses to start and conclusion 
            %  of the scanning session.

            bound = min(this.max_frames, size(b, 1));
            b = b(1:bound, :);
        end
        function mat = task_dtseries(this, sub, task)
            %  Args:
            %      this mlraut.HCP
            %      sub {mustBeTextScalar} = this.current_subject
            %      task {mustBeTextScalar} = this.current_task
            %  Returns:
            %      mat (numeric):  time x grayordinate from BOLDData

            arguments
                this mlraut.HCP
                sub {mustBeTextScalar} = this.current_subject
                task {mustBeTextScalar} = this.current_task
            end
            this.current_subject = sub;
            this.current_task = task;

            mat = this.bold_data_.task_dtseries();
        end
        function ic = task_niigz(this)
            ic = this.bold_data_.task_niigz();
        end
        function ic = task_signal_reference(this)
            ic = this.bold_data_.task_signal_reference();
        end

        function this = HCP(opts)
            %%
            %  Args:
            %      opts.max_frames double = Inf
            %      opts.subjects cell {mustBeText} = {}
            %      opts.tasks cell {mustBeText} = {}

            arguments
                opts.max_frames double = Inf
                opts.subjects cell {mustBeText} = {}
                opts.tasks cell {mustBeText} = {}
            end
            this.max_frames = opts.max_frames;
            this.subjects_ = opts.subjects;
            this.tasks_ = opts.tasks;
            this.bold_data_ = mlraut.BOLDData(this);
            this.cohort_data_ = mlraut.CohortData.create(this);

            %% DEPRECATED
            % if isempty(this.subjects_)
            %     this.subjects_ = cellfun(@(x) basename(x), ...
            %         glob(fullfile(this.root_dir, '*'))', UniformOutput=false);
            % end
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
        cifti_last_
        cohort_data_
        subjects_
        tasks_
    end

    methods (Access = protected)
        function that = copyElement(this)
            that = copyElement@matlab.mixin.Copyable(this);
            if ~isempty(this.bold_data_)
                that.bold_data_ = copy(this.bold_data_); end
            if ~isempty(this.cohort_data_)
                that.cohort_data_ = copy(this.cohort_data_); end
        end
    end

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
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
