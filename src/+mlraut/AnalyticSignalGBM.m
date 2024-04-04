classdef AnalyticSignalGBM < handle & mlraut.AnalyticSignal
    %% line1
    %  line2
    %  
    %  Created 11-Apr-2023 22:47:06 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 9.14.0.2206163 (R2023a) for MACI64.  Copyright 2023 John J. Lee.
    
    properties (Constant)
        % containing valid CE
        SUBS = [
            %"RT126", ...
            "I3CR0023", "I3CR0111", "I3CR0116", "I3CR0211", "I3CR0311", ...
            "I3CR0413", "I3CR0439", "I3CR0453", "I3CR0464", "I3CR0479", ...
            "I3CR0596", "I3CR0605", ...
            "I3CR0737", "I3CR0867", "I3CR1018", ...
            "I3CR1064", "I3CR1085", "I3CR1088", "I3CR1132", "I3CR1159", ...
            "I3CR1413", "I3CR1488", "I3CR1642", "I3CR1772", ...
            "I3CR1821", "I3CR1824", "I3CR1832", "I3CR0043", ...
            "I3CR1562"]
        % "I3CR0408", "I3CR0570", "I3CR0625", "I3CR0671", "I3CR0786", "I3CR0914", "I3CR1774", 

        ISLEFT = [1 0 0 0 1 ...
                  1 1 0 0 0 ...
                  1 1 ...
                  1 1 0 ...
                  1 0 0 0 1 ...
                  1 0 1 0 ...
                  1 0 0 1 ...
                  0]
    end

    properties
        config_hemispheres = "lesionR"
    end

    methods
        function this = AnalyticSignalGBM(varargin)
            this = this@mlraut.AnalyticSignal(varargin{:});

            this.current_subject = this.subjects{1};
            %this.tasks_ = {'ses-1_task-rest_run-01_desc-preproc', 'ses-1_task-rest_run-02_desc-preproc'};
            this.max_frames = 158;
            this.plot_range = 1:158;
        end
        
        function this = call_subject(this, s)
            arguments
                this mlraut.AnalyticSignal
                s double
            end

            this.malloc();
            for t = 1:this.num_tasks
                try
                    this.current_task = this.tasks{t};

                    % BOLD
                    try
                        [bold_, isleft] = this.task_dtseries();
                        bold_ = this.build_global_signal_regressed(bold_);
                        bold_ = hilbert(bold_);
                    catch ME
                        disp([this.current_subject ' ' this.current_task ' BOLD missing or defective:']);
                        handwarning(ME)
                        continue
                    end

                    % Physio
                    try
                        physio_ = this.task_physio(bold_, flipLR=isleft);
                        assert(~isempty(physio_))
                        physio_ = this.build_global_signal_regressed(physio_);
                        physio_ = hilbert(physio_);
                    catch ME
                        disp([this.current_subject ' ' this.current_task ' physio missing or defective:']);
                        handwarning(ME)
                        continue
                    end

                    % Store BOLD, Physio, and Analytic signals
                    this.bold_signal_ = this.build_band_passed(this.build_centered_and_rescaled(bold_));
                    this.bold_signal_ = this.build_final_normalization(this.bold_signal_);
                    if ~all(physio_ == 0)
                        this.physio_signal_ = this.build_band_passed(this.build_centered_and_rescaled(physio_));
                        this.physio_signal_ = this.build_final_normalization(this.physio_signal_);
                        % <psi_p|BOLD_operator|psi_p> ~ <psi_p|psi_b>, not unitary
                        this.analytic_signal_ = this.build_band_passed( ...
                            this.build_centered_and_rescaled(conj(physio_)) .* ...
                            this.build_centered_and_rescaled(bold_));
                        this.analytic_signal_ = this.build_final_normalization(this.analytic_signal_);
                    else
                        this.physio_signal_ = physio_;
                        this.analytic_signal_ = this.bold_signal_;
                    end
                    
                    % Averages for networks
                    this.average_network_signals(this.analytic_signal_);

                    % Store reduced analytic signal, real(), imag(), abs(), angle()
                    if this.do_save
                        % Store reduced analytic signals for all s, t
                        % grid of data from s, t may be assessed with stats
                        save(this, s, t);
                    end
                    if this.do_save_ciftis
                        this.write_ciftis( ...
                            abs(this.analytic_signal_), ...
                            sprintf('abs_as_sub-%s_ses-%s_%s', this.subjects{s}, this.tasks{t}, this.tags), ...
                            do_save_dynamic=this.do_save_dynamic);
                        this.write_ciftis( ...
                            angle(this.analytic_signal_), ...
                            sprintf('angle_as_sub-%s_ses-%s_%s', this.subjects{s}, this.tasks{t}, this.tags), ...
                            do_save_dynamic=this.do_save_dynamic);

                        % analytic_signal_ - bold_signal_
                        diff_ = this.analytic_signal_ - this.bold_signal_;
                        this.write_ciftis( ...
                            abs(diff_), ...
                            sprintf('abs_diff_sub-%s_ses-%s_%s', this.subjects{s}, this.tasks{t}, this.tags), ...
                            do_save_dynamic=this.do_save_dynamic);
                        this.write_ciftis( ...
                            angle(diff_), ...
                            sprintf('angle_diff_sub-%s_ses-%s_%s', this.subjects{s}, this.tasks{t}, this.tags), ...
                            do_save_dynamic=this.do_save_dynamic);
                    end
                catch ME
                    handwarning(ME)
                end
            end
        end
        function j = jsonread(this)
            j = jsonread(this.cohort_data_.json_fqfn);
        end
        function [mat,isleft] = task_dtseries(this, varargin)
            %  Args:
            %      this mlraut.AnalyticSignal
            %      sub {mustBeTextScalar} = this.current_subject
            %      task {mustBeTextScalar} = this.current_task
            %  Returns:
            %      mat (numeric):  time x grayordinate from BOLDData

            [mat,isleft] = this.task_dtseries_gbm(varargin{:});
            mat = this.trim_frames(mat);
            mat = this.omit_late_frames(mat);
            mat = single(mat);
        end
        function [bold,isleft] = task_dtseries_gbm(this, sub, task)
            %  Args:
            %      subj (text)
            %      task (text)
            %  Returns:
            %      BOLD (numeric):  time x grayordinate
            %      isleft (logical):  flip(imaging_context, 1)

            arguments
                this mlraut.AnalyticSignalGBM
                sub {mustBeTextScalar} = this.current_subject
                task {mustBeTextScalar} = this.current_task
            end
            isleft = false;
            switch this.config_hemispheres
                case "lesionR"
                    [bold,isleft] = this.task_dtseries_lesionR(sub, task);
                case "lesionR-CE"
                    [bold,isleft] = this.task_dtseries_lesionR(sub, task);
                case "lesionR-WT"
                    [bold,isleft] = this.task_dtseries_lesionR(sub, task);
                case "nolesion"
                    bold = this.task_dtseries_nolesion(sub, task);
                case "alllesion"
                    bold = this.task_dtseries_alllesion(sub, task);
                otherwise
                    bold = this.task_dtseries(sub, task);
            end
        end
        function bold = task_dtseries_1hemi(this, sub, task, hemi)
            %% duplicates selected hemi to contralateral hemi
            %  Args:
            %      subj (text)
            %      task (text)
            %      hemi char % in {'l', 'L', 'r', 'R', ''}; selecting 'L' copies 'L' ordinates to 'R' ordinates;
            %                                               selecting '' returns this.task_dtseries().
            %  Returns:
            %      BOLD (numeric):  time x grayordinate

            arguments
                this mlraut.AnalyticSignalGBM
                sub {mustBeTextScalar} = this.current_subject
                task {mustBeTextScalar} = this.current_task
                hemi char = ''
            end

            bold = this.task_dtseries(sub, task);
            assert(size(bold, 2) == this.num_nodes, stackstr())
            l_ordinates = 1:29706;
            r_ordinates = 29707:59412;
            switch upper(hemi)
                case 'L'
                    bold(:,r_ordinates) = bold(:,l_ordinates);
                case 'R'
                    bold(:,l_ordinates) = bold(:,r_ordinates);
                otherwise
            end
        end
        function bold = task_dtseries_flipLR(this, sub, task)
            %  Args:
            %      subj (text)
            %      task (text)
            %  Returns:
            %      BOLD (numeric):  time x grayordinate

            arguments
                this mlraut.AnalyticSignalGBM
                sub {mustBeTextScalar} = this.current_subject
                task {mustBeTextScalar} = this.current_task
            end

            bold = this.task_dtseries(sub, task);
            assert(size(bold, 2) == this.num_nodes, stackstr())
            l_ordinates = 1:29706;
            r_ordinates = 29707:59412;
            buff = bold(:, l_ordinates);
            bold(:, l_ordinates) = bold(:, r_ordinates);
            bold(:, r_ordinates) = buff;
        end
        function [bold,isleft] = task_dtseries_lesionR(this, sub, task)
            isleft = contains(this.json.location, "left");
            if isleft
                bold = this.task_dtseries_flipLR(sub, task);
                return
            end
            bold = this.task_dtseries(sub, task); % right & bilateral
        end
        function bold = task_dtseries_nolesion(this, sub, task)
            isleft = contains(this.json.location, "left");
            isright = contains(this.json.location, "right");
            isbl = contains(this.json.location, "b/l");
            if isleft
                bold = this.task_dtseries_1hemi(sub, task, 'R');
                return
            end
            if isright
                bold = this.task_dtseries_1hemi(sub, task, 'L');
                return
            end
            if isbl
                ME = MException("mlraut:RunTimeException", "%s: is bilateral", stackstr());
                throw(ME);
            end
        end
        function bold = task_dtseries_alllesion(this, sub, task)
            isleft = contains(this.json.location, "left");
            isright = contains(this.json.location, "right");
            if isleft
                bold = this.task_dtseries_1hemi(sub, task, 'L');
                return
            end
            if isright
                bold = this.task_dtseries_1hemi(sub, task, 'R');
                return
            end
            bold = this.task_dtseries(sub, task); % bilateral
        end
        function ic = task_signal_mask(this)
            ic = this.task_signal_reference();
            ic = ic.blurred(7).thresh(100).binarized();
        end
    end

    methods (Static)
        function flirt_tumor_segs()
            SUBS = mlraut.AnalyticSignalGBM.SUBS;
            SEGS = {'WT', 'CE', 'TC'};
            srcdir = "/home/usr/jjlee/Tmp";
            workdir = "/home/usr/jjlee/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify"; % "sub-I3CR0023/MNINonLinear"
            pwd0 = pushd(workdir);
            %parfor (idx = 1:length(SUBS), 8)
            for idx = 1:length(SUBS)
                try                    
                    src = fullfile(srcdir, SUBS{idx}, "atlas");
                    
                    % work with Kiyun's GBM imaging in Tmp
                    pwd1 = pushd(src);
    
                    % flirt
                    % /usr/local/fsl/bin/flirt -in /Volumes/PrecunealSSD2/AnalyticSignalGBM/analytic_signal/sourcedata/I3CR0023/atlas/T1w.nii.gz -ref /Volumes/PrecunealSSD2/AnalyticSignalGBM/analytic_signal/sourcedata/I3CR0023/atlas/I3CR0023_orient-rpi_mpr1.nii.gz -out /Volumes/PrecunealSSD2/AnalyticSignalGBM/analytic_signal/sourcedata/I3CR0023/atlas/T1w_on_mpr.nii.gz -omat /Volumes/PrecunealSSD2/AnalyticSignalGBM/analytic_signal/sourcedata/I3CR0023/atlas/T1w_on_mpr.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear

                    cmd_flirt = "/usr/local/fsl/bin/flirt";
                    in = fullfile(workdir, "sub-"+SUBS{idx}, "MNINonLinear", "T1w.nii.gz");
                    gref = glob('*_orient-rpi_mpr1.nii.gz');
                    gref = gref(contains(gref, SUBS{idx}));
                    ref = gref{1};
                    out = "T1w_on_mpr.nii.gz";
                    omat = "T1w_on_mpr.mat";
                    opts = "-bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear";
                    mysystem(sprintf("%s -in %s -ref %s -out %s -omat %s %s", cmd_flirt, in, ref, out, omat, opts));
        
                    % invert xfm
                    % /usr/local/fsl/bin/convert_xfm -omat /Volumes/PrecunealSSD2/AnalyticSignalGBM/analytic_signal/sourcedata/I3CR0023/atlas/mpr_on_T1w.mat -inverse /Volumes/PrecunealSSD2/AnalyticSignalGBM/analytic_signal/sourcedata/I3CR0023/atlas/T1w_on_mpr.mat
                    cmd_convert = "/usr/local/fsl/bin/convert_xfm";
                    omat_inv = "mpr_on_T1w.mat";
                    mysystem(sprintf("%s -omat %s -inverse %s", cmd_convert, omat_inv, omat));
        
                    % applyxfm
                    % /usr/local/fsl/bin/flirt -in /Volumes/PrecunealSSD2/AnalyticSignalGBM/analytic_signal/sourcedata/I3CR0023/atlas/seg_WT_on_mpr1.nii.gz -applyxfm -init /Volumes/PrecunealSSD2/AnalyticSignalGBM/analytic_signal/sourcedata/I3CR0023/atlas/mpr_on_T1w.mat -out /Volumes/PrecunealSSD2/AnalyticSignalGBM/analytic_signal/sourcedata/I3CR0023/atlas/WT_on_T1w.nii.gz -paddingsize 0.0 -interp trilinear -ref /Volumes/PrecunealSSD2/AnalyticSignalGBM/analytic_signal/sourcedata/I3CR0023/atlas/T1w.nii.gz
                    for sidx = 1:length(SEGS)
                        try
                            gin_seg = glob(convertStringsToChars(sprintf("*seg*%s_on_orient-rpi_mpr1.nii.gz", SEGS{sidx})));
                            if isempty(gin_seg)
                                gin_seg = glob(convertStringsToChars(sprintf("*%s_on_orient-rpi_mpr1.nii.gz", SEGS{sidx})));                           ; 
                            end
                            if isempty(gin_seg); continue; end
                            in_seg = gin_seg{1};
                            out_on_T1w = SEGS{sidx} + "_on_T1w.nii.gz";
                            opts1 = "-paddingsize 0.0 -interp nearestneighbour";
                            mysystem(sprintf("%s -in %s -applyxfm -init %s -out %s %s -ref %s", cmd_flirt, in_seg, omat_inv, out_on_T1w, opts1, in));

                            % erode segmentations
                            ic = mlfourd.ImagingContext2(out_on_T1w);
                            ic = ic.thresh(1-eps('single'));
                            ic = ic.binarized();
                            ic.filename = out_on_T1w;
                            ic.filepath = myfileparts(in);
                            ic.save();
                        catch ME
                            handwarning(ME)
                        end
                    end

                    popd(pwd1);
                catch ME
                    handwarning(ME)
                end
            end            
            popd(pwd0);
        end
        function prepare_tumor_segs()
            SUBS = mlraut.AnalyticSignalGBM.SUBS;
            pwd0 = pushd(fullfile("/data/nil-bluearc/shimony/jjlee/Kiyun/rsFC_PreProc"));
            %parfor (idx = 1:length(SUBS), 16)
            for idx = 1:length(SUBS)
                try
                    g = glob(convertStringsToChars(SUBS(idx)+"*"));
                    assert(~isempty(g));
                    ori = g{1};
                    targ = fullfile("/home/usr/jjlee/Tmp", SUBS(idx), "atlas");
                    ensuredir(targ);
                    copyfile(fullfile(ori, "atlas", "*"), targ);
    
                    % arrange folder atlas
                    pwd1 = pushd(targ);
                    mysystem("gzip *.nii");
                    ensuredir("4dfp");
                    movefile("*_t4", "4dfp");
                    movefile("*.4dfp.*", "4dfp");

                    % movefile will often fail; need protection within try/catch
                    %ensuredir("log");
                    %movefile("*.log", "log");
                    %ensuredir("lst");
                    %movefile("*.lst", "lst");
                    %ensuredir("pve");
                    %movefile("*_pve*", "pve");
                    %ensuredir("syn");
                    %movefile("*Warp*", "syn");
                    %movefile("*0GenericAffine*", "syn");
    
                    % .4dfp.* -> .nii.gz
                    pwd2 = pushd("4dfp");
                    g = asrow(glob("*.4dfp.hdr"));
                    for gidx = 1:length(g)
                        fp = mybasename(g{gidx});
                        try
                            mysystem(sprintf("nifti_4dfp -n %s %s", fp, fp));
                            mysystem(sprintf("gzip %s.nii", fp));
                        catch ME
                            handwarning(ME)
                        end
                    end
                    mysystem(sprintf("mv -f *.nii.gz %s", targ))
                    popd(pwd2);
    
                    % AFNI 3dresample *_mpr1.nii.gz
                    try
                        g = asrow(glob("*_mpr1.nii.gz"));
                        for gidx = 1:length(g)
                            mlpipeline.Bids.afni_3dresample(g{gidx});
                        end
                    catch ME
                        handwarning(ME)
                    end
                    try
                        g = asrow(glob("*WT_on_mpr1.nii.gz"));
                        for gidx = 1:length(g)
                            mlpipeline.Bids.afni_3dresample(g{gidx});
                        end
                    catch ME
                        handwarning(ME)
                    end
                    try
                        g = asrow(glob("*CE_on_mpr1.nii.gz"));
                        for gidx = 1:length(g)
                            mlpipeline.Bids.afni_3dresample(g{gidx});
                        end
                    catch ME
                        handwarning(ME)
                    end
                    try
                        g = asrow(glob("*TC_on_mpr1.nii.gz"));
                        for gidx = 1:length(g)
                            mlpipeline.Bids.afni_3dresample(g{gidx});
                        end
                    catch ME
                        handwarning(ME)
                    end
    
                    popd(pwd1);
                catch ME
                    handwarning(ME)
                end
            end
            popd(pwd0);
        end
        function serialcall(opts)
            arguments
                opts.config_hemispheres {mustBeTextScalar} = "" % "lesionR" "nolesion" "alllesion" ""
            end

            root_dir = '/home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify';
            out_dir = '/home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalGBM/analytic_signal/matlabout';
            ensuredir(out_dir);
            tasks = {'ses-1_task-rest_run-01_desc-preproc', 'ses-1_task-rest_run-02_desc-preproc'};

            g = glob(fullfile(root_dir, 'sub-*'));
            %g = flip(g);
            g = cellfun(@(x) basename(x), g, UniformOutput=false);
            %g = g(1:end);
            leng = length(g);
            for idxg = 1:leng
                try
                    this = mlraut.AnalyticSignalGBM(subjects=g(idxg), ...
                        root_dir=root_dir, out_dir=out_dir, tasks=tasks);
                    this.config_hemispheres = opts.config_hemispheres; 
                    call(this);
                catch ME
                    handwarning(ME)
                end
            end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
