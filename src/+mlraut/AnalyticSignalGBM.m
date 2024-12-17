classdef AnalyticSignalGBM < handle & mlraut.AnalyticSignalHCP
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

        %% set config_hemispheres := "lesionR" to accumulate collection of right-sided tumors.
        %  set config_hemispheres := "lesionR_CE" to use contrast-enhanced
        %  set config_hemispheres := "lesionR_WT" to use whole tumor
        %  set config_hemispheres := "nolesion"
        %  set config_hemispheres := "alllesion"
        %  set config_hemispheres := ""

        config_hemispheres = "lesionR"
    end

    methods
        function this = AnalyticSignalGBM(varargin)
            this = this@mlraut.AnalyticSignalHCP(varargin{:});

            this.current_subject = this.subjects{1};
            %this.tasks_ = {'ses-1_task-rest_run-01_desc-preproc', 'ses-1_task-rest_run-02_desc-preproc'};
            %this.max_frames = 158;
        end

        function build_conc(this)

            %% for AnalyticSignalGBM, concat BOLD runs to have 320 frames available

            ciftis = mglob(fullfile(this.cohort_data.mninonlinear_dir, "Results", "*", "ses-*.dtseries.nii"));
            niftis = mglob(fullfile(this.cohort_data.mninonlinear_dir, "Results", "*", "ses-*.nii.gz"));
            niftis = niftis(~endsWith(niftis, "_avgt.nii.gz"));

            task_dir = this.cohort_data.task_dir;
            ensuredir(task_dir);

            fn = char(fullfile(task_dir, this.current_task + "_Atlas_s0.dtseries.nii"));
            if ~isfile(fn)
                cdata = [];
                len = 0;
                for c = ciftis
                    c1 = cifti_read(c);
                    cdata = [cdata, c1.cdata]; %#ok<AGROW>
                    len = len + c1.diminfo{2}.length;
                end
                c1.cdata = cdata;
                c1.diminfo{2}.length = len;
                cifti_write(c1, fn);
            end

            fn = char(fullfile(task_dir, this.current_task + ".nii.gz"));
            if ~isfile(fn)
                img = [];
                for n = niftis
                    n1 = mlfourd.ImagingFormatContext2(n);
                    if isempty(img)
                        img = [img, n1.img]; %#ok<AGROW>
                    else
                        img = cat(4, img, n1.img);
                    end
                end
                n1.img = img;
                n1.fqfilename = fn;
                save(n1);

                n2 = mlfourd.ImagingContext2(n1);
                n2 = n2.timeAveraged();
                save(n2);
                n2.fileprefix = "SBRef_dc";
                save(n2);
            end

            %% flip maps on T1w for use with sets of central & right-sided GBMs

            mnl_dir = this.cohort_data.mninonlinear_dir;
            ce = mlfourd.ImagingContext2(fullfile(mnl_dir, "CE_on_T1w.nii.gz"));
            ce = flip(ce, 1);
            ce.save;
            wt = mlfourd.ImagingContext2(fullfile(mnl_dir, "WT_on_T1w.nii.gz"));
            wt = flip(wt, 1);
            wt.save;
            t1w = mlfourd.ImagingContext2(fullfile(mnl_dir, "T1w.nii.gz"));
            t1w = flip(t1w, 1);
            t1w.save;
        end

        function build_angles_gt_0(this)

            %% shift phases to start at zero for use with "Videen-style" color spaces in wb_view

            wb_dir = this.out_dir;
            toglob = fullfile(wb_dir, "angle*_avgt.dscalar.nii");
            mg = mglob(toglob);
            for c = mg
                c = char(c);
                c1 = cifti_read(c);
                min_ = min(c1.cdata, [], "all");
                c1.cdata = c1.cdata - min_;
                cifti_write(c1, strrep(c, '_avgt.dscalar.nii', '_shifted_avgt.dscalar.nii'));
            end
        end

        function this = call(this, opts)
            %% CALL all subjects

            arguments
                this mlraut.AnalyticSignalHCP
                opts.do_qc logical = false
            end

            % exclude subjects
            %this.subjects = this.subjects(~contains(this.subjects, '_7T'));
            %this.subjects = this.subjects(~contains(this.subjects, 'sub-'));

            out_dir_ = this.out_dir;
            for s = 1:this.num_sub
                try
                    this.current_subject = this.subjects{s};
                    if ~contains(out_dir_, this.current_subject)
                        proposed_dir = fullfile(out_dir_, this.current_subject);
                        this.out_dir = proposed_dir;
                        ensuredir(proposed_dir);
                    end
                    this.build_conc();  % new for AnalyticSignalGBM
                    this.call_subject(s);
                catch ME
                    handexcept(ME)
                end
            end
        end
        
        function j = jsonread(this)
            j = jsonread(this.cohort_data_.json_fqfn);
        end

        function [bold,isleft] = task_dtseries(this, sub, task)
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
                    bold = this.task_dtseries_simple(sub, task);
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
            bold = this.task_dtseries_simple(sub, task); % bilateral
        end

        function bold = task_dtseries_1hemi(this, sub, task, hemi)
            %% duplicates selected hemi to contralateral hemi
            %  Args:
            %      subj (text)
            %      task (text)
            %      hemi char % in {'l', 'L', 'r', 'R', ''}; selecting 'L' copies 'L' ordinates to 'R' ordinates;
            %                                               selecting '' returns this.task_dtseries_simple().
            %  Returns:
            %      BOLD (numeric):  time x grayordinate

            arguments
                this mlraut.AnalyticSignalGBM
                sub {mustBeTextScalar} = this.current_subject
                task {mustBeTextScalar} = this.current_task
                hemi char = ''
            end

            bold = this.task_dtseries_simple(sub, task);
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

            bold = this.task_dtseries_simple(sub, task);
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
            bold = this.task_dtseries_simple(sub, task); % right & bilateral
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

        function mat = task_dtseries_simple(this, varargin)
            mat = task_dtseries_simple@mlraut.AnalyticSignal(this, varargin{:});
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
                        rout_dir=out_dir, tasks=tasks);
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
