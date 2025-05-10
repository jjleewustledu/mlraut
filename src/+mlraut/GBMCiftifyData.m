classdef GBMCiftifyData < handle & mlraut.CohortData
    %% line1
    %  line2
    %  
    %  Created 08-Feb-2024 18:48:01 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    

    properties        
    end

    properties (Dependent)
        CE_fqfn
        ciftify_subject_fmri_log_fqfn
        datashare_dir
        json_fqfn
        num_frames_to_trim
        out_dir
        rsFC_PreProc_loc
        root_dir
        stats_fqfn
        task_dtseries_fqfn
        task_niigz_fqfn
        task_signal_reference_fqfn   
        t1w_fqfn
        tr
        wmparc_fqfn
        WT_fqfn

        map_rt_i3cr
        table_excluded
        table_gbm  % excludes table_excluded
        table_rt_i3cr
    end

    methods %% GET
        function g = get.CE_fqfn(this)
            g = fullfile(this.mninonlinear_dir, "CE_on_T1w.nii.gz");  % mm voxels
            assert(isfile(g), stackstr())
        end
        function g = get.ciftify_subject_fmri_log_fqfn(this)
            pth = fileparts(this.task_dtseries_fqfn);
            g = fullfile(pth, "ciftify_subject_fmri.log");
        end
        function g = get.datashare_dir(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "AnalyticSignalGBM", "GBM_datashare");
            assert(isfolder(g))
        end
        function g = get.json_fqfn(this)
            g = fullfile(this.root_dir, this.sub, "gbm.json");  % mm voxels
        end
        function g = get.map_rt_i3cr(this)
            if isempty(this.map_rt_i3cr_)
                i3cr = this.table_rt_i3cr.I3CR;
                found = find(contains(i3cr, 'not', IgnoreCase=true));
                for f = asrow(found)
                    i3cr{f} = 'I3CR_unknown';
                end
                this.map_rt_i3cr_ = containers.Map([this.table_rt_i3cr.RT; i3cr], [i3cr; i3cr]);
            end
            g = this.map_rt_i3cr_;
        end
        function g = get.num_frames_to_trim(this)
            g = 0;
            % g = this.json.num_frames_to_trim;
        end
        function g = get.out_dir(this)
            if ~isemptytext(this.out_dir_)
                g = this.out_dir_;
                return
            end

            g = fullfile(getenv("SINGULARITY_HOME"), "AnalyticSignalGBM", "analytic_signal", "matlabout");
            assert(isfolder(g))
            this.out_dir_ = g;
            % /home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalGBM/analytic_signal/matlabout
        end
        function     set.out_dir(this, s)
            assert(istext(s))
            %ensuredir(s);
            this.out_dir_ = s;
        end
        function g = get.rsFC_PreProc_loc(~)
            g = "jjlee@linux1.neuroimage.wustl.edu:" + ...
                fullfile(filesep, "data", "nil-bluearc", "shimony", "bidhan", "rsFC_PreProc");
        end
        function g = get.root_dir(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "AnalyticSignalGBM", "analytic_signal", "dockerout", "ciftify");
            assert(isfolder(g))
            % /home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify
        end
        function g = get.stats_fqfn(this)
            g = this.task_dtseries_fqfn;
        end
        function g = get.table_excluded(this)
            if isempty(this.table_excluded_)
                ld = load(fullfile(this.datashare_dir, "excluded.mat"));
                this.table_excluded_ = ld.excluded;
            end
            g = this.table_excluded_;
        end
        function g = get.table_gbm(this)
            if isempty(this.table_gbm_)
                ld = load(fullfile(this.datashare_dir, "GBMClinicalDatabasesorted.mat"));
                this.table_gbm_ = ld.GBMClinicalDatabasesorted;
            end
            g = this.table_gbm_;
        end
        function g = get.table_rt_i3cr(this)
            if isempty(this.table_rt_i3cr_)
                ld = load(fullfile(this.datashare_dir, "RT_I3CR.mat"));
                this.table_rt_i3cr_ = ld.RT_I3CR;
            end
            g = this.table_rt_i3cr_;
        end
        function g = get.task_dtseries_fqfn(this)
            mg = mglob(fullfile(this.task_dir, sprintf("%s*_Atlas_s0.dtseries.nii", this.task)));
            assert(~isemptytext(mg), stackstr())
            g = mg(1);
        end
        function g = get.task_niigz_fqfn(this)
            mg = mglob(fullfile(this.task_dir, sprintf("%s*.nii.gz", this.task)));
            assert(~isemptytext(mg), stackstr())
            g = mg(1);
        end
        function g = get.task_signal_reference_fqfn(this)
            mg = mglob(fullfile(this.task_dir, sprintf("%s*.nii.gz", this.task)));
            assert(~isemptytext(mg), stackstr())
            g = mg(1);
        end
        function g = get.t1w_fqfn(this)
            g = fullfile(this.mninonlinear_dir, "T1w.nii.gz");  % mm voxels
            assert(isfile(g), stackstr())
        end
        function g = get.tr(this)
            g = this.json.tr;
        end
        function g = get.wmparc_fqfn(this)
            g = fullfile(this.mninonlinear_dir, "wmparc.nii.gz");  % mm voxels
            assert(isfile(g), stackstr())
        end
        function g = get.WT_fqfn(this)
            g = fullfile(this.mninonlinear_dir, "WT_on_T1w.nii.gz");  % mm voxels
            assert(isfile(g), stackstr())
        end
    end

    methods
        function this = GBMCiftifyData(varargin)            
            this = this@mlraut.CohortData(varargin{:});

            error("mlraut:DeprecationError", stackstr());
        end

        function build_gbm_json(this)
            
            %% parse cifify subject fMRI log

            fileID = fopen(this.ciftify_subject_fmri_log_fqfn, 'r');
            if fileID == -1
                error('mlraut:IOError', 'Failed to open the log file.');
            end

            numTRs = [];
            TR = [];
            while ~feof(fileID)
                line = fgetl(fileID);
                position1 = strfind(line, 'Number of TRs: ');
                position2 = strfind(line, 'TR(ms): ');

                if ~isempty(position1)
                    numTRs = str2double(line(position1 + length('Number of TRs: '):end));
                end
                if ~isempty(position2)
                    TR = str2double(line(position2 + length('TR(ms): '):end));
                    break; % Exit the loop once the numbers are found
                end
            end

            fclose(fileID);
            if isempty(numTRs) || isempty(TR)
                error("mlraut:RuntimeError", ...
                    "Number of TRs or TR not found in %s.", this.ciftify_subject_fmri_log_fqfn);
            end

            %% assemble gbm.json

            j.tr = TR;
            j.num_frames_ori = numTRs;
            j.num_frames_to_trim = 0;
            hemi = lower(this.find_table_value("hemi"));
            loc = lower(this.find_table_value("brain_location"));
            j.location = sprintf("%s %s", hemi, loc);
            if ~isfile(this.json_fqfn)
                jsonwrite(j, this.json_fqfn);
            else
                error("mlraut:IOError", ...
                    "%s aborted because %s already exists.", stackstr(), this.json_fqfn);
            end
        end

        function val = find_table_value(this, var_name)            
            if any(contains(this.table_gbm.Properties.VariableNames, var_name, IgnoreCase=true))
                this_i3cr = strrep(this.ihcp_.current_subject, 'sub-', '');
                if contains(this.map_rt_i3cr.keys, this_i3cr)
                    this_i3cr = this.map_rt_i3cr(this_i3cr);  % rt -> i3cr id
                end
                selected = strcmp(this.table_gbm.I3CRID, this_i3cr);
                val = string(this.table_gbm{selected, var_name});
                assert(isscalar(val))
            end
        end

        function g = surf_gii_fqfn(this, hemis)
            arguments
                this mlraut.GBMCiftifyData
                hemis {mustBeTextScalar} = "L"
            end

            if startsWith(hemis, "L", IgnoreCase=true)
                hemis = "L";
            elseif startsWith(hemis, "R", IgnoreCase=true)
                hemis = "R";
            end

            mg = mglob(fullfile(this.mninonlinear_dir, "fsaverage_LR32k", "*."+hemis+".sphere_MSMall.32k_fs_LR.surf.gii"));
            % 995174.L.midthickness_MSMAll.164k_fs_LR.surf.gii
            if isemptytext(mg)
                mg = mglob(fullfile(this.mninonlinear_dir, "fsaverage_LR32k", "*."+hemis+".sphere.32k_fs_LR.surf.gii"));
                % 995174.L.midthickness.164k_fs_LR.surf.gii
            end
            assert(~isemptytext(mg), stackstr())
            g = mg(end);
        end
    end

    methods (Static)
        function tf = ISLEFT()
            ld = load(fullfile(getenv("SINGULARITY_HOME"), ...
                "AnalyticSignalGBM", "GBM_datashare", "ISLEFT.mat"));
            tf = asrow(double(ld.ISLEFT));
        end

        function flirt_tumor_segs()
            %% use after prepare_tumor_segs

            SUBS = mlraut.AnalyticSignalGBM.SUBS;
            SEGS = {'WT', 'CE'};
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
            %% use before flirt_tumor_segs;
            %  build tumor segs in neuroimage machines, ~/Tmp

            SUBS = mlraut.GBMCiftifyData.SUBS;
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
    
                    popd(pwd1);
                catch ME
                    handwarning(ME)
                end
            end
            popd(pwd0);
        end

        function s = SUBS()
            ld = load(fullfile(getenv("SINGULARITY_HOME"), ...
                "AnalyticSignalGBM", "GBM_datashare", "SUBS.mat"));
            s = asrow(string(ld.SUBS));
        end

    end

    %% PROTECTED

    properties (Access = protected)
        out_dir_
        map_rt_i3cr_   % containers.Map supports text -> char, rt -> i3cr, i3cr -> i3cr
        table_excluded_  % chars
        table_gbm_  % strings
        table_rt_i3cr_  % chars
    end

    methods (Access = protected)
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
