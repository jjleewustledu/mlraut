classdef GBMCiftifyData2 < handle & mlraut.CohortData
    %% line1
    %  line2
    %  
    %  Created 08-Feb-2024 18:48:01 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    

    properties
        num_frames_to_trim = 0
    end

    properties (Dependent)
        CE_fqfn
        CE_ic
        ciftify_subject_fmri_log_fqfn
        datashare_dir
        edema_fqfn
        edema_ic
        json_fqfn
        out_dir
        rsFC_PreProc_loc
        root_dir
        stats_fqfn
        task_dtseries_fqfn
        task_niigz_fqfn
        task_ref_niigz_fqfn
        task_ref_dscalar_fqfn
        thickness_dscalar_fqfn
        t1w_fqfn
        tr
        wmparc_fqfn
        WT_fqfn
        WT_ic

        map_rt_i3cr
        table_excluded
        table_gbm  % excludes table_excluded
        table_rt_i3cr
    end

    methods %% GET
        function g = get.CE_fqfn(this)
            g = fullfile(this.mninonlinear_dir, "CE_on_T1w.nii.gz");  % mm voxels
            % assert(isfile(g), stackstr() + ": " + g + "not found")
        end
        function g = get.CE_ic(this)
            % assert(isfile(this.CE_fqfn));
            g = mlfourd.ImagingContext2(this.CE_fqfn);
        end
        function g = get.ciftify_subject_fmri_log_fqfn(this)
            pth = fileparts(this.task_dtseries_fqfn);
            g = fullfile(pth, "ciftify_subject_fmri.log");
        end
        function g = get.datashare_dir(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "AnalyticSignalGBM", "GBM_datashare");
            assert(isfolder(g))
        end
        function g = get.edema_fqfn(this)
            g = fullfile(this.mninonlinear_dir, "edema_on_T1w.nii.gz");  % mm voxels
            % assert(isfile(g), stackstr() + ": " + g + "not found")
        end
        function g = get.edema_ic(this)
            if ~isfile(this.edema_fqfn)
                edema = this.WT_ic & ~this.CE_ic;
                edema.fileprefix = "edema_on_T1w";
                save(edema);
            end
            % assert(isfile(this.edema_fqfn), stackstr() + ": " + this.edema_fqfn + "not found")
            g = mlfourd.ImagingContext2(this.edema_fqfn);
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
        function g = get.out_dir(this)
            if ~isemptytext(this.out_dir_)
                g = this.out_dir_;
                return
            end

            % /home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalGBM/analytic_signal/matlabout
            g = fullfile(getenv("SINGULARITY_HOME"), "AnalyticSignalGBM", "analytic_signal", "matlabout");
            assert(isfolder(g))
            this.out_dir_ = g;
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
            % disp(this.task_dir)
            mg = mglob(fullfile(this.task_dir, sprintf("%s*_Atlas_s0.dtseries.nii", this.task)));
            assert(~isemptytext(mg), stackstr())
            g = mg(1);
        end
        function g = get.task_niigz_fqfn(this)
            mg = mglob(fullfile(this.task_dir, sprintf("%s*.nii.gz", this.task)));
            assert(~isemptytext(mg), stackstr())
            g = mg(1);
        end
        function g = get.task_ref_niigz_fqfn(this)
            mg = mglob(fullfile(this.task_dir, sprintf("%s*.nii.gz", this.task)));
            assert(~isemptytext(mg), stackstr())
            g = mg(1);
        end
        function g = get.task_ref_dscalar_fqfn(this)
            % disp(this.task_dir)
            mg = mglob(fullfile(this.task_dir, sprintf("%s*_Atlas_s0_avgt.dtseries.nii", this.task)));
            if isempty(mg)
                g = this.ihcp_.cifti.average_times(this.task_dtseries_fqfn);
            end
        end
        function g = get.thickness_dscalar_fqfn(this)
            g = fullfile(this.mninonlinear_dir, "fsaverage_LR32k", this.sub + ".thickness.32k_fs_LR.dscalar.nii");
            assert(isfile(g), stackstr())
        end
        function g = get.t1w_fqfn(this)
            g = fullfile(this.mninonlinear_dir, "T1w.nii.gz");  % mm voxels
            assert(isfile(g), stackstr() + ": " + g + "not found")
        end
        function g = get.tr(this)
            try
                g = this.json.tr;
            catch ME
                handwarning(ME)
                g = 2.4049;  % most likely tr for GBM
            end
        end
        function g = get.wmparc_fqfn(this)
            g = fullfile(this.mninonlinear_dir, "wmparc.nii.gz");  % mm voxels
            assert(isfile(g), stackstr() + ": " + g + "not found")
        end
        function g = get.WT_fqfn(this)
            g = fullfile(this.mninonlinear_dir, "WT_on_T1w.nii.gz");  % mm voxels
            % assert(isfile(g), stackstr() + ": " + g + "not found")
        end
        function g = get.WT_ic(this)
            % assert(isfile(this.WT_fqfn));
            g = mlfourd.ImagingContext2(this.WT_fqfn);
        end
    end

    methods
        function this = GBMCiftifyData2(ihcp, out_dir)
            arguments
                ihcp mlraut.HCP
                out_dir {mustBeTextScalar} = ""
            end

            this = this@mlraut.CohortData(ihcp);

            this.out_dir_ = out_dir;

            if ~isfile(this.json_fqfn)
                this.build_gbm_json();
            end
        end

        function build_edema(this)
            if isfile(fullfile(this.mninonlinear_dir, 'edema_on_T1w.nii.gz'))
                return
            end
            ce = mlfourd.ImagingContext2(fullfile(this.mninonlinear_dir, 'CE_on_T1w.nii.gz'));
            wt = mlfourd.ImagingContext2(fullfile(this.mninonlinear_dir, 'WT_on_T1w.nii.gz'));
            edema = wt & ~ce;
            edema.fileprefix = "edema_on_T1w";
            save(edema);
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

            try
                j.tr = TR;
                j.num_frames_ori = numTRs;
                j.num_frames_to_trim = this.num_frames_to_trim;
                try
                    hemi = lower(this.find_table_value("hemi"));
                catch 
                    hemi = "unknown";
                end
                try
                    loc = lower(this.find_table_value("brain_location"));
                catch
                    loc = "unknown";
                end
                j.location = sprintf("%s %s", hemi, loc);
            catch ME
                handwarning(ME)
            end
            if ~isfile(this.json_fqfn)
                jsonwrite(j, this.json_fqfn);
            else
                warning("mlraut:IOWarning", ...
                    "%s aborted because %s already exists.", stackstr(), this.json_fqfn);
            end
        end

        function val = find_table_value(this, var_name)
            if any(contains(this.table_gbm.Properties.VariableNames, var_name, IgnoreCase=true))
                this_i3cr = strrep(this.ihcp_.current_subject, 'sub-', '');
                this_i3cr = strrep(this_i3cr, '/', '');
                if contains(this.map_rt_i3cr.keys, this_i3cr)
                    this_i3cr = this.map_rt_i3cr(this_i3cr);  % rt -> i3cr id
                end
                selected = strcmp(this.table_gbm.I3CRID, this_i3cr);
                val = string(this.table_gbm{selected, var_name});
                % disp(this_i3cr)
                % disp(var_name)
                % disp(val)
                assert(isscalar(val))
            end
        end

        function g = surf_gii_fqfn(this, hemis)
            arguments
                this mlraut.GBMCiftifyData2
                hemis {mustBeTextScalar} = "L"
            end

            if startsWith(hemis, "L", IgnoreCase=true)
                hemis = "L";
            elseif startsWith(hemis, "R", IgnoreCase=true)
                hemis = "R";
            end

            mg = mglob(fullfile(this.mninonlinear_dir, "fsaverage_LR32k", "*."+hemis+".sphere_MSMall.32k_fs_LR.surf.gii"));
            % 995174.L.midthickness_MSMAll.32k_fs_LR.surf.gii
            if isemptytext(mg)
                mg = mglob(fullfile(this.mninonlinear_dir, "fsaverage_LR32k", "*."+hemis+".sphere.32k_fs_LR.surf.gii"));
                % 995174.L.midthickness.32k_fs_LR.surf.gii
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

        function in_on_T1w = apply_ciftify_warps(in, sid, opts)
            %% warps any NIfTI in the space of T1w/T1w.nii.gz to space of MNINonLinear

            arguments
                in {mustBeFile}
                sid {mustBeTextScalar}
                opts.mni_nonlin_dir = ""
                opts.workdir {mustBeFolder} = "/home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify"
            end            
            [~, in_bname] = fileparts(in);
            [~, in_bname] = fileparts(in_bname);
            if isemptytext(opts.mni_nonlin_dir)
                opts.mni_nonlin_dir = strrep(fileparts(in), "T1w", "MNINonLinear");
            end
            native_T1w = fullfile(fileparts(in), "T1w.nii.gz");
            mni_nonlin_T1w = fullfile(opts.mni_nonlin_dir, "T1w.nii.gz");
            xfms_dir = fullfile(opts.mni_nonlin_dir, "xfms");

            % needed for sensible return
            in_on_T1w = "";

            pwd0 = pushd(opts.mni_nonlin_dir);
            try

                %% transform from space of 001 to space of T1w/T1w

                dockerout = fileparts(opts.workdir);
                zzo = fullfile(fileparts(in), "001.nii.gz");
                zzo_on_native_T1w = fullfile(fileparts(in), "001_on_T1w.nii.gz");
                zzo_on_native_T1w_mat = fullfile(fileparts(in), "001_on_T1w.mat");
                if ~isfile(zzo_on_native_T1w_mat)

                    % 001.mgz -> 001.nii.gz
                    mgz = fullfile(dockerout, "freesurfer", "sub-" + sid, "mri", "orig", "001.mgz");
                    if ~isfile(mgz)
                        return
                    end
                    mysystem(sprintf("mri_convert %s %s", mgz, zzo));

                    % 001.nii.gz -> T1w/T1w.nii.gz
                    flirt_opts = "-bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 6 -interp trilinear";
                    mysystem(sprintf("flirt -in %s -ref %s -out %s -omat %s %s", zzo, native_T1w, zzo_on_native_T1w, zzo_on_native_T1w_mat, flirt_opts));
                end
                intermed = fullfile(fileparts(in), in_bname + "_on_intermed.nii.gz");
                mysystem(sprintf("flirt -in %s -applyxfm -init %s -out %s -ref %s -interp nearestneighbour", in, zzo_on_native_T1w_mat, intermed, native_T1w));

                %% use ciftify's warps

                premat = fullfile(xfms_dir, "T1w2StandardLinear.mat");
                warp1 = fullfile(xfms_dir, "T1w2Standard_warp_noaffine.nii.gz");
                warp2 = fullfile(xfms_dir, "combined_warp.nii.gz");
                mysystem(sprintf("convertwarp --ref=%s --warp1=%s --premat=%s --out=%s", ...
                    mni_nonlin_T1w, warp1, premat, warp2));

                in_on_T1w = fullfile(opts.mni_nonlin_dir, strrep(in_bname, "seg_", "") + "_on_T1w.nii.gz");
                mysystem(sprintf("applywarp --ref=%s --in=%s --warp=%s --interp=nn --out=%s", ...
                    mni_nonlin_T1w, intermed, warp2, in_on_T1w));

            catch ME
                handwarning(ME)
            end
            popd(pwd0);
        end

        function apply_ciftify_warps_to_segs(workdir)
            arguments
                workdir {mustBeText} = ""
            end

            if ~isfolder(workdir)
                workdir = "/home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify";
            end
            
            native_T1ws = mglob(fullfile(workdir, "sub-*", "T1w", "T1w.nii.gz"));
            native_T1w_dirs = fileparts(native_T1ws);

            parfor (idx = 1:length(native_T1w_dirs), 20)
                ndir = native_T1w_dirs(idx);
                sid = mlraut.GBMCiftifyData2.path_to_sid(ndir);
                for seg = ["seg_CE.nii", "seg_WT.nii", "seg_NCT.nii", "mask_WT.nii"]
                    try
                        fqfn = fullfile(ndir, seg);
                        if ~isfile(fqfn)
                            continue
                        end
                        in_on_T1w = mlraut.GBMCiftifyData2.apply_ciftify_warps(fqfn, sid, workdir=workdir);
                        fprintf("====== %s generated %s =====\n", stackstr(), in_on_T1w);
                    catch ME
                        handwarning(ME)
                    end
                end
            end
        end

        function sid = path_to_sid(pth)
            parts = strsplit(pth, filesep);
            sub_fold = parts(contains(parts, "sub-"));
            sid = strrep(sub_fold, "sub-", "");
        end

        function find_native_segs()
            preproc = "/data/nil-bluearc/shimony/bidhan/rsFC_PreProc";

            native_T1ws = mglob(fullfile("/home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify/sub-*/T1w/T1w.nii.gz"));
            native_T1w_dirs = fileparts(native_T1ws);

            for ndir = native_T1w_dirs
                [~,sfold] = fileparts(fileparts(ndir));
                sid = strrep(sfold, "sub-", "");
                atlas_dir = mglob(fullfile(preproc, "*"+sid+"*", "atlas"));
                if isempty(atlas_dir) 
                    continue
                end

                for seg = ["seg_CE.nii", "seg_WT.nii", "seg_NCT.nii", "mask_WT.nii"]
                    try
                        fqfn = fullfile(atlas_dir, seg);
                        if isfile(fqfn)
                            copyfile(fqfn, ndir)
                        end
                    catch ME
                        handwarning(ME)
                    end
                end
            end
        end

        function flirt_tumor_segs()
            %% use after prepare_tumor_segs

            srcdir = "/home/usr/jjlee/Singularity/AnalyticSignalGBM/analytic_signal/tmp";
            workdir = "/home/usr/jjlee/mnt/Workaround/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify";  
            % ciftify on NIL ZFS is unwritable because of docker-ciftify

            SUBS = [mglob(fullfile(srcdir, "I3CR*")), mglob(fullfile(srcdir, "RT*"))];  % length = 341
            for s = 1:length(SUBS)
                pth_ = char(SUBS(s));
                if strcmp(pth_(end), filesep)
                    [~,fold_] = fileparts(pth_(1:end-1));
                else
                    [~,fold_] = fileparts(pth);
                end
                SUBS(s) = string(fold_);  % trailing folder, no fileseps                
            end
            SEGS = ["CE", "WT"];

            parfor (idx = 1:length(SUBS), 20)
            % for idx = 64:64  % length(SUBS)
                try

                    %% work with Kiyun's GBM imaging in srcdir; continue if data intermediates aren't ready

                    pwd1 = pushd(fullfile(srcdir, SUBS(idx), "atlas"));
    
                    mni_nonlin_dir = fullfile(workdir, "sub-"+SUBS(idx), "MNINonLinear");
                    ref_T1w = fullfile(mni_nonlin_dir, "T1w.nii.gz");
                    if ~isfile(ref_T1w)  % ciftify is not available
                        continue; 
                    end

                    %wt = fullfile(mni_nonlin_dir, "WT_on_T1w.nii.gz");
                    %if isfile(wt)  % already done
                    %    continue;
                    %end

                    %% mpr1

                    native_mpr1 = mglob("*_orient-rpi_mpr1.nii.gz");
                    native_mpr1 = native_mpr1(contains(native_mpr1, SUBS(idx)));
                    if isempty(native_mpr1)  % missing mpr1 from srcdir (viz., GBM repository)
                        continue; 
                    end
                    native_mpr1 = native_mpr1(1);

                    %% mskt

                    native_mskt = mglob("*_orient-rpi_mpr1_mskt.nii.gz");
                    native_mskt = native_mskt(contains(native_mskt, SUBS(idx)));

                    %% rigidly transform native mpr1 to ciftify's native T1w

                    intermed_dir = fullfile(workdir, "sub-"+SUBS(idx), "T1w");
                    intermed_T1w = fullfile(intermed_dir, "T1w_brain.nii.gz");

                    % flirt intermediate T1w brain as r.b.
                    intermed_out = "intermed_T1w_on_mpr.nii.gz";
                    intermed_omat = "intermed_T1w_on_mpr.mat";
                    opts = "-bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 6 -interp trilinear";
                    if ~isempty(native_mskt) && ~isemptytext(native_mskt)
                        native_mskt = native_mskt(1);
                        mysystem(sprintf("flirt -in %s -ref %s -refweight %s -out %s -omat %s %s", intermed_T1w, native_mpr1, native_mskt, intermed_out, intermed_omat, opts));
                    else 
                        % missing mskt from srcdir (viz., GBM repository)
                        mysystem(sprintf("flirt -in %s -ref %s -out %s -omat %s %s", intermed_T1w, native_mpr1, intermed_out, intermed_omat, opts));
                    end

                    % invert r.b. transform from intermediate
                    intermed_omat2 = "mpr_on_T1w_intermed.mat";
                    mysystem(sprintf("convert_xfm -omat %s -inverse %s", intermed_omat2, intermed_omat));

                    %% use ciftify's warps

                    xfms_dir = fullfile(mni_nonlin_dir, "xfms");

                    premat = fullfile(xfms_dir, "T1w2StandardLinear.mat");
                    warp1 = fullfile(xfms_dir, "T1w2Standard_warp_noaffine.nii.gz");
                    warpout = "combined_warp.nii.gz";
                    mysystem(sprintf("convertwarp --ref=%s --warp1=%s --premat=%s --out=%s", ...
                        ref_T1w, warp1, premat, warpout));

                    %% transform to intermediate; warp to MNI

                    mpr_on_intermed = "mpr_on_intermed_T1w.nii.gz";
                    mysystem(sprintf("flirt -in %s -applyxfm -init %s -out %s -ref %s", native_mpr1, intermed_omat2, mpr_on_intermed, intermed_T1w));
                    mpr_on_T1w = strrep(native_mpr1, "mpr1.nii.gz", "mpr1_on_T1w.nii.gz");
                    mysystem(sprintf("applywarp --ref=%s --in=%s --warp=%s --interp=nn --out=%s", ...
                        ref_T1w, mpr_on_intermed, warpout, mpr_on_T1w));

                    for sidx = 1:length(SEGS)
                        try
                            
                            % ensure segmentation is viable
                            seg_ori = mglob(sprintf("seg_%s_on_orient-rpi_mpr1.nii.gz", SEGS(sidx)));
                            if ~isfile(seg_ori); continue; end

                            % apply r.b. transform to intermediate
                            seg_on_intermed = SEGS(sidx) + "_on_intermed_T1w.nii.gz";
                            opts1 = "-paddingsize 0.0 -interp nearestneighbour";
                            mysystem(sprintf("flirt -in %s -applyxfm -init %s -out %s %s -ref %s", seg_ori, intermed_omat2, seg_on_intermed, opts1, intermed_T1w));

                            % apply warp:  intermediate T1w => T1w ref
                            seg_on_T1w = fullfile(mni_nonlin_dir, SEGS(sidx) + "_on_T1w.nii.gz");
                            mysystem(sprintf("applywarp --ref=%s --in=%s --warp=%s --interp=nn --out=%s", ...
                                ref_T1w, seg_on_intermed, warpout, seg_on_T1w));

                            % set zero all segmentation values < 1
                            ic = mlfourd.ImagingContext2(seg_on_T1w);
                            ic = ic.thresh(1-eps('single'));
                            ic = ic.binarized();
                            ic.fqfilename = seg_on_T1w;
                            ic.save();
                            fprintf("=============== %s completed %s ===============\n", stackstr(), ic.fqfilename)
                        catch ME
                            handwarning(ME)
                        end
                    end
                    popd(pwd1);
                catch ME
                    handwarning(ME)
                end
            end 
        end

        function ce = migrate_ce(src_atl_dir, targ_atl_dir, opts)
            %% contrast enhancing

            arguments
                src_atl_dir {mustBeFolder}
                targ_atl_dir {mustBeTextScalar}
                opts.subid {mustBeTextScalar} = "*"
            end

            ast = fullfile(src_atl_dir, "seg_CE_on_mpr1.4dfp.*");
            if ~isempty(mglob(ast))
                ensuredir(targ_atl_dir);
                copyfile(ast, targ_atl_dir);
                pwd0 = pushd(targ_atl_dir);
                fp = mybasename(strrep(ast, ".*", ".hdr"));
                mysystem(sprintf("nifti_4dfp -n %s %s", fp, fp));
                mysystem(sprintf("gzip %s.nii", fp));
                mlpipeline.Bids.afni_3dresample(fp+".nii.gz");
                ce = mlfourd.ImagingContext2(fullfile(targ_atl_dir, "seg_CE_on_orient-rpi_mpr1.nii.gz"));
                delete("*.4dfp.*");
                popd(pwd0);
                return
            end
            
            mask_ast = strrep(ast, "seg", "mask");
            if ~isempty(mglob(mask_ast))

                mlraut.GBMCiftifyData2.migrate_mskt(src_atl_dir, targ_atl_dir, subid=opts.subid);

                ensuredir(targ_atl_dir);
                copyfile(mask_ast, targ_atl_dir);
                pwd0 = pushd(targ_atl_dir);
                fp = mybasename(strrep(mask_ast, ".*", ".hdr"));
                mysystem(sprintf("nifti_4dfp -n %s %s", fp, fp));
                mysystem(sprintf("gzip %s.nii", fp));
                mlpipeline.Bids.afni_3dresample(fp+".nii.gz");
                ce = mlfourd.ImagingContext2(fullfile(targ_atl_dir, "mask_CE_on_orient-rpi_mpr1.nii.gz"));
                ce = ce.binarized();
                ce = ~ce;
                msk = mlfourd.ImagingContext2(fullfile(targ_atl_dir, opts.subid+"_orient-rpi_mpr1_mskt.nii.gz"));
                msk = msk.thresh(msk.dipmax/2);
                msk = msk.binarized();
                ce = ce .* msk;
                ce.fileprefix = "seg_CE_on_orient-rpi_mpr1";
                ce.save();
                delete("*.4dfp.*");
                popd(pwd0);
                return
            end
            
            fprintf("%s: could not find any %s\n", stackstr(), ast);
            error("mlraut:FileNotFound", stackstr());
        end

        function mpr = migrate_mpr(src_atl_dir, targ_atl_dir, opts)
            arguments
                src_atl_dir {mustBeFolder}
                targ_atl_dir {mustBeTextScalar}
                opts.subid {mustBeTextScalar} = "*"
            end

            ast = fullfile(src_atl_dir, opts.subid+"_mpr1.4dfp.*");
            if ~isempty(mglob(ast))
                ensuredir(targ_atl_dir);
                copyfile(ast, targ_atl_dir);
                pwd0 = pushd(targ_atl_dir);
                fp = mybasename(strrep(ast, ".*", ".hdr"));
                mysystem(sprintf("nifti_4dfp -n %s %s", fp, fp));
                mysystem(sprintf("gzip %s.nii", fp));
                mlpipeline.Bids.afni_3dresample(fp+".nii.gz");
                mpr = mlfourd.ImagingContext2(fullfile(targ_atl_dir, strrep(fp, "mpr1", "orient-rpi_mpr1")+".nii.gz"));
                delete("*.4dfp.*");
                popd(pwd0);
                return
            end
            
            fprintf("%s: could not find any %s\n", stackstr(), ast);
            error("mlraut:FileNotFound", stackstr());
        end

        function mskt = migrate_mskt(src_atl_dir, targ_atl_dir, opts)
            arguments
                src_atl_dir {mustBeFolder}
                targ_atl_dir {mustBeTextScalar}
                opts.subid {mustBeTextScalar} = "*"
            end

            ast = fullfile(src_atl_dir, opts.subid+"_mpr1_mskt.4dfp.*");
            if ~isempty(mglob(ast))
                ensuredir(targ_atl_dir);
                copyfile(ast, targ_atl_dir);
                pwd0 = pushd(targ_atl_dir);
                fp = mybasename(strrep(ast, ".*", ".hdr"));
                mysystem(sprintf("nifti_4dfp -n %s %s", fp, fp));
                mysystem(sprintf("gzip %s.nii", fp));
                mlpipeline.Bids.afni_3dresample(fp+".nii.gz");
                mskt = mlfourd.ImagingContext2(fullfile(targ_atl_dir, strrep(fp, "mpr1_mskt", "orient-rpi_mpr1_mskt")+".nii.gz"));
                delete("*.4dfp.*");
                popd(pwd0);
                return
            end
            
            fprintf("%s: could not find any %s\n", stackstr(), ast);
            error("mlraut:FileNotFound", stackstr());
        end

        function wt = migrate_wt(src_atl_dir, targ_atl_dir, opts)
            %% whole tumor

            arguments
                src_atl_dir {mustBeFolder}
                targ_atl_dir {mustBeTextScalar}
                opts.subid {mustBeTextScalar} = "*"
            end

            ast = fullfile(src_atl_dir, "seg_WT_on_mpr1.4dfp.*");
            if ~isempty(mglob(ast))
                ensuredir(targ_atl_dir);
                copyfile(ast, targ_atl_dir);
                pwd0 = pushd(targ_atl_dir);
                fp = mybasename(strrep(ast, ".*", ".hdr"));
                mysystem(sprintf("nifti_4dfp -n %s %s", fp, fp));
                mysystem(sprintf("gzip %s.nii", fp));
                mlpipeline.Bids.afni_3dresample(fp+".nii.gz");
                wt = mlfourd.ImagingContext2(fullfile(targ_atl_dir, "seg_WT_on_orient-rpi_mpr1.nii.gz"));
                delete("*.4dfp.*");
                popd(pwd0);
                return
            end
            
            mask_ast = strrep(ast, "seg", "mask");
            if ~isempty(mglob(mask_ast))

                mlraut.GBMCiftifyData2.migrate_mskt(src_atl_dir, targ_atl_dir, subid=opts.subid);

                ensuredir(targ_atl_dir);
                copyfile(mask_ast, targ_atl_dir);
                pwd0 = pushd(targ_atl_dir);
                fp = mybasename(strrep(mask_ast, ".*", ".hdr"));
                mysystem(sprintf("nifti_4dfp -n %s %s", fp, fp));
                mysystem(sprintf("gzip %s.nii", fp));
                mlpipeline.Bids.afni_3dresample(fp+".nii.gz");
                wt = mlfourd.ImagingContext2(fullfile(targ_atl_dir, "mask_WT_on_orient-rpi_mpr1.nii.gz"));
                wt = wt.binarized();
                wt = ~wt;
                msk = mlfourd.ImagingContext2(fullfile(targ_atl_dir, opts.subid+"_orient-rpi_mpr1_mskt.nii.gz"));
                msk = msk.thresh(msk.dipmax/2);
                msk = msk.binarized();
                wt = wt .* msk;
                wt.fileprefix = "seg_WT_on_orient-rpi_mpr1";
                wt.save();
                delete("*.4dfp.*");
                popd(pwd0);
                return
            end
            
            fprintf("%s: could not find any %s\n", stackstr(), ast);
            error("mlraut:FileNotFound", stackstr());
        end

        function prepare_tumor_segs(src_dir, targ_dir)
            %% use before flirt_tumor_segs;
            %  build tumor segs on neuroimage machines, AnalyticSignalGBM/analytic_signal/tmp

            arguments
                src_dir {mustBeFolder} = "/data/nil-bluearc/shimony/jjlee/Kiyun/rsFC_PreProc"
                targ_dir {mustBeFolder} = "/home/usr/jjlee/Singularity/AnalyticSignalGBM/analytic_signal/tmp"
            end

            import mlraut.GBMCiftifyData2.migrate_ce
            import mlraut.GBMCiftifyData2.migrate_wt
            import mlraut.GBMCiftifyData2.migrate_mpr
            import mlraut.GBMCiftifyData2.migrate_mskt
            
            % rsFC_PreProc has mixed I3CR and RT ids; set pwd there
            pwd0 = pushd(src_dir);
            sub_folds = [asrow(mglob("I3CR*")), asrow(mglob("RT*"))];
            sub_folds = sub_folds(isfolder(sub_folds));
            sub_folds = strrep(sub_folds, filesep, "");

            % parfor (sidx = 1:length(sub_folds), 16)
            for sidx = 1:length(sub_folds)
                subid = sub_folds(sidx);
                fprintf("%s: working on %s\n", stackstr(), subid);
                try                        
                    src_atl_dir = fullfile(src_dir, subid, "atlas");
                    if contains(subid, "_MR_Pre")
                        subid = extractBefore(subid, "_MR_Pre");
                    end
                    targ_atl_dir = fullfile(targ_dir, subid, "atlas");

                    % skip if missing original subject atlas
                    if ~isfolder(src_atl_dir); continue; end

                    % skip if there are existing results
                    if isfile(fullfile(targ_atl_dir, "seg_WT_on_orient-rpi_mpr1.nii.gz")) && ...
                            isfile(fullfile(targ_atl_dir, subid+"_orient-rpi_mpr1.nii.gz"))
                        continue; 
                    end

                    % arrange targ_atl_dir
                    mpr = migrate_mpr(src_atl_dir, targ_atl_dir, subid=subid);
                    wt = migrate_wt(src_atl_dir, targ_atl_dir, subid=subid);
                    ce = migrate_ce(src_atl_dir, targ_atl_dir, subid=subid);
                    mskt = migrate_mskt(src_atl_dir, targ_atl_dir, subid=subid);
                catch ME
                    handwarning(ME)
                end
            end

            popd(pwd0);

        end

        function prepare_mskts(src_dir, targ_dir)
            %% use before flirt_tumor_segs;
            %  build tumor segs on neuroimage machines, AnalyticSignalGBM/analytic_signal/tmp

            arguments
                src_dir {mustBeFolder} = "/data/nil-bluearc/shimony/jjlee/Kiyun/rsFC_PreProc"
                targ_dir {mustBeFolder} = "/home/usr/jjlee/Singularity/AnalyticSignalGBM/analytic_signal/tmp"
            end

            import mlraut.GBMCiftifyData2.migrate_mskt
            
            % rsFC_PreProc has mixed I3CR and RT ids; set pwd there
            pwd0 = pushd(src_dir);
            sub_folds = [asrow(mglob("I3CR*")), asrow(mglob("RT*"))];
            sub_folds = sub_folds(isfolder(sub_folds));
            sub_folds = strrep(sub_folds, filesep, "");

            parfor (sidx = 1:length(sub_folds), 20)
            % for sidx = 1:length(sub_folds)
                subid = sub_folds(sidx);
                fprintf("%s: working on %s\n", stackstr(), subid);
                try                        
                    src_atl_dir = fullfile(src_dir, subid, "atlas");
                    if contains(subid, "_MR_Pre")
                        subid = extractBefore(subid, "_MR_Pre");
                    end
                    targ_atl_dir = fullfile(targ_dir, subid, "atlas");

                    % skip if missing original subject atlas
                    if ~isfolder(src_atl_dir); continue; end

                    % arrange targ_atl_dir
                    migrate_mskt(src_atl_dir, targ_atl_dir, subid=subid);
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
        map_rt_i3cr_   % containers.Map supports text -> char, rt -> i3cr, i3cr -> i3cr
        out_dir_
        table_excluded_  % chars
        table_gbm_  % strings
        table_rt_i3cr_  % chars
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
