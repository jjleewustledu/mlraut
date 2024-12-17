classdef Ciftify < handle & mlsystem.IHandle
    %% See also:
    %  https://unfmontreal.github.io/Dcm2Bids/docs/tutorial/first-steps/
    %  https://edickie.github.io/ciftify/#/01_installation
    %  https://zenodo.org/record/3369937#.ZDMelezMIeY
    %  https://www.sciencedirect.com/science/article/pii/S1053811919303714?via%3Dihub#appsec1
    %  
    %  Created 20-Jun-2023 14:03:04 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 9.14.0.2286388 (R2023a) Update 3 for MACI64.  Copyright 2023 John J. Lee.
    
    properties
        surferdir = '/data/nil-bluearc/shimony/bidhan/FSCompleted/FSout'
        scansdir = '/data/nil-bluearc/shimony/bidhan/rsFC_PreProc'
        scansdir2 = '/data/nil-bluearc/shimony/Park/Proj_DTI'
        scansdir3 = '/data/nil-bluearc/shimony/ADaniel/Tumor2/rsfMRITumor'
    end

    properties (Constant)
        HOME = '/vgpool02/data2/jjlee'
    end

    properties (Dependent)
        freesurfer
        fsaverage
        not_missing_surfer % e.g., I3CR0000_NoT2
        rawdatadir
        tmp_dcm2bids
        workdir
    end

    methods %% GET
        function g = get.freesurfer(this)
            g = fullfile(this.workdir, 'freesurfer');
        end
        function g = get.fsaverage(this)
            g = fullfile(this.surferdir, 'fsaverage');
        end
        function g = get.not_missing_surfer(this)
            if ~isempty(this.not_missing_surfer_)
                g = this.not_missing_surfer_;
                return
            end

            ld = load(fullfile(this.freesurfer, "not_missing_surfer.mat"));
            this.not_missing_surfer_ = ...
                cellfun(@(x) strip(x, "right", filesep), ld.not_missing_surfer, UniformOutput=false);
            g = this.not_missing_surfer_;
        end
        function g = get.rawdatadir(this)
            g = fullfile(this.workdir, 'rawdata');
        end
        function g = get.tmp_dcm2bids(this)
            g = fullfile(this.workdir, "tmp_dcm2bids");
        end
        function g = get.workdir(~)
           g = fullfile(getenv('SINGULARITY_HOME'), 'AnalyticSignalHCP', 'analytic_signal');
        end
    end

    methods
        function this = Ciftify(varargin)
        end

        function prep_filesystem_for_fmriprep(this)
            arguments
                this mlraut.Ciftify
            end
            wrkdir = this.workdir;

            pwd0 = pushd(wrkdir);

            globbed = glob('derivatives-*');
            globbed = cellfun(@(x) strip(x, "right", filesep), globbed, UniformOutput=false);
            for idx = 1:length(globbed)
                try
                    % aufbau bids assets
                    s = sprintf("rsync -a %s %s", ".bidsignore", globbed{idx});
                    s1 = sprintf("rsync -a %s %s", "dataset_description.json", globbed{idx});
                    system(s);
                    system(s1);
                catch ME
                    handwarning(ME)
                    fprintf("%s: err while working with %s\n ", stackstr(), globbed{idx})
                end
            end

            popd(pwd0);
        end
        function prep_filesystem_for_surfer(this)
            arguments
                this mlraut.Ciftify
            end
            wrkdir = this.workdir;
            rawdir = this.rawdatadir;
            not_missing_surfer1 = this.not_missing_surfer;

            pwd0 = pushd(rawdir);

            globbed = glob('{I3CR*,RT*}');
            globbed = cellfun(@(x) strip(x, "right", filesep), globbed, UniformOutput=false);
            for idx = 1:length(globbed)
                try
                    % check for missing data
                    gfold = globbed{idx};
                    if ~isfolder(fullfile(rawdir, gfold, "SCANS")) || ...
                            ~any(contains(not_missing_surfer1, gfold))
                        continue
                    end

                    % aufbau filesystem tree
                    re = regexp(gfold, "(?<subId>[ICRT0-9]+)(|_\w+)", "names");
                    subId = re.subId;
                    source = fullfile(this.freesurfer, gfold);
                    destination = fullfile(wrkdir, "dockerout-"+subId, "freesurfer", "sub-"+subId);
                    ensuredir(destination);
                    s = sprintf("rsync -raL %s/* %s", source, destination);
                    system(s);
                catch ME
                    handwarning(ME)
                    fprintf("%s: err while working with %s\n ", stackstr(), globbed{idx})
                end
            end

            popd(pwd0);
        end
        function par_dcm2bids(this, exclusions)
            %% 1st activate anaconda: "conda activate dcm2bids"

            arguments
                this mlraut.Ciftify
                exclusions {mustBeText} = "default" % = ["I3CR0002", "I3CR1843", "I3CR0006"]
            end
            wrkdir = this.workdir;
            rawdir = this.rawdatadir;
            code_json = fullfile(wrkdir, "code", "dcm2bids_config.json");
            not_missing_surfer1 = this.not_missing_surfer;

            pwd0 = pushd(rawdir);

            globbed = glob('{I3CR,RT}*');
            globbed = cellfun(@(x) strip(x, "right", filesep), globbed, UniformOutput=false);
            globbed = globbed(~contains(globbed, exclusions));
            parfor idx = 3:length(globbed)
                try
                    % check for missing data
                    gfold = globbed{idx};
                    if ~isfolder(fullfile(rawdir, gfold, "SCANS")) || ...
                            ~any(contains(not_missing_surfer1, gfold))
                        continue
                    end
                    
                    % aufbau filesystem tree
                    re = regexp(gfold, "(?<subId>[ICRT0-9]+)(|_\w+)", "names");
                    subId = re.subId;
                    derivdir = fullfile(wrkdir, "derivatives-"+subId);
                    ensuredir(derivdir);
                    outdir = fullfile(wrkdir, "dockerout-"+subId);
                    ensuredir(outdir);

                    % aufbau bids
                    pwd1 = pushd(derivdir);
                    s = sprintf("dcm2bids -d %s/%s/SCANS -p sub-%s -s 1 -c %s -o %s --forceDcm2niix", ...
                        rawdir, gfold, subId, code_json, derivdir);
                    s1 = sprintf("mv tmp_dcm2bids/sub-* %s", this.tmp_dcm2bids);
                    mysystem(s);
                    mysystem(s1);

                    popd(pwd1);
                catch ME
                    handwarning(ME)
                    fprintf("%s: err while working with %s\n ", stackstr(), globbed{idx})
                    if contains(ME.message, "dcm2bids")
                        fprintf("%s: ensure execution of 'conda activate dcm2bids'\n", stackstr())
                    end
                end
            end

            popd(pwd0);
        end
    end

    methods (Static)
        function propcluster()
            c = parcluster;
            c.AdditionalProperties.EmailAddress = '';
            c.AdditionalProperties.EnableDebug = 1;
            c.AdditionalProperties.GpusPerNode = 0;
            c.AdditionalProperties.MemUsage = '4000'; % in MB; GBM T1w & BOLD data ~ 131 MB; Matlab ~ 3.5 GB
            c.AdditionalProperties.Node = '';
            c.AdditionalProperties.Partition = '';
            c.AdditionalProperties.WallTime = '12:00:00'; % 24 h
            c.saveProfile
            disp(c.AdditionalProperties)
        end
        function [j,c] = parcluster(n_cpus, opts)
            %% #PARCLUSTER
            %  See also https://sites.wustl.edu/chpc/resources/software/matlab_parallel_server/
            
            arguments
                n_cpus double = 1
                opts.to_skip double = 0
            end

            c = parcluster;
            disp(c.AdditionalProperties)            
            mladni.CHPC3.setenvs();

            ld = load(fullfile(getenv('SINGULARITY_HOME'), 'AnalyticSignalHCP', 'analytic_signal', 'parcluster_globbed.mat'));
            globbed = ld.globbed;
            if 1 == n_cpus
                j = c.batch(@mlraut.Ciftify.fmriprep_ciftify_singularity, 1, {globbed{1}, 1}, ...
                    'CurrentFolder', '.', 'AutoAddClientPath', false);
                return
            end
            for ji = opts.to_skip+1:n_cpus
                j{ji} = c.batch(@mlraut.Ciftify.fmriprep_ciftify_singularity, 1, {globbed{ji}, 1}, ...
                    'CurrentFolder', '.', 'AutoAddClientPath', false);  %#ok<AGROW>
            end
        end
        function parlocal(to_glob, n_cpus)
            arguments
                to_glob {mustBeTextScalar} = pwd % derivatives dir
                n_cpus double = 8
            end
            assert(contains(to_glob, "deriv"), "%s: to_glob->%s", stackstr(), to_glob)

            globbed = glob(to_glob);
            globbed = cellfun(@(x) strip(x, "right", filesep), globbed, UniformOutput=false);
            if 1 == n_cpus
                for gi = 1:length(globbed)
                    try
                        t = mlraut.Ciftify.fmriprep_ciftify_docker(globbed{gi});
                        fprintf("%s: elapsed time %s", stackstr(), t);
                    catch ME
                        handexcept(ME)
                    end
                end
                return
            end
            parpool(n_cpus)
            parfor (gi = 1:length(globbed), n_cpus)
                try
                    t = mlraut.Ciftify.fmriprep_ciftify_docker(globbed{gi});
                    fprintf("%s: elapsed time %s", stackstr(), t);
                catch ME
                    handexcept(ME)
                end
            end
        end
        function t = fmriprep_ciftify_docker(derivdir, n_cpus)
            %% docker run -ti --rm -v ${AS_HOME}/derivatives:/data:ro -v ${AS_HOME}/dockerout:/out -v ${AS_HOME}/license.txt:/fs_license.txt tigrlab/fmriprep_ciftify:v1.3.2-2.3.3 /data /out participant --n_cpus 12 --verbose --debug --fs-license /fs_license.txt 

            arguments
                derivdir {mustBeFolder}
                n_cpus {mustBeScalarOrEmpty} = 1
            end
            %mladni.CHPC3.setenvs();            
            HOME= mlraut.Ciftify.HOME;
            t0 = tic;
            try
                outdir = strrep(derivdir, "derivatives", "dockerout");
                license = fullfile(HOME, "AnalyticSignalOAS", "analytic_signal", "license.txt");
                image = "tigrlab/fmriprep_ciftify:v1.3.2-2.3.3";
                s = sprintf("docker run -ti --rm -v %s:/data:ro -v %s:/out -v %s:/fs_license.txt %s /data /out participant --n_cpus %g --verbose --debug --fs-license /fs_license.txt", ...
                    derivdir, outdir, license, image, n_cpus);
                system(s)
            catch ME
                handwarning(ME)
            end
            t = toc(t0);  
        end
        function t = fmriprep_ciftify_singularity(derivdir, n_cpus)
            %% singularity run --cleanenv /my_images/fmriprep_cifitfy-1.1.2-2.1.0.simg \
            %     path/to/data/dir path/to/output/dir \
            %     participant \
            %     --participant-label label
            %  See also:  mladni.FDG.deepmrseg_apply()
            
            arguments
                derivdir {mustBeTextScalar}
                n_cpus {mustBeScalarOrEmpty} = 1
            end
            mladni.CHPC3.setenvs();            
            t0 = tic;
            try
                derivdir = fullfile(getenv("SINGULARITY_HOME"), "AnalyticSignalHCP", "analytic_signal", mybasename(derivdir));
                outdir = strrep(derivdir, "derivatives", "dockerout");
                license = fullfile(getenv("SINGULARITY_HOME"), "AnalyticSignalHCP", "analytic_signal", "license.txt");
                image = fullfile(getenv("SINGULARITY_HOME"), "AnalyticSignalHCP", "tigrlab_fmriprep_ciftify_v1.3.2-2.3.3-2019-08-16-0a84f7a43b38.simg");
                s = sprintf("singularity run --cleanenv -B %s:/data -B %s:/out -B %s:/fs_license.txt %s /data /out participant --n_cpus %g --fs-license /fs_license.txt", ...
                    derivdir, outdir, license, image, n_cpus);
                system(s)
            catch ME
                handwarning(ME)
            end
            t = toc(t0);  
        end
    end

    %% PROTECTED

    properties (Access = protected)
        not_missing_surfer_
    end

    methods (Access = protected)
        function that = copyElement(this)
            that = copyElement@matlab.mixin.Copyable(this);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
