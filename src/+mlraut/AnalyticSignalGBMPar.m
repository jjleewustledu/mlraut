classdef AnalyticSignalGBMPar < handle & mlraut.AnalyticSignalGBM
    %% line1
    %  line2
    %  
    %  Created 29-Jun-2023 11:22:03 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 9.14.0.2286388 (R2023a) Update 3 for MACI64.  Copyright 2023 John J. Lee.
    

    properties (Constant)
        unrecoverable_ce_wt = [ ...
            "sub-I3CR1149", "sub-I3CR0483", "sub-I3CR1318", "sub-I3CR0632", "sub-I3CR0639", "sub-I3CR0964", "sub-I3CR1026" ...
        ]  % see also Singularity/AnalyticSignaGBM/GBM_datashare/README_JJL.docx
    end

    methods (Static)
        function out = link_kiyun(sub)
            tic
            mlraut.CHPC3.setenvs();
            mnidir = fullfile(getenv("SINGULARITY_HOME"), ...
                "AnalyticSignalGBM", "analytic_signal", "dockerout", "ciftify", sub, "MNINonLinear");
            defectsdir = fullfile(mnidir, "Defects");
            mkdir(defectsdir);
            
            moveExisting(fullfile(mnidir, "*flip-1.nii.gz"), defectsdir);
            moveExisting(fullfile(mnidir, "CE_on_T1w.nii.gz"), defectsdir);
            moveExisting(fullfile(mnidir, "edema_on_T1w.nii.gz"), defectsdir);
            moveExisting(fullfile(mnidir, "WT_on_T1w.nii.gz"), defectsdir);
            try
                cek = fullfile(mnidir, "CE_on_T1w_kiyun.nii.gz");
                ce = fullfile(mnidir, "CE_on_T1w.nii.gz");
                system(sprintf("ln -s %s %s", cek, ce));
            catch ME
                handwarning(ME)
            end
            try
                wtk = fullfile(mnidir, "WT_on_T1w_kiyun.nii.gz");
                wt = fullfile(mnidir, "WT_on_T1w.nii.gz");
                system(sprintf("ln -s %s %s", wtk, wt));
            catch ME
                handwarning(ME)
            end

            out = toc;
        end

        function this = mean_twistor(physio)
            %% see also mlraut.Lee2024.build_mean_for_gbm()

            arguments
                physio {mustBeTextScalar} = 'physio_iFV'
            end

            this = mlraut.AnalyticSignalGBMPar( ...
                subjects={'sub-I3CR0023'}, ...
                do_7T=false, ...
                do_resting=true, ...
                do_task=false, ...
                do_global_signal_regression=false, ...
                do_save=false, ...
                do_save_dynamic=false, ...
                do_save_ciftis=false, ...
                hp_thresh=0.01, ...
                lp_thresh=0.1, ...
                tags="AnalyticSignalGBMPar-mean-twistor");

            mats = asrow(glob(fullfile(this.out_dir, physio, '*', 'sub-*_ses-*AnalyticSignalGBM*.mat')));
            n = length(mats);
            nx = this.num_nodes;

            X_ = zeros(1, nx);
            Y_ = zeros(1, nx);
            Z_ = zeros(1, nx);
            T_ = zeros(1, nx);

            for mat = mats
                % tic
                ld = load(mat{1});
                this_subset = ld.this_subset;
                try
                    psi = this_subset.bold_signal;
                catch ME
                    if strcmp(ME.identifier, 'MATLAB:nonExistentField')
                        psi = this_subset.analytic_signal.*this_subset.physio_signal;  % overly normalized in this_subset
                    else
                        rethrow(ME)
                    end
                end
                phi = this_subset.physio_signal;
                X = mean((psi.*conj(phi) + phi.*conj(psi))/sqrt(2), 1);
                Y = mean((psi.*conj(phi) - phi.*conj(psi))/sqrt(2i), 1);
                Z = mean((psi.*conj(psi) - phi.*conj(phi))/sqrt(2), 1);
                T = mean((psi.*conj(psi) + phi.*conj(phi))/sqrt(2), 1);

                X_ = X_ + X/n;
                Y_ = Y_ + Y/n;
                Z_ = Z_ + Z/n;
                T_ = T_ + T/n;
                % toc
            end

            this.write_ciftis( ...
                X_, sprintf('X_as_sub-all_ses-all_%s', this.tags), ...
                do_save_dynamic=false);
            this.write_ciftis( ...
                Y_, sprintf('Y_as_sub-all_ses-all_%s', this.tags), ...
                do_save_dynamic=false);
            this.write_ciftis( ...
                Z_, sprintf('Z_as_sub-all_ses-all_%s', this.tags), ...
                do_save_dynamic=false);
            this.write_ciftis( ...
                T_, sprintf('T_as_sub-all_ses-all_%s', this.tags), ...
                do_save_dynamic=false);
            this.write_ciftis( ...
                T_+Z_, sprintf('T+Z_as_sub-all_ses-all_%s', this.tags), ...
                do_save_dynamic=false);
        end
        
        function parcall(cores, opts)
            arguments
                cores {mustBeScalarOrEmpty} = 32
                opts.config_hemispheres {mustBeTextScalar} = ""   % "lesionR-iFV" "lesionR-CE" "nolesion" "alllesion" ""
            end

            root_dir = '/home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify';
            out_dir = '/home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalGBM/analytic_signal/matlabout/lesionR-CE';
            ensuredir(out_dir);
            tasks = {'ses-1_task-rest_run-01_desc-preproc', 'ses-1_task-rest_run-02_desc-preproc'};

            g = glob(fullfile(root_dir, 'sub-*'));
            %g = flip(g);
            g = cellfun(@(x) basename(x), g, UniformOutput=false);
            %g = g(1:end);
            leng = length(g);
            parfor (idxg = 1:leng, cores)
                try
                    this = mlraut.AnalyticSignalGBMPar(subjects=g(idxg), ...
                        root_dir=root_dir, out_dir=out_dir, tasks=tasks, ...
                        source_physio='CE');
                    this.config_hemispheres = opts.config_hemispheres; %#ok<PFBNS>
                    call(this);
                catch ME
                    handwarning(ME)
                end
            end
        end

        function parcall_physios(cores)   
            %% left insula, no midline shift

            arguments
                cores {mustBeScalarOrEmpty} = 3
            end

            root_dir = fullfile( ...
                getenv('SINGULARITY_HOME'), ...
                'AnalyticSignalGBM/analytic_signal/dockerout/ciftify');
            cd(root_dir);

            g = convertStringsToChars("sub-" + mlraut.GBMCiftifyData2.SUBS);
            leng = length(g);
            for idxg = 1:1
            % parfor (idxg = 1:leng, cores)
                try
                    tic

                    out_dir = fullfile( ...
                        getenv('SINGULARITY_HOME'), ...
                        'AnalyticSignalGBM/analytic_signal/matlabout/physio_iFV');
                    ensuredir(out_dir);
                    this = mlraut.AnalyticSignalGBM( ...
                        subjects=g(idxg), ...
                        tasks={'ses-1_task-rest_run-all_desc-preproc'}, ...
                        do_save=true, ...
                        do_save_ciftis=true, ...
                        out_dir=out_dir, ...
                        source_physio='iFV-brightest');
                    call(this);

                    out_dir = fullfile( ...
                        getenv('SINGULARITY_HOME'), ...
                        'AnalyticSignalGBM/analytic_signal/matlabout/physio_CE');
                    ensuredir(out_dir);
                    ce = mlfourd.ImagingContext2( ...
                        fullfile(root_dir, g(idxg), 'MNINonLinear', 'CE_on_T1w.nii.gz'));
                    this = mlraut.AnalyticSignalGBM( ...
                        subjects=g(idxg), ...
                        tasks={'ses-1_task-rest_run-all_desc-preproc'}, ...
                        do_save=true, ...
                        do_save_ciftis=true, ...
                        out_dir=out_dir, ...
                        roi=ce);
                    call(this);
                    
                    % out_dir = fullfile( ...
                    %     getenv('SINGULARITY_HOME'), ...
                    %     'AnalyticSignalGBM/analytic_signal/matlabout/physio_WT');
                    % ensuredir(out_dir);
                    % wt = mlfourd.ImagingContext2( ...
                    %     fullfile(root_dir, g(idxg), 'MNINonLinear', 'WT_on_T1w.nii.gz'));
                    % this = mlraut.AnalyticSignalGBM( ...
                    %     subjects=g(idxg), ...
                    %     tasks={'ses-1_task-rest_run-all_desc-preproc'}, ...
                    %     do_save=true, ...
                    %     do_save_ciftis=true, ...
                    %     out_dir=out_dir, ...
                    %     roi=wt);
                    % call(this);

                    toc
                catch ME
                    handwarning(ME)
                end
            end

            % Elapsed time is ___ seconds.
        end

        %% running call on single server or cluster

        function server_call_sub(sub, opts)
            %% for servers

            arguments
                sub {mustBeText} = "sub-I3CR0015"
                opts.source_physio {mustBeText} = "iFV-brightest"
            end
            sub = {convertStringsToChars(sub)};

            root_dir = fullfile( ...
                getenv('SINGULARITY_HOME'), 'AnalyticSignalGBM', 'analytic_signal', 'dockerout', 'ciftify');
            cd(root_dir);
            out_dir = fullfile( ...
                getenv('SINGULARITY_HOME'), 'AnalyticSignalGBM', 'analytic_signal', 'matlabout');

            this = mlraut.AnalyticSignalGBMPar( ...
                subjects=sub, ...
                do_7T=false, ...
                do_plot_networks=false, ...
                do_resting=true, ...
                do_task=false, ...
                do_global_signal_regression=false, ...
                do_save=true, ...
                do_save_dynamic=false, ...
                do_save_ciftis=false, ...
                do_save_subset=false, ...
                force_band=false, ...
                hp_thresh=0.01, ...
                lp_thresh=0.1, ...
                out_dir=out_dir, ...
                source_physio=opts.source_physio, ...
                v_physio=50);
            call(this);
        end

        function server_create_heatmap(cores, opts)
            %% for servers

            arguments
                cores {mustBeScalarOrEmpty} = 1
                opts.N_sub {mustBeScalarOrEmpty} = []
                opts.flip_globbed logical = false
                opts.physio {mustBeTextScalar} = "CE"
            end

            root_dir = fullfile( ...
                getenv('SINGULARITY_HOME'), 'AnalyticSignalGBM', 'analytic_signal', 'dockerout', 'ciftify');
            matlabout_dir = fullfile( ...
                getenv('SINGULARITY_HOME'), 'AnalyticSignalGBM', 'analytic_signal', 'matlabout');
            subs_list = mglob(fullfile(matlabout_dir, sprintf('sub-I3CR*/*%s*.mat', opts.physio)));
            [~,subs_list] = fileparts(fileparts(subs_list));
            subs_list = unique(subs_list);
            map_filename = opts.physio + "_on_T1w.nii.gz";

            ifc = [];
            heatmap = zeros(91, 109, 91, 'double');
            for s1 = asrow(subs_list)
                try
                    fqfn = fullfile(root_dir, s1, "MNINonLinear", map_filename);
                    if ~isfile(fqfn)
                        continue
                    end
                    ifc = mlfourd.ImagingFormatContext2(fqfn);
                    heatmap = heatmap + double(ifc.img);
                catch ME
                    handwarning(ME)
                end
            end

            assert(~isempty(ifc))
            ifc.img = heatmap;
            ifc.filepath = matlabout_dir;
            ifc.fileprefix = mybasename(map_filename) + "_heatmap";
            save(ifc);
        end

        function server_find_misregistered(opts)
            %% for servers

            arguments
                opts.map_filename {mustBeTextScalar} = "CE_on_T1w.nii.gz"
            end

            root_dir = fullfile( ...
                getenv('SINGULARITY_HOME'), 'AnalyticSignalGBM', 'analytic_signal', 'dockerout', 'ciftify');
            subs_file = fullfile( ...
                getenv('SINGULARITY_HOME'), 'AnalyticSignalGBM', 'analytic_signal', 'matlabout', 'subs_that_complete.mat');
            ld = load(subs_file);
            subs_list = ld.subs;

            fprintf("%s for %s:\n", stackstr(), opts.map_filename)
            for s1 = asrow(subs_list)
                try
                    % find map
                    fqfn = fullfile(root_dir, "sub-"+s1, "MNINonLinear", opts.map_filename);
                    if ~isfile(fqfn)
                        continue
                    end
                    ifc = mlfourd.ImagingFormatContext2(fqfn);
                    
                    % build extra-cerebral mask
                    fqfn1 = fullfile(fileparts(fqfn), "T1w.nii.gz");
                    ic1 = mlfourd.ImagingContext2(fqfn1);
                    ic1 = ic1.numlt(0.001);

                    % tally extra-cerebral voxels
                    extra_vxls = sum(ifc.img .* ic1.imagingFormat.img > 0, "all");
                    if extra_vxls > 0
                        fprintf("%s: has %g extra voxels\n", s1, extra_vxls)
                    end
                catch ME
                    handwarning(ME)
                end
            end
            fprintf("\n")
        end

        function [j,c] = cluster_batch_call(globbing_mat, opts)
            %% for clusters running Matlab parallel server
            %  globbing_mat (text): on local machine calling parallel server
            %  sub_indices (double): nonempty for selecting subjects
            %  globbing_var (text) = "globbed": the object name of interest in globbing_mat

            arguments
                globbing_mat {mustBeFile} = ...
                    fullfile( ...
                    filesep, 'Users', 'jjlee', 'Singularity', 'AnalyticSignalGBM', 'analytic_signal', 'dockerout', 'ciftify', ...
                    'mlraut_AnalyticSignalGBMPar_globbing.mat')
                opts.sub_indices double = []
                opts.globbing_var = "globbed"
                opts.source_physio = "CE"
            end
            ld = load(globbing_mat);
            globbed = convertStringsToChars(ld.(opts.globbing_var));
            globbed = asrow(globbed);
            if ~isempty(opts.sub_indices)
                globbed = globbed(opts.sub_indices);
            end

            c = mlraut.CHPC3.propcluster();
            disp(c.AdditionalProperties)
            for g = globbed
                try
                    j = c.batch( ...
                        @mlraut.AnalyticSignalGBMPar.construct_and_call, ...
                        1, ...
                        {g(1), "source_physio", opts.source_physio}, ...
                        'CurrentFolder', '.', ...
                        'AutoAddClientPath', false);
                catch ME
                    handwarning(ME)
                end
            end

            % j = c.batch(@mlraut.AnalyticSignalHCPAgingPar.construct_and_call, 1, {}, 'CurrentFolder', '.', 'AutoAddClientPath', false);
        end

        function [j,c] = cluster_batch_call_2(globbing_mat, opts)
            %% for clusters running Matlab parallel server
            %  globbing_mat (text): on local machine calling parallel server
            %  globbing_var (text) = "globbed": the object name of interest in globbing_mat

            arguments
                globbing_mat {mustBeFile} = ...
                    fullfile( ...
                    filesep, 'Users', 'jjlee', 'Singularity', 'AnalyticSignalGBM', 'analytic_signal', 'dockerout', 'ciftify', ...
                    'mlraut_AnalyticSignalGBMPar_globbing.mat')
                opts.globbing_var = "globbed"
                opts.source_physio = "CE_on_T1w"
                opts.v_physio double = 50
            end
            ld = load(globbing_mat);
            globbed = convertStringsToChars(ld.(opts.globbing_var));
            globbed = asrow(globbed);

            % % find matlabout/sub-* which may need recovery
            % out_dir = fullfile( ...
            %     getenv('SINGULARITY_HOME'), 'AnalyticSignalGBM', 'analytic_signal', 'matlabout');
            % sub_dirs = mglob(fullfile(out_dir, 'sub-*'));
            % deselect = false(size(sub_dirs));
            % for idx = 1:length(sub_dirs)
            %     deselect(idx) = ~isempty(mglob(fullfile(sub_dirs(idx), '*iFV-brightest*.mat')));
            % end
            % sub_dirs(deselect) = [];
            % [~,subs] = fileparts(sub_dirs);
            % 
            % % select the subset of globbed that contains subs
            % select = false(size(globbed));
            % for idx = 1:length(subs)
            %     select = select | contains(globbed, subs(idx));
            % end
            % globbed = globbed(select);

            % pad and reshape globbed
            Ncol = 10;
            Nrow = ceil(numel(globbed)/Ncol);
            padding = repmat("", [1, Ncol*Nrow - numel(globbed)]);
            globbed = [globbed, padding];
            globbed = reshape(globbed, Nrow, Ncol);
            globbed = convertStringsToChars(globbed);
            fprintf("%s:globbed:\n", stackstr())
            disp(size(globbed))
            disp(ascol(globbed))

            % contact cluster slurm

            c = mlraut.CHPC3.propcluster();
            disp(c.AdditionalProperties)
            for irow = 1:Nrow
                try
                    j = c.batch( ...
                        @mlraut.AnalyticSignalGBMPar.construct_and_call, ...
                        1, ...
                        {globbed(irow, :), "source_physio", opts.source_physio, "v_physio", opts.v_physio}, ...
                        'Pool', Ncol, ...
                        'CurrentFolder', tempdir, ...
                        'AutoAddClientPath', false);
                catch ME
                    handwarning(ME)
                end
            end
        end

        function [j,c] = cluster_batch_call_fultz(globbing_mat, opts)
            %% for clusters running Matlab parallel server
            %  globbing_mat (text): on local machine calling parallel server
            %  sub_indices (double): nonempty for selecting subjects
            %  globbing_var (text) = "globbed": the object name of interest in globbing_mat

            arguments
                globbing_mat {mustBeFile} = ...
                    fullfile( ...
                    filesep, 'Users', 'jjlee', 'Singularity', 'AnalyticSignalGBM', 'analytic_signal', 'dockerout', 'ciftify', ...
                    'mlraut_AnalyticSignalGBMPar_globbing.mat')
                opts.sub_indices double = []
                opts.globbing_var = "globbed"
            end
            ld = load(globbing_mat);
            globbed = convertStringsToChars(ld.(opts.globbing_var));
            globbed = asrow(globbed);
            if ~isempty(opts.sub_indices)
                globbed = globbed(opts.sub_indices);
            end

            c = mlraut.CHPC3.propcluster_4gb_1h();
            disp(c.AdditionalProperties)
            for g = globbed
                for phys = {'iFV-brightest', 'CE_on_T1w', 'gray'}
                    try
                        j = c.batch( ...
                            @mlraut.AnalyticSignalGBMPar.construct_and_call_fultz, ...
                            1, ...
                            {g(1), "source_physio", phys{1}}, ...
                            'CurrentFolder', '.', ...
                            'AutoAddClientPath', false);
                    catch ME
                        handwarning(ME)
                    end
                end
            end

            % j = c.batch(@mlraut.AnalyticSignalHCPAgingPar.construct_and_call, 1, {}, 'CurrentFolder', '.', 'AutoAddClientPath', false);
        end

        function cluster_batch_link_kiyun()

            c = mlraut.CHPC3.propcluster_16gb_1h();
            disp(c.AdditionalProperties)

            % j = c.batch( ...
            %     @mlraut.CHPC3.parallel_example, ...
            %     1, ...
            %     {3}, ...
            %     'CurrentFolder', '.', ...
            %     'AutoAddClientPath', false);
 
            subs = {'sub-I3CR0266', 'sub-I3CR0639', 'sub-I3CR1023', 'sub-I3CR1837'};

            for s = subs
                j = c.batch( ...
                    @mlraut.AnalyticSignalGBMPar.link_kiyun, ...
                    1, ...
                    {s}, ...
                    'CurrentFolder', '.', ...
                    'AutoAddClientPath', false);
            end
        end

        function durations = construct_and_call(subjects, opts)
            %% must be called by batch()

            arguments
                subjects cell = {'sub-I3CR1488'}
                opts.tags {mustBeTextScalar} = "AnalyticSignalGBMPar"
                opts.out_dir {mustBeFolder} = fullfile( ...
                    filesep, 'scratch', 'jjlee', 'Singularity', 'AnalyticSignalGBM', 'analytic_signal', 'matlabout')
                opts.source_physio {mustBeTextScalar} = "CE"  % "iFV" "CE" "edema"
                opts.v_physio double = 50
            end
            durations = nan(1, length(subjects));
            diary(fullfile(opts.out_dir, subjects{1}, "diary.log"));

            parfor sidx = 1:length(subjects)

                tic;

                % populate tasks for each subject
                root_dir = fullfile( ...
                    filesep, 'scratch', 'jjlee', 'Singularity', 'AnalyticSignalGBM', 'analytic_signal', 'dockerout', 'ciftify');
                globbed_tasks = asrow(glob(fullfile(root_dir, subjects{sidx}, 'MNINonLinear', 'Results', '*rest*')));
                globbed_tasks = mybasename(globbed_tasks);

                % setup
                mlraut.CHPC3.setenvs();
                setenv("VERBOSITY", "1");
                ensuredir(opts.out_dir); %#ok<*PFBNS>
                ensuredir(fullfile(opts.out_dir, subjects{sidx}));

                % construct & call
                fprintf("constructing mlraut.AnalyticSignalGBMPar for %s\n", subjects{sidx});
                this = mlraut.AnalyticSignalGBMPar( ...
                    subjects=subjects{sidx}, ...
                    tasks=globbed_tasks, ...
                    do_resting=true, ...
                    do_task=false, ...
                    do_global_signal_regression=false, ...
                    do_save=true, ...
                    do_save_dynamic=false, ...
                    do_save_ciftis=false, ...
                    do_save_subset=false, ...
                    hp_thresh=0.01, ...
                    lp_thresh=0.1, ...
                    out_dir=opts.out_dir, ...
                    source_physio=opts.source_physio, ...
                    tags=opts.tags, ...
                    v_physio=opts.v_physio);                
                call(this);

                durations(sidx) = toc;
                fprintf("duration for %s: %s seconds", durations(sidx), subjects{sidx});
            end

            diary("off");
        end

        function duration = construct_and_call_fultz(subjects, opts)
            %% must be called by batch()

            arguments
                subjects cell = {'sub-I3CR1488'}
                opts.tasks cell = {}  % {'ses-1_task-rest_run-all_desc-preproc'}
                opts.tags {mustBeTextScalar} = "fultz"
                opts.out_dir {mustBeFolder} = fullfile( ...
                    filesep, 'scratch', 'jjlee', 'Singularity', 'AnalyticSignalGBM', 'analytic_signal', 'fultz')
                opts.source_physio {mustBeTextScalar} = "iFV-brightest"  % "iFV" "CE" "edema"
            end

            % populate opts.tasks as needed
            if isempty(opts.tasks)
                root_dir = fullfile( ...
                    filesep, 'scratch', 'jjlee', 'Singularity', 'AnalyticSignalGBM', 'analytic_signal', 'dockerout', 'ciftify');
                assert(~isempty(subjects))
                globbed = asrow(glob(fullfile(root_dir, subjects{1}, 'MNINonLinear', 'Results', '*rest*')));
                opts.tasks = mybasename(globbed);
            end
            
            mlraut.CHPC3.setenvs();
            setenv("VERBOSITY", "1");
            ensuredir(opts.out_dir);
            ensuredir(fullfile(opts.out_dir, subjects{1}));
            diary(fullfile(opts.out_dir, subjects{1}, "diary.log"));
            tic;
            disp("constructing mlraut.AnalyticSignalGBMPar")
            this = mlraut.AnalyticSignalGBMPar( ...
                subjects=subjects, ...
                tasks=opts.tasks, ...
                do_global_signal_regression=false, ...
                do_save=true, ...
                out_dir=opts.out_dir, ...   
                v_physio=50, ...
                plot_range=1:69, ...
                source_physio=opts.source_physio, ...
                hp_thresh=[], ...
                lp_thresh=0.1, ...
                rescaling="iqr", ...
                tag=opts.tags);
            disp("calling this")
            call(this);
            duration = toc;
            fprintf("tic-toc duration: %s seconds", duration);
            diary("off");
        end
    end

    methods
        function this = AnalyticSignalGBMPar(varargin)
            this = this@mlraut.AnalyticSignalGBM(varargin{:});
        end

        function this = call(this)
            %% CALL all subjects

            % exclude subjects

            out_dir_ = this.out_dir;
            for s = 1:this.num_sub
                this.current_subject = this.subjects{s};
                if ~contains(out_dir_, this.current_subject)
                    proposed_dir = fullfile(out_dir_, this.current_subject);
                    ensuredir(proposed_dir);
                    this.out_dir = proposed_dir;
                end 
                this.call_subject();
            end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
