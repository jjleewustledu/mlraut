classdef AnalyticSignalHCPPar < handle & mlraut.AnalyticSignalHCP
    %% line1
    %  line2
    %  
    %  Created 13-Apr-2023 02:11:52 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 9.14.0.2206163 (R2023a) for MACI64.  Copyright 2023 John J. Lee.
    
    methods (Static)
        function ihcp = mean_bold_ryans_way()
            out_dir = "/Volumes/PrecunealSSD2/AnalyticSignalHCP";
            mat = fullfile(out_dir, ...
                "sub-996782_ses-rfMRI-REST1-RL_proc-RV-gsr1-butter2-lp0p05-hp0p01-scaleiqr-Test-AnalyticSignalHCP-setupAnalyticSignalHCP.mat");            
            pwd0 = pushd(out_dir);

            ld = load(mat);
            ihcp = ld.this;

            % init
            num_bins = 40;
            bin_lims = asrow(linspace(-pi, pi, num_bins + 1));
            %midpoints = conv(bin_lims,[.5,.5],'valid'); 
            bins_ = zeros(num_bins, ihcp.num_nodes);

            % wrapped physio
            wrap_pupils = angle(ihcp.physio_signal);  % already entered, filtered, rescaled, analytic; Nt x 1

            % average bold by phase bins
            for b = 2:num_bins+1
                net_sigs2 = ihcp.bold_signal;  % already gsr, centered, filtered, rescaled, analytic; Nt x Ngo
                not_selected = wrap_pupils < bin_lims(b-1) | wrap_pupils > bin_lims(b);
                net_sigs2(not_selected, :) = nan;
                bins_(b-1,:) = mean(net_sigs2, 1, "omitnan");
            end

            % smooth bins
            smooth_width = 3;
            bins3 = bins_;
            bins2 = repmat(bins_, smooth_width, 1);
            for i = num_bins+1:2*num_bins
                bins3(i-num_bins, :) = nanmean(bins2(i-smooth_width:i+smooth_width, :), 1);
            end

            bins3 = real(bins3);
            ihcp.write_cifti(bins3, stackstr() + ".dtseries.nii");

            popd(pwd0);
        end

        function this = mean_twistor(nmats)
            arguments
                nmats {mustBeNumeric} = [1, inf]  % use finite for validity checks
            end

            warning("off", "MATLAB:class:LoadDefinitionUpdated");

            % \emph{this} supplies utilities
            out_dir = fullfile(getenv("SINGULARITY_HOME"), "AnalyticSignalHCP");
            this = mlraut.AnalyticSignalHCPPar( ...
                subjects={'996782'}, ...
                tasks={'rfMRI_REST1_RL'}, ...
                hp_thresh=0.01, ...
                lp_thresh=0.1, ...
                out_dir=out_dir, ...
                source_physio="iFV", ...
                tags="ASHCPPar-mean-twistor");
            call(this);  % malloc & construct delegates
            this.out_dir = out_dir;  % for group averages

            % init
            src_dir = fullfile(getenv("HOME"), "shared_shimony", "AnalyticSignalHCP");
            mats = asrow(glob(fullfile(src_dir, '*/sub-*_ses-*ASHCPPar*.mat')));
            if nmats(2) > length(mats)
                nmats(2) = length(mats);
            end
            mats = mats(nmats(1):nmats(2));  % 3690 mat files
            nbin = this.twistors.num_bins_angles;
            ngo = this.num_nodes;
            X_ = zeros(nbin, ngo);
            reY_ = zeros(nbin, ngo);
            imY_ = zeros(nbin, ngo);
            Z_ = zeros(nbin, ngo);
            T_ = zeros(nbin, ngo);
            r_ = zeros(1, ngo);
            bold_ = zeros(nbin, ngo);
            replvs_ = zeros(nbin, ngo);
            implvs_ = zeros(nbin, ngo);
            errs = 0;

            % abbrev.
            bin = @this.bin_by_physio_angle;

            for mat = mats
                tic
                try
                    ld = load(mat{1});
                    psi = ld.this_subset.bold_signal;
                    phi = ld.this_subset.physio_signal;
                    X = bin(this.X(psi, phi), phi);
                    Y = bin(this.Y(psi, phi), phi);
                    Z = bin(this.Z(psi, phi), phi);
                    T = bin(this.T(psi, phi), phi);
                    r = asrow(this.connectivity(psi, phi));
                    bold = bin(psi, phi);
                    plvs = bin(this.phase_locked_values(psi, phi), phi);

                    X_ = X_ + real(X);
                    reY_ = reY_ + real(Y);
                    imY_ = imY_ + imag(Y);
                    Z_ = Z_ + real(Z);
                    T_ = T_ + real(T);
                    r_ = r_ + real(r);
                    bold_ = bold_ + real(bold);
                    replvs_ = replvs_ + real(plvs);
                    implvs_ = implvs_ + imag(plvs);
                catch ME
                    fprintf("while working with %s: ", mat{1});
                    handwarning(ME)
                    errs = errs + 1;
                end
                fprintf("time working with %s: ", mat{1});
                toc
            end
            fprintf("errs->%g\n", errs)
            nmats_corr = nmats(end) - nmats(1) + 1 - errs;
            X_ = X_/nmats_corr;
            reY_ = reY_/nmats_corr;
            imY_ = imY_/nmats_corr;
            Z_ = Z_/nmats_corr;
            T_ = T_/nmats_corr;
            r_ = r_/nmats_corr;
            bold_ = bold_/nmats_corr;
            replvs_ = replvs_/nmats_corr;
            implvs_ = implvs_/nmats_corr;


            dt = this.twistors_.theta_berry/this.twistors_.num_bins_angles;
            units_t = "RADIAN";
            ntag = sprintf('n%g', nmats_corr);
            this.cifti.write_cifti( ...
                X_, sprintf('X_as_sub-%s_ses-%s_%s', ntag, ntag, this.tags), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                reY_, sprintf('reY_as_sub-%s_ses-%s_%s', ntag, ntag, this.tags), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                imY_, sprintf('imY_as_sub-%s_ses-%s_%s', ntag, ntag, this.tags), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                Z_, sprintf('Z_as_sub-%s_ses-%s_%s', ntag, ntag, this.tags), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                T_, sprintf('T_as_sub-%s_ses-%s_%s', ntag, ntag, this.tags), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                r_, sprintf('comparator_as_sub-%s_ses-%s_%s', ntag, ntag, this.tags), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                bold_, sprintf('bold_as_sub-%s_ses-%s_%s', ntag, ntag, this.tags), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                replvs_, sprintf('replvs_as_sub-%s_ses-%s_%s', ntag, ntag, this.tags), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                implvs_, sprintf('implvs_as_sub-%s_ses-%s_%s', ntag, ntag, this.tags), dt=dt, units_t=units_t);

            warning("on", "MATLAB:class:LoadDefinitionUpdated");
        end

        %% running call on single server

        function server_call(cores, opts)
            %% 33 GiB memory needed per instance of this running on a single process

            arguments
                cores {mustBeScalarOrEmpty} = 8
                opts.N_sub {mustBeScalarOrEmpty} = 100
                opts.flip_globbed logical = true
            end

            root_dir = '/home/usr/jjlee/mnt/CHPC_hcpdb/packages/unzip/HCP_1200';
            %root_dir = '/home/usr/jjlee/mnt/CHPC_scratch/Singularity/HcpAging/HCPAgingRec/fmriresults01';
            out_dir = '/home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalHCP';
            %out_dir = '/vgpool02/data2/jjlee/AnalyticSignalHcpAging';

            g = glob(fullfile(root_dir, '*'));
            if opts.flip_globbed
                g = flip(g); % examine more recent ones
            end
            g = cellfun(@(x) basename(x), g, UniformOutput=false);
            g = g(~contains(g, 'manifests'));
            g = g(101:100+opts.N_sub);
            leng = length(g);
            as = [];
            %for idxg = 1:1
            %parfor (idxg = 1:2, 2)
            parfor (idxg = 1:leng, cores)
                try
                    % if isfolder(fullfile(out_dir, g{idxg}))
                    %     continue
                    % end
                    as = mlraut.AnalyticSignalHCPPar( ...
                        subjects=g(idxg), ...  
                        tasks={'rfMRI_REST1_LR', 'rfMRI_REST1_RL', 'rfMRI_REST2_LR', 'rfMRI_REST2_RL'}, ...
                        do_save=true, ...
                        do_save_subset=true, ...
                        hp_thresh=0.01, ...
                        lp_thresh=0.1, ...
                        out_dir=out_dir, ...
                        source_physio=[ ...
                        "iFV", "sFV", "3rdV", "latV", "centrumsemiovale", ...
                        "ctx", "cbm", "str", "thal", ...
                        "RV-mean", "RV-std"], ...
                        tags="ASHCPPar-server-call");
                    call(as);
                catch ME
                    handwarning(ME)
                end
            end
        end
    
        %% running call, means, models on cluster

        function [j,c,msg,id] = cluster_construct_and_call(globbing_mat, opts)
            %% for clusters running Matlab parallel server

            arguments
                globbing_mat {mustBeFile} = ...
                    fullfile( ...
                    getenv('SINGULARITY_HOME'), 'AnalyticSignalHCP', 'mlraut_AnalyticSignalHCPPar_globbing.mat')
                opts.sub_indices double = 1:2  % total ~ 1:1113
                opts.globbing_var = "globbed"
                opts.Ncol {mustBeInteger} = 2
                opts.account_name char = 'joshua_shimony'
                opts.skip_existing logical = true
            end
            ld = load(globbing_mat);
            globbed = convertStringsToChars(ld.(opts.globbing_var));
            globbed = asrow(globbed);
            if ~isempty(opts.sub_indices)
                globbed = globbed(opts.sub_indices);
            end

            % pad and reshape globbed
            Nrow = ceil(numel(globbed)/opts.Ncol);
            padding = repmat("", [1, opts.Ncol*Nrow - numel(globbed)]);
            globbed = [globbed, padding];
            globbed = reshape(globbed, Nrow, opts.Ncol);
            globbed = convertStringsToChars(globbed);
            fprintf("%s:globbed:\n", stackstr())
            disp(size(globbed))
            disp(ascol(globbed))

            % contact cluster slurm

            warning('off', 'MATLAB:legacy:batchSyntax');
            warning('off', 'parallel:convenience:BatchFunctionNestedCellArray');
            warning('off', 'MATLAB:TooManyInputs');

            c = mlraut.CHPC3.propcluster(opts.account_name, mempercpu='64gb', walltime='3:00:00');
            disp(c.AdditionalProperties)
            for irow = 1:Nrow
                try
                    j = c.batch( ...
                        @mlraut.AnalyticSignalHCPPar.construct_and_call, ...
                        1, ...
                        {globbed(irow, :), 'skip_existing', opts.skip_existing}, ...
                        'Pool', opts.Ncol, ...
                        'CurrentFolder', '/scratch/jjlee/Singularity/AnalyticSignalHCP', ...
                        'AutoAddClientPath', false);
                catch ME
                    handwarning(ME)
                end
            end

            [msg,id] = lastwarn();

            % j = c.batch(@mlraut.AnalyticSignalHCPPar.construct_and_call, 1, {}, 'CurrentFolder', '.', 'AutoAddClientPath', false);
        end

        function [j,c,msg,id] = cluster_construct_means(globbing_mat, opts)
            %% for clusters running Matlab parallel server

            arguments
                globbing_mat {mustBeText} = ...
                    fullfile( ...
                    getenv('SINGULARITY_HOME'), 'AnalyticSignalHCP', 'mlraut_AnalyticSignalHCPPar_globbing.mat')
                opts.globbing_var = "globbed"
                opts.sub_range = []  % total ~ 1:1113
                opts.new_physio {mustBeText} = "RV-std"
                opts.anatomy {mustBeTextScalar} = "ctx"
                opts.tag {mustBeText} = "locuscer"
                opts.test_range = []  % 1:2
                opts.Ncol {mustBeInteger} = 32
                opts.account_name char = 'joshua_shimony' % 'aristeidis_sotiras'
                opts.instance {mustBeTextScalar} = "mean_Z_sup_alpha_instance"
            end

            % gather globbed
            if isscalar(globbing_mat) && endsWith(globbing_mat, ".mat") && ~isemptytext(opts.globbing_var)
                ld = load(globbing_mat);
                globbed = ld.(opts.globbing_var);
            else
                globbed = globbing_mat;
            end
            globbed = asrow(globbed);
            if ~isempty(opts.sub_range)
                globbed = globbed(opts.sub_range);
            end

            % pad and reshape globbed
            Nrow = ceil(numel(globbed)/opts.Ncol);
            padding = repmat("", [1, opts.Ncol*Nrow - numel(globbed)]);
            globbed = [globbed, padding];
            globbed = reshape(globbed, Nrow, opts.Ncol);
            globbed = convertStringsToChars(globbed);
            fprintf("%s:globbed:\n", stackstr())
            disp(size(globbed))
            disp(globbed)

            % contact cluster slurm

            warning('off', 'MATLAB:legacy:batchSyntax');
            warning('off', 'parallel:convenience:BatchFunctionNestedCellArray');
            warning('off', 'MATLAB:TooManyInputs');

            c = mlraut.CHPC3.propcluster(opts.account_name, mempercpu='48gb', walltime='24:00:00');
            disp(c.AdditionalProperties)

            for col_idx = 1:opts.Ncol
                try
                    j = c.batch( ...
                        @mlraut.AnalyticSignalHCPPar.construct_means, ...
                        1, ...
                        {globbed(:, col_idx), ...
                            'col_idx', col_idx, 'new_physio', opts.new_physio, 'test_range', opts.test_range, ...
                            'tag', opts.tag, 'anatomy', opts.anatomy, 'instance', opts.instance}, ...
                        'CurrentFolder', '/scratch/jjlee/Singularity/AnalyticSignalHCP', ...
                        'AutoAddClientPath', false);
                catch ME
                    handwarning(ME)
                end
            end

            warning('on', 'MATLAB:legacy:batchSyntax');
            warning('on', 'parallel:convenience:BatchFunctionNestedCellArray');
            warning('on', 'MATLAB:TooManyInputs');

            [msg,id] = lastwarn();
        end

        function [j,c,msg,id] = cluster_construct_models(globbing_mat, opts)
            %% for clusters running Matlab parallel server

            arguments
                globbing_mat {mustBeFile} = ...
                    fullfile( ...
                    getenv('SINGULARITY_HOME'), 'AnalyticSignalHCP', 'mlraut_AnalyticSignalHCPPar_globbing.mat')
                opts.globbing_var = "globbed"
                opts.sub_range = []  % total ~ 1:1113
                opts.new_physio {mustBeText} = "iFV"
                opts.tag {mustBeText} = ""
                opts.test_range = []  % 1:2
                opts.Ncol {mustBeInteger} = 32
                opts.account_name char = 'joshua_shimony'
            end
            ld = load(globbing_mat);
            globbed = convertStringsToChars(ld.(opts.globbing_var));
            globbed = asrow(globbed);
            if ~isempty(opts.sub_range)
                globbed = globbed(opts.sub_range);
            end

            % pad and reshape globbed
            Nrow = ceil(numel(globbed)/opts.Ncol);
            padding = repmat("", [1, opts.Ncol*Nrow - numel(globbed)]);
            globbed = [globbed, padding];
            globbed = reshape(globbed, Nrow, opts.Ncol);
            globbed = convertStringsToChars(globbed);
            fprintf("%s:globbed:\n", stackstr())
            disp(size(globbed))
            disp(ascol(globbed))

            % contact cluster slurm

            warning('off', 'MATLAB:legacy:batchSyntax');
            warning('off', 'parallel:convenience:BatchFunctionNestedCellArray');
            warning('off', 'MATLAB:TooManyInputs');

            c = mlraut.CHPC3.propcluster(opts.account_name, mempercpu='40gb', walltime='24:00:00');
            disp(c.AdditionalProperties)

            for col_idx = 1:opts.Ncol
                try
                    j = c.batch( ...
                        @mlraut.AnalyticSignalHCPPar.construct_models, ...
                        1, ...
                        {globbed(:, col_idx), ...
                            'col_idx', col_idx, 'new_physio', opts.new_physio, 'test_range', opts.test_range, ...
                            'tag', opts.tag}, ...
                        'CurrentFolder', '/scratch/jjlee/Singularity/AnalyticSignalHCP', ...
                        'AutoAddClientPath', false);
                catch ME
                    handwarning(ME)
                end
            end

            warning('on', 'MATLAB:legacy:batchSyntax');
            warning('on', 'parallel:convenience:BatchFunctionNestedCellArray');
            warning('on', 'MATLAB:TooManyInputs');

            [msg,id] = lastwarn();
        end

        function durations = construct_and_call(subjects, opts)
            arguments
                subjects cell = {'996782'}
                opts.tasks cell = {'rfMRI_REST1_LR', 'rfMRI_REST1_RL', 'rfMRI_REST2_LR', 'rfMRI_REST2_RL'}
                opts.tags {mustBeTextScalar} = "ASHCPPar"
                opts.out_dir {mustBeFolder} = fullfile(getenv("HOME"), "Singularity", "AnalyticSignalHCP")
                opts.skip_existing logical = true
            end
            
            subjects = subjects(~cellfun(@isempty, subjects));  % Remove empty cells
            durations = nan(1, length(subjects));

            %parfor sidx = 1:length(subjects)
            for sidx = 1:length(subjects)

                tic;
            
                % setup
                mlraut.CHPC3.setenvs();
                ensuredir(opts.out_dir); %#ok<*PFBNS>
                ensuredir(fullfile(opts.out_dir, subjects{sidx}));

                if opts.skip_existing
                    g = mglob(fullfile(opts.out_dir, subjects{sidx}, "*.mat"));
                    if length(g) > 1
                        continue
                    end
                end

                try
                    % construct & call
                    as = mlraut.AnalyticSignalHCPPar( ...
                        subjects=subjects{sidx}, ...
                        tasks=opts.tasks, ...
                        do_save=true, ...
                        do_save_subset=true, ...
                        hp_thresh=0.01, ...
                        lp_thresh=0.1, ...
                        out_dir=opts.out_dir, ...
                        source_physio=[ ...
                        "iFV", "sFV", "3rdV", "latV", "centrumsemiovale", ...
                        "ctx", "cbm", "str", "thal", ...
                        "RV-mean", "RV-std"], ...
                        tags=opts.tags);
                    call(as);
                catch ME
                    handwarning(ME)
                end

                durations(sidx) = toc;
            end
        end
    
        function durations = construct_means(subjects, opts)
            arguments
                subjects cell = {'996782'}
                opts.tasks cell = {'rfMRI_REST1_LR', 'rfMRI_REST1_RL', 'rfMRI_REST2_LR', 'rfMRI_REST2_RL'}
                opts.tags {mustBeTextScalar} = "ASHCPPar-construct-means"
                opts.out_dir {mustBeFolder} = "/scratch/jjlee/Singularity/AnalyticSignalHCP"
                opts.col_idx {mustBeScalarOrEmpty} = nan
                opts.new_physio {mustBeText} = "iFV"
                opts.anatomy {mustBeTextScalar} = "ctx"
                opts.test_range = []
                opts.tag {mustBeText} = ""
                opts.instance {mustBeTextScalar} = "mean_twistor_instance"
            end

            tic;

            % setup
            mlraut.CHPC3.setenvs();
            ensuredir(opts.out_dir); %#ok<*PFBNS>

            mat_fqfn = [];
            idx = 0;
            while isempty(mat_fqfn) && idx < length(subjects)
                idx = idx + 1;
                if endsWith(opts.out_dir, subjects{idx})
                    out_dir_ = fileparts(opts.out_dir);
                else
                    out_dir_ = opts.out_dir;
                end
                patt = fullfile( ...
                    out_dir_, ...
                    subjects{idx}, ...
                    sprintf("sub-%s_ses-*_proc-*ASHCPPar*.mat", subjects{idx}));
                % fprintf("%s: %s\n", stackstr(), patt);
                mat_fqfn = mglob(patt);
                % fprintf("mat_fqfn: \n")
                % disp(mat_fqfn)
            end
            as = mlraut.AnalyticSignalHCPPar.load(mat_fqfn(1), class="mlraut.AnalyticSignalHCPPar");
            as.out_dir = opts.out_dir;
            as.(opts.instance)( ...
                subjects, col_idx=opts.col_idx, new_physio=opts.new_physio, test_range=opts.test_range, ...
                tag=opts.tag, anatomy=opts.anatomy);

            durations = toc;
        end
    
        function durations = construct_models(subjects, opts)
            arguments
                subjects cell = {'996782'}
                opts.tasks cell = {'rfMRI_REST1_LR', 'rfMRI_REST1_RL', 'rfMRI_REST2_LR', 'rfMRI_REST2_RL'}
                opts.tags {mustBeTextScalar} = "ASHCPPar-construct-models"
                opts.out_dir {mustBeFolder} = "/scratch/jjlee/Singularity/AnalyticSignalHCP"
                opts.col_idx {mustBeInteger}
                opts.new_physio {mustBeText} = "iFV"
                opts.test_range = []
                opts.tag {mustBeText} = ""
            end
            subjects = convertCharsToStrings(subjects);

            tic;

            % setup
            mlraut.CHPC3.setenvs();
            ensuredir(opts.out_dir); %#ok<*PFBNS>

            for sub = asrow(subjects)
                mat_fqfns = mglob(fullfile( ...
                    opts.out_dir, sub, sprintf("sub-%s_ses-*_proc-*ASHCPPar*.mat", sub)));
                for mat = asrow(mat_fqfns)
                    ld = load(mat);
                    this_subset = ld.this_subset;
                    psi = this_subset.bold_signal;
                    phi = hilbert(this_subset.physio_supplementary(opts.new_physio));

                    %% generative model
                    [mdl.xi,mdl.alpha_est,mdl.residual_err,mdl.mse] = mlraut.AnalyticSignalHCPPar.generative_model( ...
                        psi, phi, Niter=20);
                    mdl.residual = psi - mdl.xi;

                    %% greyordinate MSE
                    mdl.mesh.mse_go = mean(abs(mdl.residual).^2, 1);

                    %% greyordinate R^2

                    % Calculate means for each greyordinate
                    psi_mean = mean(psi, 1); % [1 x N_greyords]

                    % Total sum of squares for each greyordinate
                    TSS = sum((psi - psi_mean).^2, 1); % [1 x N_greyords]

                    % Residual sum of squares for each greyordinate
                    RSS = sum((mdl.residual).^2, 1); % [1 x N_greyords]

                    % R-squared for each greyordinate
                    R_squared = 1 - (RSS ./ TSS); % [1 x N_greyords]

                    % Handle edge case where TSS = 0 (constant signal)
                    R_squared(TSS == 0) = NaN;
                    mdl.mesh.R2_go = R_squared;

                    %% save
                    mdl_fqfn = strrep(mat, "-ASHCPPar", "-mdl-ASHCPPar");
                    save(mdl_fqfn, "mdl", "-v7.3");
                end
            end

            durations = toc;
        end
    
        function gather_means_dtseries(out_dir, opts)
            %% generates only dtseries, for viewing phase variations

            arguments
                out_dir {mustBeFolder} = pwd
                opts.measures {mustBeText} = ["reZsup0", "reZsup1", "reZsup2", "imZsup0", "imZsup1", "imZsup2", "rezeta", "imzeta", "T", "X", "Y", "Z", "plvs", "bold"]
                opts.physio {mustBeText} = "RV-std"
                opts.gsr logical = false
                opts.ddt logical = true
                opts.tag {mustBeTextScalar} = "ASHCPPar*_par*"
                opts.new_tag {mustBeTextScalar} = "meanfield"
                opts.Nphase_out {mustBeScalarOrEmpty} = 80
                opts.theta_berry {mustBeScalarOrEmpty} = 4*pi  % mlraut.Twistors.theta_berry
                opts.Fs {mustBeScalarOrEmpty} = 1/0.72
            end
            if ~isemptytext(opts.new_tag) && ~startsWith(opts.new_tag, "-") && ~startsWith(opts.new_tag, "_")
                opts.new_tag = "-" + opts.new_tag;
            end

            for measure = opts.measures
                try
                    globbed = asrow(mglob( ...
                        fullfile(out_dir, ...
                        sprintf("%s_as_sub-*_ses-*_proc-%s-gsr%g-ddt%g*%s*.nii", ...
                        measure, opts.physio, opts.gsr, opts.ddt, opts.tag))));

                    % for all repeated measures
                    g_last = "/path/to/a_measurement.nii";
                    cii = [];
                    mats = {};
                    re = [];
                    n = 0;
                    weights = [];
                    for g = globbed
                        try
                            g_last = g;
                            cii = cifti_read(g);
                            if any(~isfinite(cii.cdata), "all")
                                continue
                            end
                            mats = [mats, cii.cdata'];  %#ok<*AGROW> % {Nt x Ngo} x Nmats; Nt x Ngo is preferred for internal repr.
                            re = regexp(g, "\S_sub-n(?<subn>\d+)_ses-n(?<sesn>\d+)_\S", "names");
                            subn = str2double(re.subn);
                            sesn = str2double(re.sesn);
                            assert(subn == sesn)
                            n = n + subn;
                            weights = [weights, subn];
                        catch
                            continue
                        end
                    end
                    if isempty(mats)
                        continue
                    end
                    weights = weights / n;

                    % apply weights to mats
                    mats_weighted = cellfun(@(mat, weight) mat * weight, mats, num2cell(weights), 'UniformOutput', false);

                    % concatenate arrays along the 3rd dimension, then sum
                    mat_avg = sum(cat(3, mats_weighted{:}), 3);  % ~ N_theta x N_go
                    Nnu = size(mat_avg, 1);
                    assert(Nnu > 1)
                    assert(mod(Nnu, opts.Nphase_out) < eps)

                    % filename
                    re = regexp(g_last, "\S_sub-n(?<subn>\d+)_ses-n(?<sesn>\d+)_\S", "names");  % last cii may be defective & re is from `last - 1`
                    fqfn = strrep(g_last, re.subn, string(n));
                    fqfn = regexprep(fqfn, "_par\d+", "");  % remove par tag
                    fqfn = regexprep(fqfn, '-\d{14}(?=\W|$)', '');  % remove datetime tag
                    if ~isemptytext(opts.new_tag)
                        fqfn = string(fqfn);
                        [pth,fp,x] = myfileparts(fqfn);
                        fp = strrep(fp, "-ASHCPPar", opts.new_tag + "-ASHCPPar");
                        fqfn = fullfile(pth, fp + x);
                    end

                    if ~isempty(opts.Nphase_out) && opts.Nphase_out < Nnu
                        % down-sample phase-ordered measure ~ 80 x N_go -> 8 x N_go

                        % manage theta_berry
                        [mat_avg,Nnu] = reconcile_theta_berry(mat_avg, opts.theta_berry, mlraut.Twistors.theta_berry);

                        Navg = Nnu / opts.Nphase_out;  % num. of frames to average to satisfy opts.Nphase_out
                        mat_avg = squeeze(mean(reshape(mat_avg, Navg, opts.Nphase_out, []), 1));  % requires Navg \in \mathbb{Z}^+
                    end

                    % save weighted average
                    cii.diminfo{2} = cifti_diminfo_make_series( ...
                        opts.Nphase_out, 0, opts.theta_berry/opts.Nphase_out, 'RADIAN');
                    cii.cdata = mat_avg';
                    cifti_write(cii, convertStringsToChars(fqfn));
                catch ME
                    handwarning(ME)
                end
            end

            function [mat_avg_1,Nnu] = reconcile_theta_berry(mat_avg, theta_berry_1, theta_berry_0)
                %% e.g., reconciles phase range [-2*pi, 2*pi] => [-pi, pi]
                %  for theta_berry_1 = 2*pi, theta_berry_0 = 4*pi.

                assert(theta_berry_1 <= theta_berry_0)
                assert(mod(theta_berry_0, theta_berry_1) == 0)
                Nberry = round(theta_berry_0 / theta_berry_1);  % 2
                Ngo = size(mat_avg, 2);  % 91282
                Ntheta = size(mat_avg, 1);  % 80
                Ntheta_1 = round(Ntheta / Nberry);  % 40

                mat_avg_1 = zeros(Ntheta_1, Ngo);
                for b = 1:Nberry
                    theta_left = (b - 1) * Ntheta_1 + 1;
                    theta_right = b * Ntheta_1;
                    mat_avg_1 = mat_avg_1 + mat_avg(theta_left:theta_right, :);
                end
                mat_avg_1 = mat_avg_1 / Nberry;                
                Nnu = size(mat_avg_1, 1);
            end
        end

        function gather_means_dscalar(out_dir, opts)
            %% generates only dscalar, averaging out all phases

            arguments
                out_dir {mustBeFolder} = pwd
                opts.measures {mustBeText} = ["T", "X", "Y", "Z", "plvs", "comparator"]
                opts.physio {mustBeText} = "RV-std"
                opts.gsr logical = false
                opts.ddt logical = true
                opts.tag {mustBeTextScalar} = "ASHCPPar*_par*"
                opts.new_tag {mustBeTextScalar} = "meanfield"
            end

            for measure = opts.measures
                try
                    globbed = asrow(mglob( ...
                        fullfile(out_dir, ...
                        sprintf("%s_as_sub-*_ses-*_proc-%s-gsr%g-ddt%g*%s*.nii", ...
                        measure, opts.physio, opts.gsr, opts.ddt, opts.tag))));

                    % for all repeated measures
                    g_last = "/path/to/a_measurement.nii";
                    cii = [];
                    mats = {};
                    re = [];
                    n = 0;
                    weights = [];
                    for g = globbed
                        try
                            g_last = g;
                            cii = cifti_read(g);
                            if any(~isfinite(cii.cdata))
                                continue
                            end
                            mats = [mats, cii.cdata'];  %#ok<*AGROW> % {Nt x Ngo} x Nmats; Nt x Ngo is preferred for internal repr.
                            re = regexp(g, "\S_sub-n(?<subn>\d+)_ses-n(?<sesn>\d+)_\S", "names");
                            subn = str2double(re.subn);
                            sesn = str2double(re.sesn);
                            assert(subn == sesn)
                            n = n + subn;
                            weights = [weights, subn];
                        catch
                            continue
                        end
                    end
                    if isempty(mats)
                        continue
                    end
                    weights = weights / n;

                    % apply weights to mats
                    mats_weighted = cellfun(@(mat, weight) mat * weight, mats, num2cell(weights), 'UniformOutput', false);

                    % concatenate arrays along the 3rd dimension, then sum
                    mat_avg = sum(cat(3, mats_weighted{:}), 3);

                    % filename
                    fqfn = strrep(g_last, re.subn, string(n));
                    fqfn = regexprep(fqfn, "_par\d+", "");  % remove par tag
                    fqfn = regexprep(fqfn, '-\d{14}(?=\W|$)', '');  % remove datetime tag
                    if ~isemptytext(opts.new_tag)
                        fqfn = string(fqfn);
                        [pth,fp,x] = myfileparts(fqfn);
                        fp = strrep(fp, "-ASCPPar", opts.new_tag + "-ASCPPar");
                        fqfn = fullfile(pth, fp + x);
                    end

                    % take mean of all samples
                    mat_avg = mean(mat_avg, 1);
                    cii.diminfo{2} = cifti_diminfo_make_scalars(1);
                    fqfn = strrep(fqfn, ".dtseries.nii", ".dscalar.nii");

                    % save weighted average
                    cii.cdata = mat_avg';
                    cifti_write(cii, convertStringsToChars(fqfn));
                catch ME
                    handwarning(ME)
                end
            end
        end
    
        %% mean-field methods

        function alpha_est = em_alpha_estimation(M1, M2, max_iter, unit_modulus)
            if nargin < 3, max_iter = 20; end
            if nargin < 4, unit_modulus = true; end

            m1_vec = M1(:);
            m2_vec = M2(:);

            % Initialize
            alpha_init = (m1_vec' * m2_vec) / (m1_vec' * m1_vec);
            if unit_modulus
                alpha_est = alpha_init / abs(alpha_init);
            else
                alpha_est = alpha_init;
            end

            sigma2 = var(abs(m2_vec - alpha_est * m1_vec));

            for iter = 1:max_iter
                % E-step: compute responsibilities
                residuals = m2_vec - alpha_est * m1_vec;
                likelihood = exp(-abs(residuals).^2 / (2*sigma2));
                weights = likelihood / sum(likelihood);

                % M-step: update parameters
                alpha_old = alpha_est;

                if unit_modulus
                    % Enhanced phase estimation using circular statistics
                    % Compute element-wise phase ratios where M1 is significant
                    magnitude_threshold = 0.1 * max(abs(m1_vec));
                    valid_mask = abs(m1_vec) > magnitude_threshold;

                    if sum(valid_mask & weights > 1e-6) > 0
                        % Use weighted circular mean for robust phase estimation
                        phase_ratios = m2_vec(valid_mask) ./ m1_vec(valid_mask);
                        valid_weights = weights(valid_mask);

                        % Weighted circular mean
                        weighted_sum = sum(valid_weights .* phase_ratios);
                        alpha_est = weighted_sum / abs(weighted_sum);
                    else
                        % Fallback to correlation-based method
                        weighted_correlation = sum(weights .* conj(m1_vec) .* m2_vec);
                        alpha_est = weighted_correlation / abs(weighted_correlation);
                    end
                else
                    % Unconstrained update
                    alpha_est = sum(weights .* conj(m1_vec) .* m2_vec) / ...
                        sum(weights .* abs(m1_vec).^2);
                end

                % Update noise variance
                sigma2 = sum(weights .* abs(m2_vec - alpha_est * m1_vec).^2) / sum(weights);

                % Convergence check
                if abs(alpha_est - alpha_old) < 1e-6
                    break;
                end
            end
        end

        function [xi_mu,alpha_est,residual_err,mse] = generative_model(psi_mu, phi_mu, opts)
            %% Generates complex BOLD timeseries from X, Y, Z, T that have refinements for the phase of 
            %  the complex physio timeseries.
            %
            %  Args:
            %      psi_mu {mustBeMatrix} : ~ complex Nt x Ngo
            %      phi_mu {mustBeVector} : ~ complex Nt x 1
            %      opts.model struct = struct( ...
            %          "physio", "iFV", "gsr", true, "ddt", true, "butter", 8, "lp_thresh", 0.1, "hp_thresh", 0.01)
            %      opts.Niter {mustBeScalarOrEmpty} = 20 : iterative if Niter < 20 else EM; no alpha_est if Niter empty
            %      opts.t_range = 200:896  % time frames
            %      opts.g_range = 1:2*32492  % cortical vertices
            %
            %  Returns:
            %      xi_mu ~ Nt x Ngo
            %      alpha_est ~ 1
            %      residual_err ~ 1
            
            arguments
                psi_mu {mustBeMatrix}
                phi_mu {mustBeVector}
                opts.model struct = struct( ...
                    "physio", "iFV", "gsr", true, "ddt", true, "butter", 8, "lp_thresh", 0.1, "hp_thresh", 0.01)
                opts.Niter {mustBeScalarOrEmpty} = []
                opts.t_range = 200:896  % time frames
                opts.g_range = 1:2*32492  % cortical vertices
            end
            mdl = opts.model;  % abbrev

            % models generate from -dpsi/dt which have one less time sample
            Nt = size(phi_mu, 1) - 1;
            phi_mu = phi_mu(1:Nt, :);

            % mean field in Euclidean coords <- averaged spectra
            overbar_X_plus_iY = mlraut.AnalyticSignalHCPPar.overbar_A_nu(mdl, measure="X") + ...
                1i * mlraut.AnalyticSignalHCPPar.overbar_A_nu(mdl, measure="Y");
            overbar_X_plus_iY = ifft(overbar_X_plus_iY);

            overbar_T_minus_Z = mlraut.AnalyticSignalHCPPar.overbar_A_nu(mdl, measure="T") - ...
                mlraut.AnalyticSignalHCPPar.overbar_A_nu(mdl, measure="Z");
            overbar_T_minus_Z = ifft(overbar_T_minus_Z);

            neg_dxi_mu_dt = phi_mu .* overbar_X_plus_iY ./ overbar_T_minus_Z;
            xi_mu = mlraut.AnalyticSignalHCPPar.reconstruct_xi(neg_dxi_mu_dt, psi_mu);
            if isempty(opts.Niter)
                alpha_est = 1;
                xi_ = xi_mu(opts.t_range, opts.g_range);
                psi_ = psi_mu(opts.t_range, opts.g_range);
                residual = psi_ - xi_;
                residual_err = norm(residual, 'fro') / norm(psi_, 'fro');
                mse = mean(abs(residual(:)).^2);
                return
            end
            [xi_mu,alpha_est,residual_err,mse] = mlraut.AnalyticSignalHCPPar.optimize_phase_factor( ...
                xi_mu, psi_mu, Niter=opts.Niter, t_range=opts.t_range, g_range=opts.g_range);
        end

        function [neg_dxi_mu_dt,alpha_est,residual_err,mse] = generative_model_momentum(neg_dpsi_mu_dt, phi_mu, opts)
            %% Generates complex BOLD timeseries from X, Y, Z, T that have refinements for the phase of 
            %  the complex physio timeseries.
            %
            %  Args:
            %      neg_dpsi_mu_dt {mustBeMatrix} : ~ complex Nt x Ngo
            %      phi_mu {mustBeVector} : ~ complex Nt x 1
            %      opts.model struct = struct( ...
            %          "physio", "iFV", "gsr", true, "ddt", true, "butter", 8, "lp_thresh", 0.1, "hp_thresh", 0.01)
            %      opts.Niter {mustBeScalarOrEmpty} = 20 : iterative if Niter < 20 else EM; no alpha_est if Niter empty
            %      opts.t_range = 200:896  % time frames
            %      opts.g_range = 1:2*32492  % cortical vertices
            %
            %  Returns:
            %      neg_dxi_mu_dt ~ Nt x Ngo
            %      alpha_est ~ 1
            %      residual_err ~ 1
            
            arguments
                neg_dpsi_mu_dt {mustBeMatrix}
                phi_mu {mustBeVector}
                opts.model struct = struct( ...
                    "physio", "iFV", "gsr", true, "ddt", true, "butter", 8, "lp_thresh", 0.1, "hp_thresh", 0.01)
                opts.Niter {mustBeScalarOrEmpty} = []
                opts.t_range = 200:896  % time frames
                opts.g_range = 1:2*32492  % cortical vertices
            end
            mdl = opts.model;  % abbrev

            % models generate from -dpsi/dt which have one less time sample
            Nt = size(phi_mu, 1) - 1;
            phi_mu = phi_mu(1:Nt, :);

            % mean field in Euclidean coords <- averaged spectra
            overbar_X_plus_iY = mlraut.AnalyticSignalHCPPar.overbar_A_nu(mdl, measure="X") + ...
                1i * mlraut.AnalyticSignalHCPPar.overbar_A_nu(mdl, measure="Y");
            overbar_X_plus_iY = ifft(overbar_X_plus_iY);

            overbar_T_minus_Z = mlraut.AnalyticSignalHCPPar.overbar_A_nu(mdl, measure="T") - ...
                mlraut.AnalyticSignalHCPPar.overbar_A_nu(mdl, measure="Z");
            overbar_T_minus_Z = ifft(overbar_T_minus_Z);

            neg_dxi_mu_dt = phi_mu .* overbar_X_plus_iY ./ overbar_T_minus_Z;
            scaling = iqr(abs(neg_dpsi_mu_dt(opts.t_range, opts.g_range)), "all") / ...
                iqr(abs(neg_dxi_mu_dt(opts.t_range, opts.g_range)), "all");
            neg_dxi_mu_dt = scaling * neg_dxi_mu_dt;
            if isempty(opts.Niter)
                alpha_est = 1;
                modelled_ = neg_dxi_mu_dt(opts.t_range, opts.g_range);
                measured_ = neg_dpsi_mu_dt(opts.t_range, opts.g_range);
                residual = measured_ - modelled_;
                residual_err = norm(residual, 'fro') / norm(measured_, 'fro');
                mse = mean(abs(residual(:)).^2);
                return
            end
            [neg_dxi_mu_dt,alpha_est,residual_err,mse] = mlraut.AnalyticSignalHCPPar.optimize_phase_factor( ...
                neg_dxi_mu_dt, neg_dpsi_mu_dt, Niter=opts.Niter, t_range=opts.t_range, g_range=opts.g_range);
        end

        function phase_factor = iterative_alpha_estimation(M1, M2, max_iter)

            % Start with inner-product estimate
            alpha_est = (M2(:)' * M1(:)) / (M1(:)' * M1(:));
            phase_factor = alpha_est / abs(alpha_est);

            % Refine by removing outliers
            for iter = 1:max_iter
                residual = M2 - phase_factor * M1;
                weights = 1 ./ (1 + abs(residual).^2 / median(abs(residual(:))).^2);

                % Weighted least squares update
                w_vec = weights(:);
                m1_vec = M1(:);
                m2_vec = M2(:);

                alpha_est = sum(w_vec .* conj(m1_vec) .* m2_vec) / sum(w_vec .* abs(m1_vec).^2);
                phase_factor = alpha_est / abs(alpha_est);
            end
        end

        function [xi1,alpha_est,residual_err,mse] = optimize_phase_factor(xi, psi, opts)
            %% Adjust complex matrix xi to match complex matrix psi, estimating a scalar phase factor between them, 
            %  such that psi ≈ α·xi + noise, and α is complex.
            %  EM if Niter >=20 else iterative.

            arguments
                xi {mustBeMatrix}
                psi {mustBeMatrix}
                opts.Niter {mustBeScalarOrEmpty} = 20
                opts.t_range = 200:896  % time frames
                opts.g_range = 1:2*32492  % cortical vertices
            end

            % Limit optimizations to cortical vertices
            xi_ = xi(opts.t_range, opts.g_range);
            psi_ = psi(opts.t_range, opts.g_range);

            % expectation maximization
            if opts.Niter < 20
                alpha_est = mlraut.AnalyticSignalHCPPar.iterative_alpha_estimation(xi_, psi_, opts.Niter);
            else
                alpha_est = mlraut.AnalyticSignalHCPPar.em_alpha_estimation(xi_, psi_, opts.Niter);
            end
            xi1 = xi * alpha_est;

            % Assess estimation quality
            residual = psi_ - xi_ * alpha_est;
            residual_err = norm(residual, 'fro') / norm(psi_, 'fro');
            mse = mean(abs(residual(:)).^2);

            % Diagnostics
            % fprintf("%s mse: %g\n", stackstr(), mse);
            % figure; histogram(abs(residual(:))); title('Residual Magnitudes');
        end

        function A_nu = overbar_A_nu(mdl, opts)
            arguments
                mdl struct
                opts.measure {mustBeTextScalar}
            end

            pth = fullfile(getenv("SINGULARITY_HOME"), "AnalyticSignalHCP", mdl.physio + "-fft");
            g1 = mglob(fullfile(pth, ...
                sprintf("re%s_as_*proc-%s*-meanfield*.dtseries.nii", opts.measure, mdl.physio)));
            assert(~isemptytext(g1))
            cii1 = cifti_read(g1(1));
            g2 = mglob(fullfile(pth, ...
                sprintf("im%s_as_*proc-%s*-meanfield*.dtseries.nii", opts.measure, mdl.physio)));
            assert(~isemptytext(g2))
            cii2 = cifti_read(g2(1));
            A_nu = cii1.cdata' + 1i * cii2.cdata';
        end       

        function xi = reconstruct_xi(neg_dxi_dt, psi, opts)
            %% reconstructs complex bold signal xi from -dxi/dt, 
            %  matching initial time and dynamic range of bold signal psi

            arguments
                neg_dxi_dt {mustBeNumeric}
                psi {mustBeNumeric}
                opts.t_range = 200:896  % time frames
                opts.g_range = 1:2*32492  % cortical vertices
            end

            lambda_dxi_dt = iqr(abs(neg_dxi_dt(opts.t_range, opts.g_range)), "all");
            lambda_psi = iqr(abs(psi(opts.t_range, opts.g_range)), "all");
            lambda_dpsi_dt = iqr(abs(-diff(psi(opts.t_range, opts.g_range))), "all");

            psi1 = psi(1, :) / lambda_psi;
            neg_dxi_dt = neg_dxi_dt * lambda_dpsi_dt / lambda_dxi_dt;
            xi_ = -cumsum([psi1; neg_dxi_dt]);
            
            % final matching of ranges
            lambda_xi = iqr(abs(xi_(opts.t_range, opts.g_range)), "all");
            xi = xi_ * lambda_psi / lambda_xi;
        end
    
    end

    methods
        function this = mean_angle_instance(this, subs, opts)
            %% accepts text array subs; find mat files for subs; writes files indexed by opts.col_idx

            arguments
                this mlraut.AnalyticSignalHCPPar
                subs {mustBeText}
                opts.col_idx {mustBeScalarOrEmpty} = nan
                opts.new_physio {mustBeText} = "iFV"
                opts.anatomy {mustBeTextScalar} = "ctx"  % "cbm", "str", "thal"
                opts.test_range = []
                opts.tag {mustBeText} = ""
            end
            subs = convertCharsToStrings(subs);

            % DEBUG
            disp(subs)

            warning("off", "MATLAB:class:LoadDefinitionUpdated");

            % \emph{this} supplies utilities

            % abbrev.                 
            transform = @this.bin_by_physio_angle;
            nbin = this.twistors.num_bins_angles;

            % init
            ngo = this.num_nodes;
            angle_ = zeros(nbin, ngo);
            nmats_correct = 0;

            % find sub*.mat
            if ~isempty(opts.test_range)
                subs = subs(opts.test_range, :);
            end
            for sub = asrow(subs)
                if isemptytext(sub)
                    continue
                end
                mats = asrow(mglob( ...
                    fullfile(this.out_dir, sub, sprintf("sub-%s_ses-*ASHCPPar*.mat", sub))));
                if any(contains(mats, "-all"))
                    mats = mats(~contains(mats, "-all"));  % don't double count
                end
                if isemptytext(mats)
                    continue
                end
                errs = 0;

                for mat = mats
                    tic
                    try
                        ld = load(mat);
                        psi = ld.this_subset.bold_signal;
                        if isemptytext(opts.new_physio)
                            phi = ld.this_subset.physio_signal;
                        elseif strcmpi(opts.new_physio, ld.this_subset.source_physio)
                            phi = ld.this_subset.physio_signal;
                        elseif strcmpi(opts.new_physio, 'vis')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 1);
                        elseif strcmpi(opts.new_physio, 'sms')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 2);
                        elseif strcmpi(opts.new_physio, 'dan')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 3);
                        elseif strcmpi(opts.new_physio, 'van')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 4);
                        elseif strcmpi(opts.new_physio, 'lim')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 5);
                        elseif strcmpi(opts.new_physio, 'fpn')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 6);
                        elseif strcmpi(opts.new_physio, 'dmn')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 7);
                        else
                            re_phi = ld.this_subset.physio_supplementary(opts.new_physio);
                            assert(~isempty(re_phi))
                            assert(all(isfinite(re_phi)))
                            phi = hilbert(re_phi);
                        end

                        angle = transform(this.angle(psi, phi), phi);
                        angle_ = angle_ + angle;
                    catch ME
                        fprintf("while working with %s: ", mat);
                        handwarning(ME)
                        errs = errs + 1;
                    end
                    fprintf("time working with %s: ", mat);
                    toc
                end
                % adjust for errs
                nmats_correct = nmats_correct + numel(mats) - errs;
            end

            % complete averaging
            assert(nmats_correct > 0)
            angle_ = angle_/nmats_correct;

            dt = this.twistors_.theta_berry/this.twistors_.num_bins_angles;
            units_t = "RADIAN";
            ntag = sprintf('n%g', nmats_correct);
            if ~isemptytext(opts.new_physio)
                if ~strcmp(opts.anatomy, "ctx")
                    ptags = strrep(this.tags, this.source_physio, sprintf("%s-%s", opts.anatomy, opts.new_physio));
                else
                    ptags = strrep(this.tags, this.source_physio, opts.new_physio);
                end
            end
            if ~isemptytext(opts.tag)
                ptags = ptags + "-" + opts.tag;
            end
            this.cifti.write_cifti( ...
                angle_, sprintf('angle_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);

            warning("on", "MATLAB:class:LoadDefinitionUpdated");
        end

        function this = mean_twistor_instance(this, subs, opts)
            %% accepts text array subs; find mat files for subs; writes files indexed by opts.col_idx

            arguments
                this mlraut.AnalyticSignalHCPPar
                subs {mustBeText}
                opts.col_idx {mustBeScalarOrEmpty} = nan
                opts.new_physio {mustBeText} = "iFV"
                opts.anatomy {mustBeTextScalar} = "ctx"  % "cbm", "str", "thal"
                opts.test_range = []
                opts.tag {mustBeText} = ""

                % for this.Z_sup_alpha()
                opts.center_coord {mustBeNumeric} = mlraut.Twistors.LOCUS_CERULEUS_COORD  % mm, for precuneus
                opts.c {mustBeScalarOrEmpty} = this.lp_thresh  % speed limit for Poincare invariance ~ 0.1 ~ 1/timescale
                opts.use_E4 logical = true  % Ward & Wells sec. 8.1
                opts.s {mustBeScalarOrEmpty} = 0  % helicity in Penrose & Rindler eq. 6.2.7; only influences ~opts.use_E4
            end
            subs = convertCharsToStrings(subs);

            % DEBUG
            disp(subs)

            warning("off", "MATLAB:class:LoadDefinitionUpdated");

            % \emph{this} supplies utilities

            % abbrev.                 
            transform = @this.bin_by_physio_angle;
            nbin = this.twistors.num_bins_angles;

            % init
            ngo = this.num_nodes;
            rezeta_ = zeros(nbin, ngo);
            imzeta_ = zeros(nbin, ngo);
            reZ_sup_0_ = zeros(nbin, ngo);
            reZ_sup_1_ = zeros(nbin, ngo);
            reZ_sup_2_ = zeros(nbin, ngo);
            imZ_sup_0_ = zeros(nbin, ngo);
            imZ_sup_1_ = zeros(nbin, ngo);
            imZ_sup_2_ = zeros(nbin, ngo);
            X_ = zeros(nbin, ngo);
            Y_ = zeros(nbin, ngo);
            Z_ = zeros(nbin, ngo);
            T_ = zeros(nbin, ngo);
            r_ = zeros(1, ngo);
            bold_ = zeros(nbin, ngo);
            % neg_dbold_dt_ = zeros(nbin, ngo);
            replvs_ = zeros(nbin, ngo);
            implvs_ = zeros(nbin, ngo);
            nmats_corr = 0;

            % find sub*.mat
            if ~isempty(opts.test_range)
                subs = subs(opts.test_range, :);
            end
            for sub = asrow(subs)
                if isemptytext(sub)
                    continue
                end
                mats = asrow(mglob( ...
                    fullfile(this.out_dir, sub, sprintf("sub-%s_ses-*ASHCPPar*.mat", sub))));
                if any(contains(mats, "-all"))
                    mats = mats(~contains(mats, "-all"));  % don't double count
                end
                if isemptytext(mats)
                    continue
                end
                errs = 0;

                for mat = mats
                    tic
                    try
                        ld = load(mat);
                        psi = ld.this_subset.bold_signal;
                        if isemptytext(opts.new_physio)
                            phi = ld.this_subset.physio_signal;
                        elseif strcmpi(opts.new_physio, ld.this_subset.source_physio)
                            phi = ld.this_subset.physio_signal;
                        elseif strcmpi(opts.new_physio, 'vis')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 1);
                        elseif strcmpi(opts.new_physio, 'sms')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 2);
                        elseif strcmpi(opts.new_physio, 'dan')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 3);
                        elseif strcmpi(opts.new_physio, 'van')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 4);
                        elseif strcmpi(opts.new_physio, 'lim')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 5);
                        elseif strcmpi(opts.new_physio, 'fpn')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 6);
                        elseif strcmpi(opts.new_physio, 'dmn')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 7);
                        else
                            re_phi = ld.this_subset.physio_supplementary(opts.new_physio);
                            assert(~isempty(re_phi))
                            assert(all(isfinite(re_phi)))
                            phi = hilbert(re_phi);
                        end

                        zeta = transform(this.zeta(psi, phi), phi);
                        Z_sup_alpha = cellfun(@(x) transform(x, phi), ...
                            this.Z_sup_alpha( ...
                                psi, phi, ...
                                center_coord=opts.center_coord, ...
                                c=opts.c, ...
                                use_E4=opts.use_E4, ...
                                s=opts.s), ...
                            UniformOutput=false);
                        X = transform(this.X(psi, phi), phi);
                        Y = transform(this.Y(psi, phi), phi);
                        Z = transform(this.Z(psi, phi), phi);
                        T = transform(this.T(psi, phi), phi);
                        r = asrow(this.connectivity(psi, phi));
                        bold = transform(psi, phi);
                        % neg_dbold_dt = transform(this.neg_dbold_dt(psi, phi), phi);
                        plvs = transform(this.phase_locked_values(psi, phi), phi);
                        
                        rezeta_ = rezeta_ + real(zeta);
                        imzeta_ = imzeta_ + imag(zeta);
                        reZ_sup_0_ = reZ_sup_0_ + real(Z_sup_alpha{1});
                        reZ_sup_1_ = reZ_sup_1_ + real(Z_sup_alpha{2});
                        reZ_sup_2_ = reZ_sup_2_ + real(Z_sup_alpha{3});
                        imZ_sup_0_ = imZ_sup_0_ + imag(Z_sup_alpha{1});
                        imZ_sup_1_ = imZ_sup_1_ + imag(Z_sup_alpha{2});
                        imZ_sup_2_ = imZ_sup_2_ + imag(Z_sup_alpha{3});
                        X_ = X_ + real(X);
                        Y_ = Y_ + real(Y);
                        Z_ = Z_ + real(Z);
                        T_ = T_ + real(T);
                        r_ = r_ + real(r);
                        bold_ = bold_ + real(bold);
                        % neg_dbold_dt_ = neg_dbold_dt_ + real(neg_dbold_dt);
                        replvs_ = replvs_ + real(plvs);
                        implvs_ = implvs_ + imag(plvs);
                    catch ME
                        fprintf("while working with %s: ", mat);
                        handwarning(ME)
                        errs = errs + 1;
                    end
                    fprintf("time working with %s: ", mat);
                    toc
                end
                % adjust for errs
                nmats_corr = nmats_corr + numel(mats) - errs;
            end

            % complete averaging
            assert(nmats_corr > 0)
            rezeta_ = rezeta_/nmats_corr;
            imzeta_ = imzeta_/nmats_corr;
            reZ_sup_0_ = reZ_sup_0_/nmats_corr;
            reZ_sup_1_ = reZ_sup_1_/nmats_corr;
            reZ_sup_2_ = reZ_sup_2_/nmats_corr;
            imZ_sup_0_ = imZ_sup_0_/nmats_corr;
            imZ_sup_1_ = imZ_sup_1_/nmats_corr;
            imZ_sup_2_ = imZ_sup_2_/nmats_corr;
            X_ = X_/nmats_corr;
            Y_ = Y_/nmats_corr;
            Z_ = Z_/nmats_corr;
            T_ = T_/nmats_corr;
            r_ = r_/nmats_corr;
            bold_ = bold_/nmats_corr;
            % neg_dbold_dt_ = neg_dbold_dt_/nmats_corr;
            replvs_ = replvs_/nmats_corr;
            implvs_ = implvs_/nmats_corr;

            dt = this.twistors_.theta_berry/this.twistors_.num_bins_angles;
            units_t = "RADIAN";
            ntag = sprintf('n%g', nmats_corr);
            if ~isemptytext(opts.new_physio)
                if ~strcmp(opts.anatomy, "ctx")
                    ptags = strrep(this.tags, this.source_physio, sprintf("%s-%s", opts.anatomy, opts.new_physio));
                else
                    ptags = strrep(this.tags, this.source_physio, opts.new_physio);
                end
            end
            if ~isemptytext(opts.tag)
                ptags = ptags + "-" + opts.tag;
            end
            this.cifti.write_cifti( ...
                rezeta_, sprintf('rezeta_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                imzeta_, sprintf('imzeta_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                reZ_sup_0_, sprintf('reZsup0_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                reZ_sup_1_, sprintf('reZsup1_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                reZ_sup_2_, sprintf('reZsup2_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                imZ_sup_0_, sprintf('imZsup0_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                imZ_sup_1_, sprintf('imZsup1_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                imZ_sup_2_, sprintf('imZsup2_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                X_, sprintf('X_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                Y_, sprintf('Y_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                Z_, sprintf('Z_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                T_, sprintf('T_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                r_, sprintf('comparator_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                bold_, sprintf('bold_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            % this.cifti.write_cifti( ...
            %     neg_dbold_dt_, sprintf('neg-dbold-dt_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                replvs_, sprintf('replvs_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                implvs_, sprintf('implvs_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);

            warning("on", "MATLAB:class:LoadDefinitionUpdated");
        end

        function this = mean_Z_sup_alpha_instance(this, subs, opts)
            %% accepts text array subs; find mat files for subs; writes files indexed by opts.col_idx

            arguments
                this mlraut.AnalyticSignalHCPPar
                subs {mustBeText}
                opts.col_idx {mustBeScalarOrEmpty} = nan
                opts.new_physio {mustBeText} = "iFV"
                opts.anatomy {mustBeTextScalar} = "ctx"  % "cbm", "str", "thal"
                opts.test_range = []
                opts.tag {mustBeText} = ""

                % for this.Z_sup_alpha()
                opts.center_coord {mustBeNumeric} = mlraut.Twistors.LOCUS_CERULEUS_COORD  % mm, for precuneus
                opts.c {mustBeScalarOrEmpty} = this.lp_thresh  % speed limit for Poincare invariance ~ 0.1 ~ 1/timescale
                opts.use_E4 logical = false  % Ward & Wells sec. 8.1
                opts.s {mustBeScalarOrEmpty} = 0  % helicity in Penrose & Rindler eq. 6.2.7; only influences ~opts.use_E4
            end
            subs = convertCharsToStrings(subs);

            % DEBUG
            disp(subs)

            warning("off", "MATLAB:class:LoadDefinitionUpdated");

            % \emph{this} supplies utilities

            % abbrev.                 
            transform = @this.bin_by_physio_angle;
            nbin = this.twistors.num_bins_angles;

            % init
            ngo = this.num_nodes;
            reZ_sup_0_ = zeros(nbin, ngo);
            reZ_sup_1_ = zeros(nbin, ngo);
            reZ_sup_2_ = zeros(nbin, ngo);
            imZ_sup_0_ = zeros(nbin, ngo);
            imZ_sup_1_ = zeros(nbin, ngo);
            imZ_sup_2_ = zeros(nbin, ngo);
            nmats_corr = 0;

            % find sub*.mat
            if ~isempty(opts.test_range)
                subs = subs(opts.test_range, :);
            end
            for sub = asrow(subs)
                if isemptytext(sub)
                    continue
                end
                mats = asrow(mglob( ...
                    fullfile(this.out_dir, sub, sprintf("sub-%s_ses-*ASHCPPar*.mat", sub))));
                if any(contains(mats, "-all"))
                    mats = mats(~contains(mats, "-all"));  % don't double count
                end
                if isemptytext(mats)
                    continue
                end
                errs = 0;

                for mat = mats
                    tic
                    try
                        ld = load(mat);
                        psi = ld.this_subset.bold_signal;
                        if isemptytext(opts.new_physio)
                            phi = ld.this_subset.physio_signal;
                        elseif strcmpi(opts.new_physio, ld.this_subset.source_physio)
                            phi = ld.this_subset.physio_signal;
                        elseif strcmpi(opts.new_physio, 'vis')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 1);
                        elseif strcmpi(opts.new_physio, 'sms')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 2);
                        elseif strcmpi(opts.new_physio, 'dan')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 3);
                        elseif strcmpi(opts.new_physio, 'van')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 4);
                        elseif strcmpi(opts.new_physio, 'lim')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 5);
                        elseif strcmpi(opts.new_physio, 'fpn')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 6);
                        elseif strcmpi(opts.new_physio, 'dmn')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 7);
                        else
                            re_phi = ld.this_subset.physio_supplementary(opts.new_physio);
                            assert(~isempty(re_phi))
                            assert(all(isfinite(re_phi)))
                            phi = hilbert(re_phi);
                        end

                        Z_sup_alpha = cellfun(@(x) transform(x, phi), ...
                            this.Z_sup_alpha( ...
                                psi, phi, ...
                                center_coord=opts.center_coord, ...
                                c=opts.c, ...
                                use_E4=opts.use_E4, ...
                                s=opts.s), ...
                            UniformOutput=false);
                        
                        reZ_sup_0_ = reZ_sup_0_ + real(Z_sup_alpha{1});
                        reZ_sup_1_ = reZ_sup_1_ + real(Z_sup_alpha{2});
                        reZ_sup_2_ = reZ_sup_2_ + real(Z_sup_alpha{3});
                        imZ_sup_0_ = imZ_sup_0_ + imag(Z_sup_alpha{1});
                        imZ_sup_1_ = imZ_sup_1_ + imag(Z_sup_alpha{2});
                        imZ_sup_2_ = imZ_sup_2_ + imag(Z_sup_alpha{3});
                    catch ME
                        fprintf("while working with %s: ", mat);
                        handwarning(ME)
                        errs = errs + 1;
                    end
                    fprintf("time working with %s: ", mat);
                    toc
                end
                % adjust for errs
                nmats_corr = nmats_corr + numel(mats) - errs;
            end

            % complete averaging
            assert(nmats_corr > 0)
            reZ_sup_0_ = reZ_sup_0_/nmats_corr;
            reZ_sup_1_ = reZ_sup_1_/nmats_corr;
            reZ_sup_2_ = reZ_sup_2_/nmats_corr;
            imZ_sup_0_ = imZ_sup_0_/nmats_corr;
            imZ_sup_1_ = imZ_sup_1_/nmats_corr;
            imZ_sup_2_ = imZ_sup_2_/nmats_corr;

            dt = this.twistors_.theta_berry/this.twistors_.num_bins_angles;
            units_t = "RADIAN";
            ntag = sprintf('n%g', nmats_corr);
            if ~isemptytext(opts.new_physio)
                if ~strcmp(opts.anatomy, "ctx")
                    ptags = strrep(this.tags, this.source_physio, sprintf("%s-%s", opts.anatomy, opts.new_physio));
                else
                    ptags = strrep(this.tags, this.source_physio, opts.new_physio);
                end
            end
            if ~isemptytext(opts.tag)
                ptags = ptags + "-" + opts.tag + ...
                    strrep(sprintf("-coord%s", mat2str(opts.center_coord)), " ", ",");
            end
            this.cifti.write_cifti( ...
                reZ_sup_0_, sprintf('reZsup0_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                reZ_sup_1_, sprintf('reZsup1_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                reZ_sup_2_, sprintf('reZsup2_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                imZ_sup_0_, sprintf('imZsup0_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                imZ_sup_1_, sprintf('imZsup1_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                imZ_sup_2_, sprintf('imZsup2_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);            

            warning("on", "MATLAB:class:LoadDefinitionUpdated");
        end
        
        function this = mean_transformed_twistor_instance(this, subs, opts)
            %% Fourier transformed.
            %  Accepts text array subs; find mat files for subs; writes files indexed by opts.col_idx

            arguments
                this mlraut.AnalyticSignalHCPPar
                subs {mustBeText}
                opts.col_idx {mustBeInteger}
                opts.new_physio {mustBeText} = "iFV"
                opts.test_range = []
                opts.tag {mustBeText} = ""
                opts.anatomy {mustBeTextScalar} = "ctx"  % "cbm", "str", "thal"
            end
            subs = convertCharsToStrings(subs);

            % DEBUG
            disp(subs)

            warning("off", "MATLAB:class:LoadDefinitionUpdated");

            % \emph{this} supplies utilities

            % abbrev.
            switch char(opts.tag)
                case 'fft'
                    transform = @(x,y) fft(x);
                    nbin = this.num_frames - 1;
                case 'ifft'
                    transform = @(x,y) ifft(x);
                    nbin = this.num_frames - 1;
                case 'radon'
                    % https://www.mathworks.com/help/releases/R2024b/images/radon-transform.html?searchHighlight=radon&s_tid=doc_srchtitle                    
                    nbin = this.num_frames - 1;
                    theta = 0:180/nbin; 
                    transform = @(I_) radon(I_, theta);
                otherwise                    
                    error("mlraut:ValueError", stackstr())
            end

            % init
            ngo = this.num_nodes;
            reX_ = zeros(nbin, ngo);
            imX_ = zeros(nbin, ngo);
            reY_ = zeros(nbin, ngo);
            imY_ = zeros(nbin, ngo);
            reZ_ = zeros(nbin, ngo);
            imZ_ = zeros(nbin, ngo);
            reT_ = zeros(nbin, ngo);
            imT_ = zeros(nbin, ngo);
            nmats_corr = 0;

            % find sub*.mat
            if ~isempty(opts.test_range)
                subs = subs(opts.test_range, :);
            end
            for sub = asrow(subs)
                if isemptytext(sub)
                    continue
                end
                mats = asrow(mglob( ...
                    fullfile(this.out_dir, sub, sprintf("sub-%s_ses-*ASHCPPar*.mat", sub))));
                if isemptytext(mats)
                    continue
                end
                errs = 0;

                for mat = mats
                    tic
                    try
                        ld = load(mat);
                        psi = ld.this_subset.bold_signal;
                        if isemptytext(opts.new_physio)
                            phi = ld.this_subset.physio_signal;
                        elseif strcmpi(opts.new_physio, ld.this_subset.source_physio)
                            phi = ld.this_subset.physio_signal;
                        elseif strcmpi(opts.new_physio, 'vis')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 1);
                        elseif strcmpi(opts.new_physio, 'sms')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 2);
                        elseif strcmpi(opts.new_physio, 'dan')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 3);
                        elseif strcmpi(opts.new_physio, 'van')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 4);
                        elseif strcmpi(opts.new_physio, 'lim')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 5);
                        elseif strcmpi(opts.new_physio, 'fpn')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 6);
                        elseif strcmpi(opts.new_physio, 'dmn')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 7);
                        else
                            re_phi = ld.this_subset.physio_supplementary(opts.new_physio);
                            assert(~isempty(re_phi))
                            assert(all(isfinite(re_phi)))
                            phi = hilbert(re_phi);
                        end
                        X = transform(this.X(psi, phi), phi);
                        Y = transform(this.Y(psi, phi), phi);
                        Z = transform(this.Z(psi, phi), phi);
                        T = transform(this.T(psi, phi), phi);

                        reX_ = reX_ + real(X);
                        imX_ = imX_ + imag(X);
                        reY_ = reY_ + real(Y);
                        imY_ = imY_ + imag(Y);
                        reZ_ = reZ_ + real(Z);
                        imZ_ = imZ_ + imag(Z);
                        reT_ = reT_ + real(T);
                        imT_ = imT_ + imag(T);
                    catch ME
                        fprintf("while working with %s: ", mat);
                        handwarning(ME)
                        errs = errs + 1;
                    end
                    fprintf("time working with %s: ", mat);
                    toc
                end
                % adjust for errs
                nmats_corr = nmats_corr + numel(mats) - errs;
            end

            % complete averaging
            reX_ = reX_/nmats_corr;
            imX_ = imX_/nmats_corr;
            reY_ = reY_/nmats_corr;
            imY_ = imY_/nmats_corr;
            reZ_ = reZ_/nmats_corr;
            imZ_ = imZ_/nmats_corr;
            reT_ = reT_/nmats_corr;
            imT_ = imT_/nmats_corr;

            df = this.Fs/(this.num_frames - 1);
            units_f = "HERTZ";
            ntag = sprintf('n%g', nmats_corr);
            if ~isemptytext(opts.new_physio)
                ptags = strrep(this.tags, this.source_physio, opts.new_physio);
            end
            if ~isemptytext(opts.tag)
                ptags = strrep(ptags, "-ASHCPPar", sprintf("-%s-ASHCPPar", opts.tag));
            end
            this.cifti.write_cifti( ...
                reX_, sprintf('reX_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=df, units_t=units_f);
            this.cifti.write_cifti( ...
                imX_, sprintf('imX_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=df, units_t=units_f);
            this.cifti.write_cifti( ...
                reY_, sprintf('reY_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=df, units_t=units_f);
            this.cifti.write_cifti( ...
                imY_, sprintf('imY_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=df, units_t=units_f);
            this.cifti.write_cifti( ...
                reZ_, sprintf('reZ_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=df, units_t=units_f);
            this.cifti.write_cifti( ...
                imZ_, sprintf('imZ_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=df, units_t=units_f);
            this.cifti.write_cifti( ...
                reT_, sprintf('reT_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=df, units_t=units_f);
            this.cifti.write_cifti( ...
                imT_, sprintf('imT_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=df, units_t=units_f);

            warning("on", "MATLAB:class:LoadDefinitionUpdated");
        end

        function this = mean_transformed_zeta(this, subs, opts)
            %% Fourier transformed.
            %  Accepts text array subs; find mat files for subs; writes files indexed by opts.col_idx

            arguments
                this mlraut.AnalyticSignalHCPPar
                subs {mustBeText}
                opts.col_idx {mustBeScalarOrEmpty} = []
                opts.new_physio {mustBeText} = ""
                opts.test_range = []
                opts.tag {mustBeText} = ""
                opts.rsn {mustBeScalarOrEmpty} = []
                opts.anatomy {mustBeTextScalar} = "ctx"  % "cbm", "str", "thal"
            end
            subs = convertCharsToStrings(subs);
            if isscalar(subs)
                out_dir_ = fullfile(this.out_dir, subs);
            else
                out_dir_ = this.out_dir;
            end

            % DEBUG
            disp(subs)

            warning("off", "MATLAB:class:LoadDefinitionUpdated");

            % \emph{this} supplies utilities

            % abbrev.
            switch char(opts.tag)
                case 'fft'
                    transform = @(x) fft(x);
                    nbin = this.num_frames - 1;
                case 'ifft'
                    transform = @(x) ifft(x);
                    nbin = this.num_frames - 1;
                case 'radon'
                    % https://www.mathworks.com/help/releases/R2024b/images/radon-transform.html?searchHighlight=radon&s_tid=doc_srchtitle                    
                    nbin = this.num_frames - 1;
                    theta = 0:180/nbin; 
                    transform = @(I_) radon(I_, theta);
                otherwise                    
                    transform = @(x) x;  % identity
            end

            % init
            ngo = this.num_nodes;
            zeta_nu_ = zeros(nbin, ngo);
            zeta_t_avggo_ = zeros(nbin, 1);
            nmats_corrected = 0;
            if isempty(opts.rsn)
                mask_rsn = [];
            else
                ld = load("~/MATLAB-Drive/arousal-waves-main/supporting_files/networks_HCP.mat");
                mask_rsn = ld.assns2 == opts.rsn;
            end

            % find sub*.mat
            if ~isempty(opts.test_range)
                subs = subs(opts.test_range, :);
            end
            for sub = asrow(subs)
                if isemptytext(sub)
                    continue
                end
                mats = asrow(mglob( ...
                    fullfile(this.out_dir, sub, sprintf("sub-%s_ses-*ASHCPPar*.mat", sub))));
                if isemptytext(mats)
                    continue
                end
                errs = 0;

                for mat = mats
                    tic
                    try
                        ld = load(mat);  % this provides methods, ld.this_subset provides data
                        psi = ld.this_subset.bold_signal;

                        if ~isempty(mask_rsn)
                            psi = psi(:, mask_rsn);
                        end

                        % assign phi
                        if isemptytext(opts.new_physio) || strcmp(opts.new_physio, ld.this_subset.source_physio)
                            phi = ld.this_subset.physio_signal;
                        elseif strcmpi(opts.new_physio, 'vis')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 1);
                        elseif strcmpi(opts.new_physio, 'sms')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 2);
                        elseif strcmpi(opts.new_physio, 'dan')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 3);
                        elseif strcmpi(opts.new_physio, 'van')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 4);
                        elseif strcmpi(opts.new_physio, 'lim')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 5);
                        elseif strcmpi(opts.new_physio, 'fpn')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 6);
                        elseif strcmpi(opts.new_physio, 'dmn')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 7);
                        else
                            re_phi = ld.this_subset.physio_supplementary(opts.new_physio);
                            assert(~isempty(re_phi))
                            assert(all(isfinite(re_phi)))
                            phi = hilbert(re_phi);
                        end
                        zeta_t = this.zeta(psi, phi);
                        zeta_t_avggo = mean(zeta_t, 2);
                        zeta_t_avggo_ = zeta_t_avggo_ + zeta_t_avggo;
                        zeta_nu = transform(zeta_t);
                        zeta_nu_ = zeta_nu_ + zeta_nu;  
                    catch ME
                        fprintf("while working with %s: ", mat);
                        handwarning(ME)
                        errs = errs + 1;
                    end
                    fprintf("time working with %s: ", mat);
                    toc
                end

                % adjust for errs
                nmats_corrected = nmats_corrected + numel(mats) - errs;
            end

            % complete averaging
            zeta_nu_ = zeta_nu_/nmats_corrected;
            zeta_t_avggo_ = zeta_t_avggo_/nmats_corrected;

            % construct identifiers
            n_tag = sprintf('n%g', nmats_corrected);
            new_tags = this.tags;
            if ~isemptytext(opts.new_physio)
                new_tags = strrep(this.tags, this.source_physio, opts.new_physio);
            end
            if ~isempty(opts.rsn)
                new_tags = strrep(new_tags, "-ASHCPPar", sprintf("-rsn%i-ASHCPPar", opts.rsn));
            end
            if ~isempty(opts.col_idx)
                par_tag = "_par" + opts.col_idx;
            else
                par_tag = "";
            end

            % save mat
            fqfn_zeta = fullfile(out_dir_, ...
                sprintf('%szeta_as_sub-%s_ses-%s_%s%s.mat', opts.tag, n_tag, n_tag, new_tags, par_tag));
            save(fqfn_zeta, "zeta_nu_", "-v7.3");
            fqfn_zetaavggo = fullfile(out_dir_, ...
                sprintf('zetaavggo_as_sub-%s_ses-%s_%s%s.mat', n_tag, n_tag, new_tags, par_tag));
            save(fqfn_zetaavggo, "zeta_t_avggo_", "-v7.3");

            % save cifti
            if isempty(mask_rsn)
                if strcmp(opts.tag, "fft")
                    dt = this.Fs/(this.num_frames - 1);
                    units_t = "HERTZ";
                else
                    dt = this.tr;
                    units_t = "SECONDS";
                end
                this.cifti.write_cifti( ...
                    real(zeta_nu_), ...
                    fullfile(out_dir_, ...
                        sprintf('re%szeta_as_sub-%s_ses-%s_%s%s', ...
                        opts.tag, n_tag, n_tag, new_tags, par_tag)), ...
                    dt=dt, units_t=units_t);
                this.cifti.write_cifti( ...
                    imag(zeta_nu_), ...
                    fullfile(out_dir_, ...
                        sprintf('im%szeta_as_sub-%s_ses-%s_%s%s', ...
                        opts.tag, n_tag, n_tag, new_tags, par_tag)), ...
                    dt=dt, units_t=units_t);
            end
            
            warning("on", "MATLAB:class:LoadDefinitionUpdated");
        end

        function h = plot_coherencyc(this, psi_multi, phi_multi, opts)
            arguments
                this mlraut.AnalyticSignalHCPPar
                psi_multi = []
                phi_multi = []
                opts.tseries {mustBeTextScalar} = "-dbold/dt"
                opts.physio_tag {mustBeTextScalar} = this.source_physio
                opts.physio_tag2 {mustBeTextScalar} = this.source_physio
            end

            warning("off", "MATLAB:plot:IgnoreImaginaryXYPart")

            % build fultz; quick if {psi,phi}_multi are supplied
            fultz = mlraut.FultzMulti( ...
                this, ...
                psi_multi=psi_multi, ...
                phi_multi=phi_multi);
            if isempty(fultz.phi_multi) || isempty(fultz.psi_multi)
                fultz.build_tseries_phi_psi_from_mat( ...
                    mat_pattern=fullfile(fileparts(this.out_dir), "**", ...
                    sprintf("sub-*_ses-*_proc-%s-*-ASHCPPar.mat", this.source_physio)));
            end

            % plot
            pwd0 = pushd(fultz.out_dir);
            h = figure;
            fultz.plot_coherencyc( ...
                tseries=opts.tseries, ...
                physio_tag=opts.physio_tag, ...
                physio_tag2=opts.physio_tag2);
            saveFigure2(h, fultz.fqfileprefix + "_plot_coherencyc");
            popd(pwd0);

            warning("on", "MATLAB:plot:IgnoreImaginaryXYPart")
        end

        function h = rugplot(this, psi, opts)
            arguments
                this mlraut.AnalyticSignalHCPPar
                psi {mustBeNumeric}
                opts.fig_fileprefix {mustBeTextScalar} = ""
                opts.is_spectral logical = false
                opts.closeFigure logical = true
                opts.title {mustBeTextScalar} = ""
                opts.complex_representation function_handle = @real
            end
            repr = opts.complex_representation;  % abbrev.
            [pth,fp,x] = fileparts(opts.fig_fileprefix);
            if ~startsWith(fp, func2str(repr))
                opts.fig_fileprefix = fullfile( ...
                    pth, ...
                    string(func2str(repr)) + string(fp) + string(x));
            end

            nd = mlraut.NetworkData(this, psi);
            img = nd.reshape_by_anatomy(psi);
            if opts.is_spectral
                img = fftshift(img, 1);
            end

            % fig position
            fig_width = 2056;
            fig_height = 1329;
            h = figure('Position', [1, 1, fig_width, fig_height]);
            movegui(h, 'northwest');  % Moves to upper left

            if opts.is_spectral
                imagesc(repr(img));
                colorbar;
                if strcmp(func2str(repr), "abs")
                    clim([0, 1000]);
                else
                    clim([-150, 150]);
                end
            else
                imagesc(repr(img));
                colorbar;
                clim([-5, 5]);
            end

            hold on

            % colors
            if strcmp(func2str(repr), "abs")
                color = [.8 .8 .8];
            else
                color = [.2 .2 .2];
            end

            % lines
            xline(8788, 'Color', color, 'LineWidth', 1);  % VIS ctx
            xline(20748, 'Color', color, 'LineWidth', 1);  % SMS
            xline(27510, 'Color', color, 'LineWidth', 1);  % DAN
            xline(34683, 'Color', color, 'LineWidth', 1);  % VAN
            xline(39219, 'Color', color, 'LineWidth', 1);  % LIM
            xline(46530, 'Color', color, 'LineWidth', 1);  % FPN
            xline(58666, 'Color', color, 'LineWidth', 1);  % DMN            
            xline(76498, 'Color', color, 'LineWidth', 1);  % cbm
            xline(79855, 'Color', color, 'LineWidth', 1);  % str
            % xline(82053, 'Color', color, 'LineWidth', 1);  % thal

            % annotation
            apos = 0.95*this.num_frames;
            text(8788-1200, apos, 'VIS', 'FontSize', 9, 'FontAngle', 'italic', 'Color', color, 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(20748-1200, apos, 'SMS', 'FontSize', 9, 'FontAngle', 'italic', 'Color', color, 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(27510-1200, apos, 'DAN', 'FontSize', 9, 'FontAngle', 'italic', 'Color', color, 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(34683-1200, apos, 'VAN', 'FontSize', 9, 'FontAngle', 'italic', 'Color', color, 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(39219-1200, apos, 'LIM', 'FontSize', 9, 'FontAngle', 'italic', 'Color', color, 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(46530-1200, apos, 'FPN', 'FontSize', 9, 'FontAngle', 'italic', 'Color', color, 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(58666-1200, apos, 'DMN', 'FontSize', 9, 'FontAngle', 'italic', 'Color', color, 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(76498-1200, apos, 'cerebellum', 'FontSize', 9, 'FontAngle', 'italic', 'Color', color, 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(79855-1200, apos, 'striatum', 'FontSize', 9, 'FontAngle', 'italic', 'Color', color, 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(82053-1200, apos, 'thalamus', 'FontSize', 9, 'FontAngle', 'italic', 'Color', color, 'HorizontalAlignment', 'left', 'Rotation', 90);

            hold off

            if opts.is_spectral

                % Calculate time values
                Nomega = size(img, 1);
                Omega = (-Nomega/2:Nomega/2-1) * 2 * pi * this.Fs / Nomega; % total freq. range in radian/s
                % time_values = linspace(0, T, Nt);

                n_ticks = 11; % Choose how many tick marks you want
                tick_indices = linspace(1, Nomega, n_ticks);
                tick_omega = linspace(Omega(1), Omega(end), n_ticks);

                % indices -> times (s)
                ax = gca;
                ax.YTick = tick_indices;
                % set(ax, ...
                %     'YTickLabel', arrayfun(@(x) sprintf("%.1f", x), tick_omega/pi, 'UniformOutput', false) + "\pi", ...
                %     'TickLabelInterpreter', 'latex');
                ax.YTickLabel = arrayfun(@(x) sprintf('%.1f', x), tick_omega, 'UniformOutput', false);

                % Add labels
                switch func2str(repr)
                    case 'real'
                        atitle = "$$\Re \, \zeta(\mathbf{x}, \omega)$$";
                    case 'imag'
                        atitle = "$$\Im \, \zeta(\mathbf{x}, \omega)$$";
                    case 'abs'
                        atitle = "$$\left| \zeta(\mathbf{x}, \omega) \right|$$";
                    otherwise
                        atitle = "$$\zeta(\mathbf{x}, \omega)$$";
                end

                xlabel("$\mathbf{x}$ (greyordinates)", Interpreter="latex")
                ylabel("$\omega$ (rad/s)", Interpreter="latex")
                if isemptytext(opts.title)
                    title(atitle, Interpreter="latex")
                else
                    title(opts.title, Interpreter="latex")
                end

                fontsize(scale=3)
            else

                % Calculate time values
                Nt = size(img, 1);
                T = Nt * this.tr; % total time in seconds
                % time_values = linspace(0, T, Nt);

                n_ticks = 11; % Choose how many tick marks you want
                tick_indices = linspace(1, Nt, n_ticks);
                tick_times = linspace(0, T, n_ticks);

                % indices -> times (s)
                ax = gca;
                ax.YTick = tick_indices;
                ax.YTickLabel = arrayfun(@(x) sprintf('%.0f', x), tick_times, 'UniformOutput', false);

                % Add labels
                switch func2str(repr)
                    case 'real'
                        atitle = "$$\Re \, \zeta(\mathbf{x}, t)$$";
                    case 'imag'
                        atitle = "$$\Im \, \zeta(\mathbf{x}, t)$$";
                    case 'abs'
                        atitle = "$$\left| \zeta(\mathbf{x}, t) \right|$$";
                    otherwise
                        atitle = "$$\zeta(\mathbf{x}, t)$$";
                end

                xlabel("$\mathbf{x}$ (greyordinates)", Interpreter="latex")
                ylabel("time (s)", Interpreter="latex")
                if isemptytext(opts.title)
                    title(atitle, Interpreter="latex")
                else
                    title(opts.title, Interpreter="latex")
                end

                fontsize(scale=3)
            end

            % save figure
            if ~isemptytext(opts.fig_fileprefix)
                saveFigure2(h, opts.fig_fileprefix, ext={'.png'}, closeFigure=opts.closeFigure);
            end
        end

        function h = rugplot_angles(this, psi, phi, opts)
            arguments
                this mlraut.AnalyticSignalHCPPar
                psi {mustBeNumeric}
                phi {mustBeNumeric}
                opts.fig_fileprefix {mustBeTextScalar} = ""
                opts.tag {mustBeTextScalar} = "psi"
                opts.is_spectral logical = false
                opts.closeFigure logical = true
                opts.title {mustBeTextScalar} = ""
                opts.complex_representation function_handle = @real
            end
            if isemptytext(opts.fig_fileprefix)
                opts.fig_fileprefix = extractBefore(this.mat_fqfn, "-ASHCPPar") + "-ASHCPPar-rugplot-angles-" + opts.tag;
            end
            repr = opts.complex_representation;  % abbrev.

            % fig position
            fig_width = 2056;
            fig_height = fig_width;  % 1329;
            h = figure('Position', [1, 1, fig_width, fig_height]);
            movegui(h, 'northwest');  % Moves to upper left

            % generate imagesc
            binned_psi = this.bin_by_physio_angle(psi, phi, num_bins_angles=800);
            nd = mlraut.NetworkData(this, binned_psi);
            img = nd.reshape_by_anatomy(binned_psi);
            imagesc(repr(img));
            colormap(ry_bc_bl_no_white(1024));
            colorbar;
            clim([-5, 5]);

            hold on

            % colors
            color = [.5 .5 .5];

            % lines
            xline(8788, 'Color', color, 'LineWidth', 1);  % VIS ctx
            xline(20748, 'Color', color, 'LineWidth', 1);  % SMS
            xline(27510, 'Color', color, 'LineWidth', 1);  % DAN
            xline(34683, 'Color', color, 'LineWidth', 1);  % VAN
            xline(39219, 'Color', color, 'LineWidth', 1);  % LIM
            xline(46530, 'Color', color, 'LineWidth', 1);  % FPN
            xline(58666, 'Color', color, 'LineWidth', 1);  % DMN            
            xline(76498, 'Color', color, 'LineWidth', 1);  % cbm
            xline(79855, 'Color', color, 'LineWidth', 1);  % str
            % xline(82053, 'Color', color, 'LineWidth', 1);  % thal

            % annotation
            apos = 0.95*size(img, 1);
            text(8788-1200, apos, 'VIS', 'FontSize', 9, 'FontAngle', 'italic', 'Color', color, 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(20748-1200, apos, 'SMS', 'FontSize', 9, 'FontAngle', 'italic', 'Color', color, 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(27510-1200, apos, 'DAN', 'FontSize', 9, 'FontAngle', 'italic', 'Color', color, 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(34683-1200, apos, 'VAN', 'FontSize', 9, 'FontAngle', 'italic', 'Color', color, 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(39219-1200, apos, 'LIM', 'FontSize', 9, 'FontAngle', 'italic', 'Color', color, 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(46530-1200, apos, 'FPN', 'FontSize', 9, 'FontAngle', 'italic', 'Color', color, 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(58666-1200, apos, 'DMN', 'FontSize', 9, 'FontAngle', 'italic', 'Color', color, 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(76498-1200, apos, 'cerebellum', 'FontSize', 9, 'FontAngle', 'italic', 'Color', color, 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(79855-1200, apos, 'striatum', 'FontSize', 9, 'FontAngle', 'italic', 'Color', color, 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(82053-1200, apos, 'thalamus', 'FontSize', 9, 'FontAngle', 'italic', 'Color', color, 'HorizontalAlignment', 'left', 'Rotation', 90);

            hold off

            % Calculate angle values
            Nangles = size(img, 1);
            n_ticks = 9; % Choose how many tick marks you want
            tick_indices = linspace(1, Nangles, n_ticks);
            tick_angles = linspace(0, this.twistors_.theta_berry/pi, n_ticks);

            % indices -> angles (rad)
            ax = gca;
            ax.YTick = tick_indices;
            ax.YTickLabel = arrayfun(@(x) sprintf('%.1f', x), tick_angles, 'UniformOutput', false);

            % Add labels
            switch func2str(repr)
                case 'real'
                    atitle = "$$\Re \, \zeta(\mathbf{x}, t)$$";
                case 'imag'
                    atitle = "$$\Im \, \zeta(\mathbf{x}, t)$$";
                case 'abs'
                    atitle = "$$\left| \zeta(\mathbf{x}, t) \right|$$";
                otherwise
                    atitle = "$$\zeta(\mathbf{x}, t)$$";
            end

            xlabel("$x$ (greyordinates)", Interpreter="latex")
            ylabel("wrapped phase of $\phi$ (rad/$\pi$)", Interpreter="latex")
            if isemptytext(opts.title)
                title(atitle, Interpreter="latex")
            else
                title(opts.title, Interpreter="latex")
            end

            fontsize(scale=3)

            % save figure
            if ~isemptytext(opts.fig_fileprefix)
                saveFigure2(h, opts.fig_fileprefix, ext={'.png'}, closeFigure=opts.closeFigure);
            end
        end

        function h1 = rugplot_single_fft_zeta(this, opts)
            arguments
                this mlraut.AnalyticSignalHCPPar
                opts.fig_fileprefix {mustBeTextScalar} = ""
                opts.closeFigure logical = true
                opts.mat_fqfn {mustBeFile} = this.mat_fqfn
                opts.complex_representation function_handle = @real
            end
            if isemptytext(opts.fig_fileprefix)
                [pth, fp] = fileparts(opts.mat_fqfn);
                opts.fig_fileprefix = fullfile(pth, "fftzeta_as_" + fp);
            end

            ld = load(opts.mat_fqfn);
            psi = ld.this_subset.bold_signal;
            phi = ld.this_subset.physio_signal;
            zeta_t = this.zeta(psi, phi);
            zeta_nu = fft(zeta_t);

            h1 = this.rugplot(zeta_nu, ...
                is_spectral=true, fig_fileprefix=opts.fig_fileprefix, closeFigure=opts.closeFigure, ...                
                complex_representation=opts.complex_representation);
        end

        function [h1,h2] = rugplot_single_zeta(this, opts)
            arguments
                this mlraut.AnalyticSignalHCPPar
                opts.fig_fileprefix {mustBeTextScalar} = ""
                opts.closeFigure logical = true
                opts.mat_fqfn {mustBeFile} = this.mat_fqfn
                opts.complex_representation function_handle = @real
            end
            if isemptytext(opts.fig_fileprefix)
                [pth, fp] = fileparts(opts.mat_fqfn);
                opts.fig_fileprefix = fullfile(pth, "zeta_as_" + fp);
            end

            ld = load(opts.mat_fqfn);
            psi = ld.this_subset.bold_signal;
            phi = ld.this_subset.physio_signal;
            zeta_t = this.zeta(psi, phi);
            % zeta_t_avggo = mean(zeta_t, 2);

            h1 = this.rugplot(zeta_t, ...
                is_spectral=false, fig_fileprefix=opts.fig_fileprefix, closeFigure=opts.closeFigure, ...                
                complex_representation=opts.complex_representation);
            h2 = [];
            % h2 = figure;
            % plot(real(zeta_t_avggo)); title("$\langle \zeta(t) \rangle_{go}$", Interpreter="latex");
            % fontsize(scale=1.618)
            % if opts.closeFigure
            %     close(h2)
            % end
        end


        function this = AnalyticSignalHCPPar(varargin)
            this = this@mlraut.AnalyticSignalHCP(varargin{:});
        end
    end

    methods (Hidden)
        function [phi_nu,selected] = t_to_nu(this, phi_t, opts)
            arguments
                this mlraut.AnalyticSignalHCPPar
                phi_t {mustBeNumeric}
                opts.num_bins_angles {mustBeScalarOrEmpty} = []
            end
            assert(size(phi_t, 2) == 1)
            if isempty(opts.num_bins_angles) && ~isempty(this.twistors)
                opts.num_bins_angles = this.twistors.num_bins_angles;
            elseif isempty(opts.num_bins_angles)
                opts.num_bins_angles = 80;
            end
            num_bins_angles_ = opts.num_bins_angles;

            phi_nu = zeros(num_bins_angles_, 1);
            binlim = asrow(linspace(-pi, pi, num_bins_angles_ + 1));
            alpha = angle(phi_t);  % wrapped [-pi, pi]

            % average phi by phase bins
            selected = false(size(phi_t, 1), num_bins_angles_);
            for b = 2:num_bins_angles_+1
                selected(:, b-1) = binlim(b-1) < alpha & alpha < binlim(b);
                phi_nu(b-1,:) = mean(phi_t(selected(:, b-1), :), 1, "omitnan");
            end
        end

        function psi_t = nu_to_t(this, psi_nu, opts)
            arguments
                this mlraut.AnalyticSignalHCPPar
                psi_nu {mustBeNumeric}
                opts.phi_t {mustBeNumeric}
                opts.selected {mustBeNumericOrLogical}
                opts.num_bins_angles {mustBeScalarOrEmpty} = []
            end
            if isempty(opts.num_bins_angles) && ~isempty(this.twistors)
                opts.num_bins_angles = this.twistors.num_bins_angles;
            elseif isempty(opts.num_bins_angles)
                opts.num_bins_angles = 80;
            end
            num_bins_angles_ = opts.num_bins_angles;            

            % approx. reconstitute psi_t
            psi_t = zeros(size(opts.phi_t, 1), size(psi_nu, 2));
            for b = 1:num_bins_angles_
                Nsamples = sum(opts.selected(:, b), 1);
                psi_t(opts.selected(:, b), :) = repmat(psi_nu(b,:), [Nsamples, 1]);
            end
        end

        function h = rugplot_fft_zeta(this, opts)
            arguments
                this mlraut.AnalyticSignalHCPPar
                opts.fig_fileprefix {mustBeTextScalar} = ""
                opts.closeFigure logical = true
                opts.mat_fqfn {mustBeFile} = this.mat_fqfn
            end
            if isemptytext(opts.fig_fileprefix)
                [pth, fp] = fileparts(opts.mat_fqfn);
                opts.fig_fileprefix = fullfile(pth, fp);
            end

            ld = load(opts.mat_fqfn);
            zeta_nu = ld.zeta_nu_; % average, N=2,4, of zeta_nu

            nd = mlraut.NetworkData(this, zeta_nu);
            img = nd.reshape_by_anatomy(zeta_nu);
            img = fftshift(img, 1);

            % fig position
            fig_width = 2056;
            fig_height = 1329;
            h = figure('Position', [1, 1, fig_width, fig_height]);
            movegui(h, 'northwest');  % Moves to upper left

            imagesc(real(img));
            colorbar;
            clim([-500, 500])

            hold on

            % lines
            xline(8788, 'k-', 'LineWidth', 1);  % VIS ctx
            xline(20748, 'k-', 'LineWidth', 1);  % SMS
            xline(27510, 'k-', 'LineWidth', 1);  % DAN
            xline(34683, 'k-', 'LineWidth', 1);  % VAN
            xline(39219, 'k-', 'LineWidth', 1);  % LIM
            xline(46530, 'k-', 'LineWidth', 1);  % FPN
            xline(58666, 'k-', 'LineWidth', 1);  % DMN            
            xline(76498, 'k-', 'LineWidth', 1);  % cbm
            xline(79855, 'k-', 'LineWidth', 1);  % str
            % xline(82053, 'k-', 'LineWidth', 1);  % thal

            % annotation
            apos = 0.95*this.num_frames;
            text(8788-1200, apos, 'VIS', 'FontSize', 9, 'FontAngle', 'italic', 'Color', [.2 .2 .2], 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(20748-1200, apos, 'SMS', 'FontSize', 9, 'FontAngle', 'italic', 'Color', [.2 .2 .2], 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(27510-1200, apos, 'DAN', 'FontSize', 9, 'FontAngle', 'italic', 'Color', [.2 .2 .2], 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(34683-1200, apos, 'VAN', 'FontSize', 9, 'FontAngle', 'italic', 'Color', [.2 .2 .2], 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(39219-1200, apos, 'LIM', 'FontSize', 9, 'FontAngle', 'italic', 'Color', [.2 .2 .2], 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(46530-1200, apos, 'FPN', 'FontSize', 9, 'FontAngle', 'italic', 'Color', [.2 .2 .2], 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(58666-1200, apos, 'DMN', 'FontSize', 9, 'FontAngle', 'italic', 'Color', [.2 .2 .2], 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(76498-1200, apos, 'cerebellum', 'FontSize', 9, 'FontAngle', 'italic', 'Color', [.2 .2 .2], 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(79855-1200, apos, 'striatum', 'FontSize', 9, 'FontAngle', 'italic', 'Color', [.2 .2 .2], 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(82053-1200, apos, 'thalamus', 'FontSize', 9, 'FontAngle', 'italic', 'Color', [.2 .2 .2], 'HorizontalAlignment', 'left', 'Rotation', 90);

            hold off

            % Calculate time values
            Nomega = size(img, 1);
            Omega = (-Nomega/2:Nomega/2-1) * 2 * pi * this.Fs / Nomega; % total freq. range in radian/s
            % time_values = linspace(0, T, Nt);

            n_ticks = 11; % Choose how many tick marks you want
            tick_indices = linspace(1, Nomega, n_ticks);
            tick_omega = linspace(Omega(1), Omega(end), n_ticks);

            % indices -> times (s) 
            ax = gca;
            ax.YTick = tick_indices;
            % set(ax, ...
            %     'YTickLabel', arrayfun(@(x) sprintf("%.1f", x), tick_omega/pi, 'UniformOutput', false) + "\pi", ...
            %     'TickLabelInterpreter', 'latex');
             ax.YTickLabel = arrayfun(@(x) sprintf('%.1f', x), tick_omega, 'UniformOutput', false);

            % Add labels
            xlabel("$\mathbf{x}$ (greyordinates)", Interpreter="latex")
            ylabel("$\omega$ (rad/s)", Interpreter="latex")
            title("$$\Re \, \zeta(\mathbf{x}, \omega)$$", Interpreter="latex")

            fontsize(scale=3)

            % save figure
            saveFigure2(h, opts.fig_fileprefix, ext={'.fig', '.png'}, closeFigure=true);
        end

        function h = rugplot_neg_dbold_dt(this, bold_signal, opts)
            arguments
                this mlraut.AnalyticSignalHCPPar
                bold_signal {mustBeNumeric}
                opts.fig_fileprefix {mustBeTextScalar} = ""
            end

            nd = mlraut.NetworkData(this, bold_signal);
            img = nd.reshape_by_anatomy(this.neg_dbold_dt(bold_signal));

            % fig position
            fig_width = 2056;
            fig_height = 1329;
            h = figure('Position', [1, 1, fig_width, fig_height]);
            movegui(h, 'northwest');  % Moves to upper left

            imagesc(real(img ));
            colorbar;
            clim([-5, 5])

            hold on

            % lines
            xline(8788, 'k-', 'LineWidth', 1);  % VIS ctx
            xline(20748, 'k-', 'LineWidth', 1);  % SMS
            xline(27510, 'k-', 'LineWidth', 1);  % DAN
            xline(34683, 'k-', 'LineWidth', 1);  % VAN
            xline(39219, 'k-', 'LineWidth', 1);  % LIM
            xline(46530, 'k-', 'LineWidth', 1);  % FPN
            xline(58666, 'k-', 'LineWidth', 1);  % DMN            
            xline(76498, 'k-', 'LineWidth', 1);  % cbm
            xline(79855, 'k-', 'LineWidth', 1);  % str
            % xline(82053, 'k-', 'LineWidth', 1);  % thal

            % annotation
            apos = 0.95*this.num_frames;
            text(8788-1200, apos, 'VIS', 'FontSize', 9, 'FontAngle', 'italic', 'Color', [.2 .2 .2], 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(20748-1200, apos, 'SMS', 'FontSize', 9, 'FontAngle', 'italic', 'Color', [.2 .2 .2], 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(27510-1200, apos, 'DAN', 'FontSize', 9, 'FontAngle', 'italic', 'Color', [.2 .2 .2], 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(34683-1200, apos, 'VAN', 'FontSize', 9, 'FontAngle', 'italic', 'Color', [.2 .2 .2], 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(39219-1200, apos, 'LIM', 'FontSize', 9, 'FontAngle', 'italic', 'Color', [.2 .2 .2], 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(46530-1200, apos, 'FPN', 'FontSize', 9, 'FontAngle', 'italic', 'Color', [.2 .2 .2], 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(58666-1200, apos, 'DMN', 'FontSize', 9, 'FontAngle', 'italic', 'Color', [.2 .2 .2], 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(76498-1200, apos, 'cerebellum', 'FontSize', 9, 'FontAngle', 'italic', 'Color', [.2 .2 .2], 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(79855-1200, apos, 'striatum', 'FontSize', 9, 'FontAngle', 'italic', 'Color', [.2 .2 .2], 'HorizontalAlignment', 'left', 'Rotation', 90);
            text(82053-1200, apos, 'thalamus', 'FontSize', 9, 'FontAngle', 'italic', 'Color', [.2 .2 .2], 'HorizontalAlignment', 'left', 'Rotation', 90);

            hold off

            % Calculate time values
            Nt = size(img, 1);
            T = Nt * this.tr; % total time in seconds
            % time_values = linspace(0, T, Nt);

            n_ticks = 11; % Choose how many tick marks you want
            tick_indices = linspace(1, Nt, n_ticks);
            tick_times = linspace(0, T, n_ticks);

            % indices -> times (s) 
            ax = gca;
            ax.YTick = tick_indices;
            ax.YTickLabel = arrayfun(@(x) sprintf('%.0f', x), tick_times, 'UniformOutput', false);

            % Add labels
            xlabel("$\mathbf{x}$ (greyordinates)", Interpreter="latex")
            ylabel("time (s)")
            title("$$-\frac{d}{dt}BOLD(\mathbf{x}, t)$$", Interpreter="latex")

            fontsize(scale=3)

            % save figure
            if ~isemptytext(opts.fig_fileprefix)
                saveFigure2(h, opts.fig_fileprefix, ext={'.fig', '.png'}, closeFigure=true);
            end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
