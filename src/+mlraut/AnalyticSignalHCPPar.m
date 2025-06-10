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
                source_physio="iFV-brightest", ...
                tags="ASHCPPar-mean-twistor");
            call(this);  % malloc & construct delegates
            this.out_dir = out_dir;  % for group averages

            % init
            mats = asrow(glob(fullfile(out_dir, '*/sub-*_ses-*ASHCPPar*.mat')));
            if nmats(2) > length(mats)
                nmats(2) = length(mats);
            end
            mats = mats(nmats(1):nmats(2));  % 3690 mat files
            nbin = this.cifti.Nbins;
            ngo = this.num_nodes;
            X_ = zeros(nbin, ngo);
            reY_ = zeros(nbin, ngo);
            imY_ = zeros(nbin, ngo);
            Z_ = zeros(nbin, ngo);
            T_ = zeros(nbin, ngo);
            r_ = zeros(1, ngo);
            bold_ = zeros(nbin, ngo);
            plvs_ = zeros(nbin, ngo);
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
                    plvs_ = plvs_ + real(plvs);
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
            plvs_ = plvs_/nmats_corr;


            dt = 2*pi/this.cifti_.Nbins;
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
                plvs_, sprintf('plvs_as_sub-%s_ses-%s_%s', ntag, ntag, this.tags), dt=dt, units_t=units_t);

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
                        "iFV-brightest", "iFV", "iFV-quantile", "sFV", "3rdV", "latV", "centrumsemiovale", "ctx", "RV", "HRV"], ...
                        tags="ASHCPPar-server-call");
                    call(as);
                catch ME
                    handwarning(ME)
                end
            end
        end
    
        %% running call on cluster

        function [j,c,msg,id] = cluster_construct_and_call(globbing_mat, opts)
            %% for clusters running Matlab parallel server

            arguments
                globbing_mat {mustBeFile} = ...
                    fullfile( ...
                    getenv('SINGULARITY_HOME'), 'AnalyticSignalHCP', 'mlraut_AnalyticSignalHCPPar_globbing.mat')
                opts.sub_indices double = []  % total ~ 1:1113
                opts.globbing_var = "globbed"
            end
            ld = load(globbing_mat);
            globbed = convertStringsToChars(ld.(opts.globbing_var));
            globbed = asrow(globbed);
            if ~isempty(opts.sub_indices)
                globbed = globbed(opts.sub_indices);
            end

            % pad and reshape globbed
            Ncol = 8;
            Nrow = ceil(numel(globbed)/Ncol);
            padding = repmat("", [1, Ncol*Nrow - numel(globbed)]);
            globbed = [globbed, padding];
            globbed = reshape(globbed, Nrow, Ncol);
            globbed = convertStringsToChars(globbed);
            fprintf("%s:globbed:\n", stackstr())
            disp(size(globbed))
            disp(ascol(globbed))

            % contact cluster slurm

            warning('off', 'MATLAB:legacy:batchSyntax');
            warning('off', 'parallel:convenience:BatchFunctionNestedCellArray');
            warning('off', 'MATLAB:TooManyInputs');

            c = mlraut.CHPC3.propcluster('joshua_shimony', mempercpu='40gb', walltime='00:30:00');
            disp(c.AdditionalProperties)
            for irow = 1:Nrow
                try
                    j = c.batch( ...
                        @mlraut.AnalyticSignalHCPPar.construct_and_call, ...
                        1, ...
                        {globbed(irow, :)}, ...
                        'Pool', Ncol, ...
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
                globbing_mat {mustBeFile} = ...
                    fullfile( ...
                    getenv('SINGULARITY_HOME'), 'AnalyticSignalHCP', 'mlraut_AnalyticSignalHCPPar_globbing.mat')
                opts.globbing_var = "globbed"
                opts.sub_range = []  % total ~ 1:1113
                opts.new_physio {mustBeText} = "iFV-brightest"
                opts.test_range = []  % 1:2
                opts.Ncol {mustBeInteger} = 32
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

            c = mlraut.CHPC3.propcluster('joshua_shimony', mempercpu='40gb', walltime='24:00:00');
            disp(c.AdditionalProperties)

            for col_idx = 1:opts.Ncol
                try
                    j = c.batch( ...
                        @mlraut.AnalyticSignalHCPPar.construct_means, ...
                        1, ...
                        {globbed(:, col_idx), 'col_idx', col_idx, 'new_physio', opts.new_physio, 'test_range', opts.test_range}, ...
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
                opts.out_dir {mustBeFolder} = "/scratch/jjlee/Singularity/AnalyticSignalHCP"
            end
            
            subjects = subjects(~cellfun(@isempty, subjects));  % Remove empty cells
            durations = nan(1, length(subjects));

            parfor sidx = 1:length(subjects)

                tic;
            
                % setup
                mlraut.CHPC3.setenvs();
                ensuredir(opts.out_dir); %#ok<*PFBNS>
                ensuredir(fullfile(opts.out_dir, subjects{sidx}));

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
                        "iFV-brightest", "iFV", "iFV-quantile", "sFV", "3rdV", "latV", "centrumsemiovale", "ctx", "RV", "HRV"], ...
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
                opts.col_idx {mustBeInteger}
                opts.new_physio {mustBeText} = ""
                opts.test_range = []
            end

            tic;

            % setup
            mlraut.CHPC3.setenvs();
            ensuredir(opts.out_dir); %#ok<*PFBNS>

            mat_fqfn = [];
            idx = 0;
            while isempty(mat_fqfn) && idx < length(subjects)
                idx = idx + 1;
                mat_fqfn = mglob(fullfile( ...
                    opts.out_dir, subjects{idx}, sprintf("sub-%s_ses-*_proc-*ASHCPPar*.mat", subjects{idx})));
            end
            as = mlraut.AnalyticSignalHCPPar.load(mat_fqfn(1), class="mlraut.AnalyticSignalHCPPar");
            as.out_dir = opts.out_dir;
            as.mean_twistor_instance( ...
                subjects, col_idx=opts.col_idx, new_physio=opts.new_physio, test_range=opts.test_range);

            durations = toc;
        end
    
        function gather_means(out_dir, opts)
            arguments
                out_dir {mustBeFolder} = pwd
                opts.measures {mustBeText} = [ ...
                    "comparator", "bold", "neg-dbold-dt", "plvs", "X", "reY", "imY", "Z", "T"]
                opts.physio {mustBeText} = "iFV-brightest"
                opts.gsr logical = true
                opts.ddt logical = true
                opts.tag {mustBeTextScalar} = "ASHCPPar"
                opts.Nphase_out {mustBeInteger} = 8
            end

            for measure = opts.measures
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
                    g_last = g;
                    cii = cifti_read(g);
                    mats = [mats, cii.cdata'];  %#ok<*AGROW> % {Nt x Ngo} x Nmats; Nt x Ngo is preferred for internal repr.
                    re = regexp(g, "\S_sub-n(?<subn>\d+)_ses-n(?<sesn>\d+)_\S", "names");
                    subn = str2double(re.subn);
                    sesn = str2double(re.sesn);
                    assert(subn == sesn)
                    n = n + subn;
                    weights = [weights, subn];
                end
                weights = weights / n;

                % apply weights to mats
                mats_weighted = cellfun(@(mat, weight) mat * weight, mats, num2cell(weights), 'UniformOutput', false);
                
                % concatenate arrays along the 3rd dimension, then sum
                mat_avg = sum(cat(3, mats_weighted{:}), 3);

                % filename
                fqfn = strrep(g_last, re.subn, string(n));
                fqfn = regexprep(fqfn, "_par\d+", "");

                % down-sample phase:  40 x Ngo -> 8 x Ngo
                Nphase = size(mat_avg, 1);
                if Nphase > 1 && mod(Nphase, opts.Nphase_out) == 0
                    if contains(measure, "bold") || contains(measure, "plv")
                        Navg = Nphase / opts.Nphase_out;
                        mat_avg = squeeze(mean(reshape(mat_avg, Navg, opts.Nphase_out, []), 1));
                        cii.diminfo{2} = cifti_diminfo_make_series( ...
                            opts.Nphase_out, 0, 2*pi/opts.Nphase_out, 'RADIAN');
                    else
                        mat_avg = mean(mat_avg, 1);
                        cii.diminfo{2} = cifti_diminfo_make_scalars(1);
                        fqfn = strrep(fqfn, ".dtseries.nii", ".dscalar.nii");
                    end
                end

                % save weighted average
                cii.cdata = mat_avg';
                cifti_write(cii, convertStringsToChars(fqfn));
            end
        end
    end

    methods
        function this = mean_twistor_instance(this, subs, opts)
            %% accepts text array subs; find mat files for subs; writes files indexed by opts.col_idx

            arguments
                this mlraut.AnalyticSignalHCPPar
                subs {mustBeText}
                opts.col_idx {mustBeInteger}
                opts.new_physio {mustBeText} = ""
                opts.test_range = []
            end
            subs = convertCharsToStrings(subs);

            % DEBUG
            disp(subs)

            warning("off", "MATLAB:class:LoadDefinitionUpdated");

            % \emph{this} supplies utilities

            % init
            nbin = this.cifti.Nbins;
            ngo = this.num_nodes;
            X_ = zeros(nbin, ngo);
            reY_ = zeros(nbin, ngo);
            imY_ = zeros(nbin, ngo);
            Z_ = zeros(nbin, ngo);
            T_ = zeros(nbin, ngo);
            r_ = zeros(1, ngo);
            bold_ = zeros(nbin, ngo);
            neg_dbold_dt_ = zeros(nbin, ngo);
            plvs_ = zeros(nbin, ngo);
            nmats_corr = 0;

            % abbrev.
            bin = @this.bin_by_physio_angle;

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
                        elseif strcmpi(opts.new_physio, 'vis')
                            phi = ld.this_subset.HCP_signals.ctx.psi(:, 1);
                        elseif strcmpi(opts.new_physio, 'sms')
                            phi = ld.this_subset.HCP_signals.ctx.psi(:, 2);
                        elseif strcmpi(opts.new_physio, 'dan')
                            phi = ld.this_subset.HCP_signals.ctx.psi(:, 3);
                        elseif strcmpi(opts.new_physio, 'van')
                            phi = ld.this_subset.HCP_signals.ctx.psi(:, 4);
                        elseif strcmpi(opts.new_physio, 'lim')
                            phi = ld.this_subset.HCP_signals.ctx.psi(:, 5);
                        elseif strcmpi(opts.new_physio, 'fpn')
                            phi = ld.this_subset.HCP_signals.ctx.psi(:, 6);
                        elseif strcmpi(opts.new_physio, 'dmn')
                            phi = ld.this_subset.HCP_signals.ctx.psi(:, 7);
                        else
                            re_phi = ld.this_subset.physio_supplementary(opts.new_physio);
                            assert(~isempty(re_phi))
                            assert(all(isfinite(re_phi)))
                            phi = hilbert(re_phi);
                        end
                        X = bin(this.X(psi, phi), phi);
                        Y = bin(this.Y(psi, phi), phi);
                        Z = bin(this.Z(psi, phi), phi);
                        T = bin(this.T(psi, phi), phi);
                        r = asrow(this.connectivity(psi, phi));
                        bold = bin(psi, phi);
                        neg_dbold_dt = bin(this.neg_dbold_dt(psi, phi), phi);
                        plvs = bin(this.phase_locked_values(psi, phi), phi);

                        X_ = X_ + real(X);
                        reY_ = reY_ + real(Y);
                        imY_ = imY_ + imag(Y);
                        Z_ = Z_ + real(Z);
                        T_ = T_ + real(T);
                        r_ = r_ + real(r);
                        bold_ = bold_ + real(bold);
                        neg_dbold_dt_ = neg_dbold_dt_ + real(neg_dbold_dt);
                        plvs_ = plvs_ + real(plvs);
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
            X_ = X_/nmats_corr;
            reY_ = reY_/nmats_corr;
            imY_ = imY_/nmats_corr;
            Z_ = Z_/nmats_corr;
            T_ = T_/nmats_corr;
            r_ = r_/nmats_corr;
            bold_ = bold_/nmats_corr;
            neg_dbold_dt_ = neg_dbold_dt_/nmats_corr;
            plvs_ = plvs_/nmats_corr;

            dt = 2*pi/this.cifti_.Nbins;
            units_t = "RADIAN";
            ntag = sprintf('n%g', nmats_corr);
            if ~isemptytext(opts.new_physio)
                ptags = strrep(this.tags, this.source_physio, opts.new_physio);
            end
            this.cifti.write_cifti( ...
                X_, sprintf('X_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                reY_, sprintf('reY_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                imY_, sprintf('imY_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                Z_, sprintf('Z_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                T_, sprintf('T_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                r_, sprintf('comparator_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                bold_, sprintf('bold_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                neg_dbold_dt_, sprintf('neg-dbold-dt_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                plvs_, sprintf('plvs_as_sub-%s_ses-%s_%s_par%i', ntag, ntag, ptags, opts.col_idx), dt=dt, units_t=units_t);
            
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

        function this = AnalyticSignalHCPPar(varargin)
            this = this@mlraut.AnalyticSignalHCP(varargin{:});
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
