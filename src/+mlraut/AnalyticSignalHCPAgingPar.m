classdef AnalyticSignalHCPAgingPar < handle & mlraut.AnalyticSignalHCPAging
    %% nlim = 10, 689
    %  
    %  Created 07-Feb-2024 23:37:28 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    
    methods (Static)

        %% running {mean,var}_* on cluster

        function [j,c] = cluster_batch_stats(opts)
            %% for clusters running Matlab parallel server
            %  Args:
            %      opts.nlim double = []
            %      opts.twists {mustBeText} = ["connectivity", "X", "Y", "Z", "T", "angle", "unwrap"]
            %      opts.funhs function_handle = [@mlraut.AnalyticSignalHCPAgingPar.mean_twistor, ...
            %                                    @mlraut.AnalyticSignalHCPAgingPar.construct_vars]

            arguments
                opts.nlim double = []
                opts.twists {mustBeText} = ["connectivity", "X", "Y", "T", "Z", "angle", "unwrap"]
                opts.funhs function_handle = [@mlraut.AnalyticSignalHCPAgingPar.construct_means, ...
                                              @mlraut.AnalyticSignalHCPAgingPar.construct_vars]
            end
            c = mlraut.CHPC3.propcluster('joshua_shimony', mempercpu=32, walltime='120:00:00');

            for f = opts.funhs
                try
                    j = c.batch( ...
                        f, ...
                        1, ...
                        {'nlim', opts.nlim, 'twists', opts.twists}, ...
                        'CurrentFolder', '/scratch/jjlee/Singularity/AnalyticSignalHCPAging', ...
                        'AutoAddClientPath', false);
                catch ME
                    handwarning(ME)
                end
            end
        end

        function [j,c] = cluster_batch_stats_rsn(opts)
            %% for clusters running Matlab parallel server
            %  Args:
            %      opts.nlim double = []
            %      opts.twists {mustBeText} = ["connectivity", "X", "Y", "Z", "T", "angle", "unwrap"]
            %      opts.rsns {mustBeNumeric} = 1:7
            %      opts.funhs function_handle = [@mlraut.AnalyticSignalHCPAgingPar.mean_twistor_rsn, ...
            %                                    @mlraut.AnalyticSignalHCPAgingPar.var_twistor_rsn]

            arguments
                opts.nlim double = []
                opts.twists {mustBeText} = ["connectivity", "X", "Y", "T", "Z", "angle", "unwrap"]
                opts.rsns {mustBeNumeric} = 1:7
                opts.funhs function_handle = [@mlraut.AnalyticSignalHCPAgingPar.mean_twistor_rsn, ...
                                              @mlraut.AnalyticSignalHCPAgingPar.var_twistor_rsn]
            end
            c = mlraut.CHPC3.propcluster('aristeidis_sotiras', mempercpu=32, walltime='120:00:00');

            for f = opts.funhs
                try
                    j = c.batch( ...
                        f, ...
                        1, ...
                        {'nlim', opts.nlim, 'twists', opts.twists, 'rsns', opts.rsns}, ...
                        'CurrentFolder', '/scratch/jjlee/Singularity/AnalyticSignalHCPAging', ...
                        'AutoAddClientPath', false);
                catch ME
                    handwarning(ME)
                end
            end
        end

        function ret = mean_comparator(opts)
            arguments
                opts.nlim = [];
            end
            nlim = opts.nlim;

            %%

            ret = 0;
            mlraut.CHPC3.setenvs();
            this = mlraut.AnalyticSignalHCPAgingPar( ...
                subjects={'HCA9992517_V1_MR'}, ...
                tasks={'fMRI_CONCAT_ALL'}, ...
                do_7T=false, ...
                do_resting=true, ...
                do_task=false, ...
                do_save=false, ...
                do_save_dynamic=false, ...
                do_save_ciftis=false, ...
                plot_range=1:225, ...
                tags="mean-twistor");

            mats = flip(asrow(mglob(fullfile(this.out_dir, "HCA*_MR/sub-*_ses-*-iFV-brightest-scaleiqr-AnalyticSignalHCPAgingPar.mat"))));
            if ~isempty(nlim)
                mats = mats(1:nlim);
            end
            nsub = length(mats);
            nx = this.num_nodes;
            nfail = 0;
            comparator_ = zeros(1, nx);

            fprintf(stackstr() + "\n");
            for mat = mats
                try
                    fprintf("loading %s\n", mat);
                    tic

                    ld = load(mat);
                    comparator_ = comparator_ + asrow(ld.this.comparator)/nsub;
                    ret = ret + 1;
                    
                    toc
                catch ME
                    handwarning(ME)
                    nfail = nfail + 1;
                end
            end
            
            if nfail > 0
                assert(nfail < nsub)
                comparator_ = comparator_*(nsub/(nsub - nfail));
            end

            %% write summary averages

            this.out_dir = fullfile(getenv('SINGULARITY_HOME'), 'AnalyticSignalHCPAging');
            if ~isempty(nlim)
                tags = this.tags(sprintf("nlim%i", nlim));
            end

            this.write_ciftis( ...
                comparator_, ...
                sprintf('mean_comparator_as_sub-all_ses-all_%s', tags));
        end

        function durations = construct_means(opts)
            arguments
                opts.nlim = []
                opts.twists {mustBeText} = ["connectivity", "X", "Y", "T", "Z", "a", "u", "bold", "plvs"]
                opts.tags {mustBeTextScalar} = "construct-means"
                opts.out_dir {mustBeFolder} = "/scratch/jjlee/Singularity/AnalyticSignalHCPAging"
                opts.mat_patt {mustBeTextScalar} = "HCA*_MR/sub-*_ses-*_proc-iFV-brightest-*-subset-ASHCPPar.mat"
            end
            durations = nan(1, length(opts.twists));

            %%

            parfor tidx = 1:length(opts.twists)
                tic;
            
                % setup
                mlraut.CHPC3.setenvs();
                ensuredir(opts.out_dir); %#ok<*PFBNS>
                twist = opts.twists(tidx);

                % load AnalyticSignalHCPAgingPar.this_subset;
                % construct AnalyticSignalHCPAgingPar for its method functions
                mat_files = mglob(fullfile(opts.out_dir, opts.mat_patt));
                if ~isempty(opts.nlim)
                    mat_files = mat_files(1:opts.nlim);
                end
                as1 = mlraut.AnalyticSignalHCPAgingPar.load(mat_files(1));

                %% accumulate stats

                % init
                N_mat = length(mat_files);
                N_go = as1.num_nodes;
                N_fail = 0;
                E_ = zeros(1, N_go);
                E_suppl_ = as1.physio_supplementary;  % containers.Map is handle
                for kidx = 1:length(E_suppl_.keys)
                    E_suppl_(E_suppl_.keys{kidx}) = zeros(1, N_go);
                end

                % accumulate granules
                for imat = 2:length(mat_files)
                    try
                        % fprintf("loading %s\n", mat_file);
                        as = mlraut.AnalyticSignalHCPAgingPar.load(mat_files(imat));

                        E_norm = as.(twist)(as.this.bold_signal, as.this.physio_signal);
                        E_ = E_ + E_norm/N_mat;




                        % physio_suppl = copy(as.this.physio_supplementary);



                        

                    catch ME
                        handwarning(ME)
                        N_fail = N_fail + 1;
                    end
                end
                if N_fail > 0
                    assert(N_fail < N_mat)
                    E_ = E_*(N_mat/(N_mat - N_fail));
                end

                % write stats
                as.out_dir = fullfile(getenv('SINGULARITY_HOME'), 'AnalyticSignalHCPAging');
                if ~isempty(opts.nlim)
                    tags = as.tags(sprintf("nlim%i", opts.nlim));
                else
                    tags = as.tags();
                end
                as.write_cifti( ...
                    E_, ...
                    sprintf("mean_%s_as_sub-all_ses-all_proc-%s", twist, tags));

                durations(tidx) = toc;
            end
        end

        function ret = mean_twistor_rsn(opts)
            arguments
                opts.nlim = 7
                opts.rsn = 7
            end
            nlim = opts.nlim;
            rsn = opts.rsn;

            %%

            ret =  0;
            mlraut.CHPC3.setenvs();
            this = mlraut.AnalyticSignalHCPAgingPar( ...
                subjects={'HCA9992517_V1_MR'}, ...
                tasks={'fMRI_CONCAT_ALL'}, ...
                do_7T=false, ...
                do_resting=true, ...
                do_task=false, ...
                do_save=false, ...
                do_save_dynamic=false, ...
                do_save_ciftis=false, ...
                plot_range=1:225, ...
                tags="mean-twistor");

            mats = flip(asrow(mglob(fullfile(this.out_dir, "HCA*_MR/sub-*_ses-*-iFV-brightest-scaleiqr-AnalyticSignalHCPAgingPar.mat"))));
            if ~isempty(nlim)
                mats = mats(1:nlim);
            end
            nsub = length(mats);
            nx = this.num_nodes;
            nfail = 0;

            T_ = zeros(1, nx);
            X_ = zeros(1, nx);
            Y_ = zeros(1, nx);
            Z_ = zeros(1, nx);
            angle_ = zeros(1, nx);
            unwrap_ = zeros(1, nx);

            fprintf(stackstr() + "\n");
            for mat = mats
                try
                    fprintf("loading %s\n", mat);
                    tic

                    ld = load(mat);

                    ctx = ld.this.HCP_signals.ctx;
                    angle_rsn = this.angle(ctx.psi(:,rsn), ctx.phi(:,rsn));
                    t_interesting = cos(angle_rsn) > 0;

                    T_rsn = this.sample_rsn( ...
                        this.T(ld.this.bold_signal, ld.this.physio_signal), t_interesting);
                    Z_rsn = this.sample_rsn( ...
                        this.Z(ld.this.bold_signal, ld.this.physio_signal), t_interesting);
                    X_rsn = this.sample_rsn( ...
                        this.X(ld.this.bold_signal, ld.this.physio_signal), t_interesting);
                    Y_rsn = this.sample_rsn( ...
                        this.Y(ld.this.bold_signal, ld.this.physio_signal), t_interesting);
                    angle_rsn = this.sample_rsn( ...
                        this.angle(ld.this.bold_signal, ld.this.physio_signal), t_interesting);
                    unwrap_rsn = this.sample_rsn( ...
                        this.unwrap(ld.this.bold_signal, ld.this.physio_signal), t_interesting);

                    T_ = T_ + T_rsn/nsub;
                    X_ = X_ + X_rsn/nsub;
                    Y_ = Y_ + Y_rsn/nsub;
                    Z_ = Z_ + Z_rsn/nsub;
                    angle_ = angle_ + angle_rsn/nsub;
                    unwrap_ = unwrap_ + unwrap_rsn/nsub;
                    ret = ret + 1;

                    toc
                catch ME
                    handwarning(ME)
                    nfail = nfail + 1;
                end
            end
            
            if nfail > 0
                assert(nfail < nsub)
                T_ = T_*(nsub/(nsub - nfail));
                X_ = X_*(nsub/(nsub - nfail));
                Y_ = Y_*(nsub/(nsub - nfail));
                Z_ = Z_*(nsub/(nsub - nfail));
                angle_ = angle_*(nsub/(nsub - nfail));
                unwrap_ = unwrap_*(nsub/(nsub - nfail));
            end

            %% write summary averages

            this.out_dir = fullfile(getenv('SINGULARITY_HOME'), 'AnalyticSignalHCPAging');
            if ~isempty(nlim)
                tags = this.tags(sprintf("rsn%i-nlim%i", rsn, nlim));
            else
                tags = this.tags(sprintf("rsn%i", rsn));
            end

            this.write_ciftis( ...
                T_, ...
                sprintf('mean_T_as_sub-all_ses-all_%s', tags), ...
                partitions=[], ...
                do_save_dynamic=false);
            this.write_ciftis( ...
                X_, ...
                sprintf('mean_X_as_sub-all_ses-all_%s', tags), ...
                partitions=[], ...
                do_save_dynamic=false);
            this.write_ciftis( ...
                Y_, ...
                sprintf('mean_Y_as_sub-all_ses-all_%s', tags), ...
                partitions=[], ...
                do_save_dynamic=false);
            this.write_ciftis( ...
                Z_, ...
                sprintf('mean_Z_as_sub-all_ses-all_%s', tags), ...
                partitions=[], ...
                do_save_dynamic=false);
            this.write_ciftis( ...
                angle_, ...
                sprintf('mean_angle_as_sub-all_ses-all_%s', tags), ...
                partitions=[], ...
                do_save_dynamic=false);
            this.write_ciftis( ...
                unwrap_, ...
                sprintf('mean_unwrap_as_sub-all_ses-all_%s', tags), ...
                partitions=[], ...
                do_save_dynamic=false);
        end

        function ret = var_comparator(opts)
            arguments
                opts.nlim = [];
            end
            nlim = opts.nlim;

            %%

            ret = 0;
            mlraut.CHPC3.setenvs();
            this = mlraut.AnalyticSignalHCPAgingPar( ...
                subjects={'HCA9992517_V1_MR'}, ...
                tasks={'fMRI_CONCAT_ALL'}, ...
                do_7T=false, ...
                do_resting=true, ...
                do_task=false, ...
                do_save=false, ...
                do_save_dynamic=false, ...
                do_save_ciftis=false, ...
                plot_range=1:225, ...
                tags="mean-twistor");

            if ~isempty(nlim)
                tags_ = this.tags(sprintf('nlim%i', nlim));
            end
            %tags_ = this.tags(sprintf('nlim%i', 689));  % DEBUGGING
            mu = cifti_read( ...
                fullfile(this.out_dir, sprintf('mean_comparator_as_sub-all_ses-all_%s_avgt.dscalar.nii', tags_)));

            mats = flip(asrow(mglob(fullfile(this.out_dir, "HCA*_MR/sub-*_ses-*-iFV--brightest-scaleiqr-AnalyticSignalHCPAgingPar.mat"))));
            if ~isempty(nlim)
                mats = mats(1:nlim);
            end
            nsub = length(mats);
            nx = this.num_nodes;
            nfail = 0;
            comparator_ = zeros(1, nx);

            fprintf(stackstr() + "\n");
            for mat = mats
                try
                    fprintf("loading %s\n", mat);
                    tic

                    ld = load(mat);
                    comparator_ = comparator_ + ...
                        (asrow(ld.this.comparator) - asrow(mu.cdata)).^2/nsub;
                    ret = ret + 1;
                    
                    toc
                catch ME
                    handwarning(ME)
                    nfail = nfail + 1;
                end
            end
            
            if nfail > 0
                assert(nfail < nsub)
                comparator_ = comparator_*(nsub/(nsub - nfail));
            end

            %% write summary averages

            this.out_dir = fullfile(getenv('SINGULARITY_HOME'), 'AnalyticSignalHCPAging');
            if ~isempty(nlim)
                tags = this.tags(sprintf("nlim%i", nlim));
            end

            this.write_ciftis( ...
                comparator_, ...
                sprintf('var_comparator_as_sub-all_ses-all_%s', tags));
        end

        function ret = construct_vars(opts)
            arguments
                opts.nlim = []
                opts.twist {mustBeTextScalar} = "X"
            end
            nlim = opts.nlim;
            twist = opts.twist;

            assert(any(contains(["T", "X", "Y", "Z", "angle", "unwrap"], twist)))

            %%

            ret =  0;
            mlraut.CHPC3.setenvs();
            this = mlraut.AnalyticSignalHCPAgingPar( ...
                subjects={'HCA9992517_V1_MR'}, ...
                tasks={'fMRI_CONCAT_ALL'}, ...
                do_7T=false, ...
                do_resting=true, ...
                do_task=false, ...
                do_save=false, ...
                do_save_dynamic=false, ...
                do_save_ciftis=false, ...
                plot_range=1:225, ...
                tags="mean-element");

            if ~isempty(nlim)
                tags_ = this.tags(sprintf('nlim%i', nlim));
            end
            mu = cifti_read( ...
                fullfile(char(this.out_dir), ...
                sprintf('mean_%s_as_sub-all_ses-all_%s_avgt.dscalar.nii', char(twist), char(tags_))));

            mats = flip(asrow(mglob(fullfile(this.out_dir, "HCA*_MR/sub-*_ses-*-iFV--brightest-scaleiqr-AnalyticSignalHCPAgingPar.mat"))));
            if ~isempty(nlim)
                mats = mats(1:nlim);
            end
            nsub = length(mats);
            nx = this.num_nodes;
            nfail = 0;

            E_ = zeros(1, nx);

            fprintf(stackstr() + "\n");
            for mat = mats
                try
                    fprintf("loading %s\n", mat);
                    tic

                    ld = load(mat);
                    E_norm = this.sample_rsn( ...
                        this.(twist)(ld.this.bold_signal, ld.this.physio_signal));
                    E_ = E_ + abs(asrow(E_norm) - asrow(mu.cdata)).^2/nsub;

                    ret = ret + 1;

                    toc
                catch ME
                    handwarning(ME)
                    nfail = nfail + 1;
                end
            end
            
            if nfail > 0
                assert(nfail < nsub)
                E_ = E_*(nsub/(nsub - nfail));
            end

            %% write summary averages

            this.out_dir = fullfile(getenv('SINGULARITY_HOME'), 'AnalyticSignalHCPAging');
            if ~isempty(nlim)
                tags = this.tags(sprintf("nlim%i", nlim));
            end

            this.write_ciftis( ...
                E_, ...
                sprintf("var_%s_as_sub-all_ses-all_%s", twist, tags), ...
                partitions=[], ...
                do_save_dynamic=false);
        end

        function ret = var_twistor_rsn(opts)
            arguments
                opts.nlim = []
                opts.rsn = 7
            end
            nlim = opts.nlim;
            rsn = opts.rsn;

            %%

            ret = 0;
            mlraut.CHPC3.setenvs();
            this = mlraut.AnalyticSignalHCPAgingPar( ...
                subjects={'HCA9992517_V1_MR'}, ...
                tasks={'fMRI_CONCAT_ALL'}, ...
                do_7T=false, ...
                do_resting=true, ...
                do_task=false, ...
                do_save=false, ...
                do_save_dynamic=false, ...
                do_save_ciftis=false, ...
                plot_range=1:225, ...
                tags="mean-twistor");

            if ~isempty(nlim)
                tags_ = this.tags(sprintf('rsn%i-nlim%i', rsn, nlim));
            else
                tags_ = this.tags(sprinf('rsn%i', rsn));
            end
            %tags_ = this.tags(sprintf('rsn%i-nlim%i', rsn, 689));  % DEBUGGING
            mu_T = cifti_read( ...
                fullfile(this.out_dir, ...
                sprintf('mean_T_as_sub-all_ses-all_%s_avgt.dscalar.nii', tags_)));
            mu_X = cifti_read( ...
                fullfile(this.out_dir, ...
                sprintf('mean_X_as_sub-all_ses-all_%s_avgt.dscalar.nii', tags_)));
            mu_Y = cifti_read( ...
                fullfile(this.out_dir, ...
                sprintf('mean_Y_as_sub-all_ses-all_%s_avgt.dscalar.nii', tags_)));
            mu_Z = cifti_read( ...
                fullfile(this.out_dir, ...
                sprintf('mean_Z_as_sub-all_ses-all_%s_avgt.dscalar.nii', tags_)));
            mu_angle = cifti_read( ...
                fullfile(this.out_dir, ...
                sprintf('mean_angle_as_sub-all_ses-all_%s_avgt.dscalar.nii', tags_)));
            mu_unwrap = cifti_read( ...
                fullfile(this.out_dir, ...
                sprintf('mean_unwrap_as_sub-all_ses-all_%s_avgt.dscalar.nii', tags_)));

            mats = flip(asrow(mglob(fullfile(this.out_dir, "HCA*_MR/sub-*_ses-*-iFV--brightest-scaleiqr-AnalyticSignalHCPAgingPar.mat"))));
            if ~isempty(nlim)
                mats = mats(1:nlim);
            end
            nsub = length(mats);
            nx = this.num_nodes;
            nfail = 0;

            T_ = zeros(1, nx);
            X_ = zeros(1, nx);
            Y_ = zeros(1, nx);
            Z_ = zeros(1, nx);
            angle_ = zeros(1, nx);
            unwrap_ = zeros(1, nx);

            fprintf(stackstr() + "\n");
            for mat = mats
                try
                    fprintf("loading %s\n", mat);
                    tic

                    ld = load(mat);
                    ctx = ld.this.HCP_signals.ctx;

                    angle_rsn = this.angle(ctx.psi(:,rsn), ctx.phi(:,rsn));
                    t_interesting = cos(angle_rsn) > 0;

                    T_rsn = this.sample_rsn( ...
                        this.T(ld.this.bold_signal, ld.this.physio_signal), t_interesting);
                    X_rsn = this.sample_rsn( ...
                        this.X(ld.this.bold_signal, ld.this.physio_signal), t_interesting);
                    Y_rsn = this.sample_rsn( ...
                        this.Y(ld.this.bold_signal, ld.this.physio_signal), t_interesting);
                    Z_rsn = this.sample_rsn( ...
                        this.Z(ld.this.bold_signal, ld.this.physio_signal), t_interesting);
                    angle_rsn = this.sample_rsn( ...
                        this.angle(ld.this.bold_signal, ld.this.physio_signal), t_interesting);
                    unwrap_rsn = this.sample_rsn( ...
                        this.unwrap(ld.this.bold_signal, ld.this.physio_signal), t_interesting);

                    T_ = T_ + abs(T_rsn - asrow(mu_T.cdata)).^2/nsub;
                    X_ = X_ + abs(X_rsn - asrow(mu_X.cdata)).^2/nsub;
                    Y_ = Y_ + abs(Y_rsn - asrow(mu_Y.cdata)).^2/nsub;
                    Z_ = Z_ + abs(Z_rsn - asrow(mu_Z.cdata)).^2/nsub;
                    angle_ = angle_ + (angle_rsn - asrow(mu_angle.cdata)).^2/nsub;
                    unwrap_ = unwrap_ + (unwrap_rsn - asrow(mu_unwrap.cdata)).^2/nsub;
                    ret = ret + 1;

                    toc
                catch ME
                    handwarning(ME)
                    nfail = nfail + 1;
                end
            end
            
            if nfail > 0
                assert(nfail < nsub)
                T_ = T_*(nsub/(nsub - nfail));
                X_ = X_*(nsub/(nsub - nfail));
                Y_ = Y_*(nsub/(nsub - nfail));
                Z_ = Z_*(nsub/(nsub - nfail));
                angle_ = angle_*(nsub/(nsub - nfail));
                unwrap_ = unwrap_*(nsub/(nsub - nfail));
            end

            %% write summary averages

            this.out_dir = fullfile(getenv('SINGULARITY_HOME'), 'AnalyticSignalHCPAging');
            if ~isempty(nlim)
                tags = this.tags(sprintf("rsn%i-nlim%i", rsn, nlim));
            end

            this.write_ciftis( ...
                T_, ...
                sprintf('var_T_as_sub-all_ses-all_%s', tags), ...
                partitions=[], ...
                do_save_dynamic=false);
            this.write_ciftis( ...
                X_, ...
                sprintf('var_X_as_sub-all_ses-all_%s', tags), ...
                partitions=[], ...
                do_save_dynamic=false);
            this.write_ciftis( ...
                Y_, ...
                sprintf('var_Y_as_sub-all_ses-all_%s', tags), ...
                partitions=[], ...
                do_save_dynamic=false);
            this.write_ciftis( ...
                Z_, ...
                sprintf('var_Z_as_sub-all_ses-all_%s', tags), ...
                partitions=[], ...
                do_save_dynamic=false);
            this.write_ciftis( ...
                angle_, ...
                sprintf('var_angle_as_sub-all_ses-all_%s', tags), ...
                partitions=[], ...
                do_save_dynamic=false);
            this.write_ciftis( ...
                unwrap_, ...
                sprintf('var_unwrap_as_sub-all_ses-all_%s', tags), ...
                partitions=[], ...
                do_save_dynamic=false);
        end

        %% running call on single server

        function server_call(cores, opts)
            %% for servers

            arguments
                cores {mustBeScalarOrEmpty} = 2
                opts.N_sub {mustBeScalarOrEmpty} = 689
                opts.flip_globbed logical = true
            end

            % root_dir = '/home/usr/jjlee/mnt/CHPC_hcpdb/packages/unzip/HCP_1200';
            root_dir = fullfile(getenv('SINGULARITY_HOME'), 'HCPAging', 'HCPAgingRec', 'fmriresults01');

            g = glob(fullfile(root_dir, 'HCA*'));
            g = strip(g, filesep);
            g = flip(g); % examine more recent ones first
            g = cellfun(@(x) basename(x), g, UniformOutput=false);
            if opts.flip_globbed
                g = flip(g); % examine more recent ones
            end
            g = g(1:opts.N_sub);
            leng = length(g);
            %for idxg = 1:1
            parfor (idxg = 1:leng, cores)
                try
                    this = mlraut.AnalyticSignalHCPAgingPar( ...
                        subjects=g(idxg), ...
                        do_7T=true, ...
                        do_resting=true, ...
                        do_task=false, ...
                        do_global_signal_regression=true, ...
                        do_save=true, ...
                        do_save_dynamic=false, ...
                        do_save_ciftis=false, ...
                        hp_thresh=0.01, ...
                        lp_thresh=0.1, ...
                        v_physio=50, ...
                        plot_range=1:250, ...
                        tags="AnalyticSignalHCPAgingPar-parcall");
                    call(this);
                catch ME
                    handwarning(ME)
                end
            end
        end

        %% running call on cluster

        function [j,c] = cluster_construct_and_call(globbing_mat, opts)
            %% for clusters running Matlab parallel server

            arguments
                globbing_mat {mustBeFile} = ...
                    fullfile( ...
                    getenv('SINGULARITY_HOME'), 'HCPAging', 'HCPAgingRec', 'fmriresults01', ...
                    'mlraut_AnalyticSignalHCPAgingPar_globbing.mat')
                opts.sub_indices double = []  % total ~ 1:725
                opts.globbing_var = "globbed"
                opts.Ncol {mustBeInteger} = 8
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

            c = mlraut.CHPC3.propcluster('joshua_shimony', mempercpu='64gb', walltime='01:00:00');
            disp(c.AdditionalProperties)
            for irow = 1:Nrow
                try
                    j = c.batch( ...
                        @mlraut.AnalyticSignalHCPAgingPar.construct_and_call, ...
                        1, ...
                        {globbed(irow, :)}, ...
                        'Pool', opts.Ncol, ...
                        'CurrentFolder', '/scratch/jjlee/Singularity/AnalyticSignalHCPAging', ...
                        'AutoAddClientPath', false);
                catch ME
                    handwarning(ME)
                end
            end

            % j = c.batch(@mlraut.AnalyticSignalHCPAgingPar.construct_and_call, 1, {}, 'CurrentFolder', '.', 'AutoAddClientPath', false);
        end

        function durations = construct_and_call(subjects, opts)
            arguments
                subjects cell = {'HCA6002236_V1_MR'}
                opts.tasks cell = {'fMRI_CONCAT_ALL'}
                opts.tags {mustBeTextScalar} = "ASHCPAgingPar"
                opts.out_dir {mustBeFolder} = "/scratch/jjlee/Singularity/AnalyticSignalHCPAging"
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
                    as = mlraut.AnalyticSignalHCPAgingPar( ...
                        subjects=subjects{sidx}, ...
                        tasks=opts.tasks, ...
                        do_save=true, ...
                        do_save_subset=true, ...
                        hp_thresh=0.01, ...
                        lp_thresh=0.1, ...
                        out_dir=opts.out_dir, ...
                        source_physio=[ ...
                        "iFV-brightest", "iFV", "iFV-quantile", "sFV", "3rdV", "latV", "centrumsemiovale", ...
                        "ctx", "cbm", "str", "thal"], ...
                        tags=opts.tags);
                    call(as);
                catch ME
                    handwarning(ME)
                end

                durations(sidx) = toc;
            end
        end

        %% utilities

        function cdata = sample_rsn(cdata, parts)
            arguments
                cdata
                parts logical = true(size(cdata, 1), 1)
            end

            %% Nt x Nx => 1 x Nx

            cdata = cdata(parts, :);  % select t interesting
            cdata = mean(cdata, 1);  % average over t
        end
    end

    methods
        function this = AnalyticSignalHCPAgingPar(varargin)
            this = this@mlraut.AnalyticSignalHCPAging(varargin{:});
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
