classdef AnalyticSignalHCPAgingPar < handle & mlraut.AnalyticSignalHCPAging
    %% line1
    %  line2
    %  
    %  Created 07-Feb-2024 23:37:28 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    
    methods (Static)
        function ret = mean_twistor_comparator(nlim)
            arguments
                nlim = 4;
            end

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

            mats = flip(asrow(mglob(fullfile(this.out_dir, "HCA*_MR/sub-*_ses-*-iFV-AnalyticSignalHCPAgingPar.mat"))));
            mats = mats(1:nlim);
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
            tags = this.tags(sprintf("nlim%i", nlim));

            this.write_ciftis( ...
                comparator_, ...
                sprintf('mean_comparator_as_sub-all_ses-all_%s', tags));
        end

        function ret = mean_element(nlim, ele)
            arguments
                nlim = 4
                ele {mustBeTextScalar} = "T"
            end

            assert(any(contains(["T", "X", "Y", "Z", "angle", "unwrap"], ele)))

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

            mats = flip(asrow(mglob(fullfile(this.out_dir, "HCA*_MR/sub-*_ses-*-iFV-AnalyticSignalHCPAgingPar.mat"))));
            mats = mats(1:nlim);
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
                        this.(ele)(ld.this.bold_signal, ld.this.physio_signal));
                    E_ = E_ + E_norm/nsub;
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
            tags = this.tags(sprintf("nlim%i", nlim));

            this.write_ciftis( ...
                E_, ...
                sprintf("mean_%s_as_sub-all_ses-all_%s", ele, tags), ...
                partitions=[], ...
                do_save_dynamic=false);
        end

        function ret = mean_twistor(nlim, rsn)
            arguments
                nlim = 4
                rsn = 7
            end

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

            mats = flip(asrow(mglob(fullfile(this.out_dir, "HCA*_MR/sub-*_ses-*-iFV-AnalyticSignalHCPAgingPar.mat"))));
            mats = mats(1:nlim);
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
            tags = this.tags(sprintf("rsn%i-nlim%i", rsn, nlim));

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

        function ret = mean_twistor_Y(nlim, rsn)
            arguments
                nlim = 4
                rsn = 7
            end

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

            mats = flip(asrow(mglob(fullfile(this.out_dir, "HCA*_MR/sub-*_ses-*-iFV-AnalyticSignalHCPAgingPar.mat"))));
            mats = mats(1:nlim);
            nsub = length(mats);
            nx = this.num_nodes;
            nfail = 0;

            Y_ = zeros(1, nx);

            fprintf(stackstr() + "\n");
            for mat = mats
                try
                    fprintf("loading %s\n", mat);
                    tic

                    ld = load(mat);

                    ctx = ld.this.HCP_signals.ctx;
                    angle_rsn = this.angle(ctx.psi(:,rsn), ctx.phi(:,rsn));
                    u_interesting = sin(angle_rsn) > 0;

                    Y_rsn = this.sample_rsn( ...
                        this.Y(ld.this.bold_signal, ld.this.physio_signal), u_interesting);

                    Y_ = Y_ + Y_rsn/nsub;
                    ret = ret + 1;

                    toc
                catch ME
                    handwarning(ME)
                    nfail = nfail + 1;
                end
            end
            
            if nfail > 0
                assert(nfail < nsub)
                Y_ = Y_*(nsub/(nsub - nfail));
            end

            %% write summary averages

            this.out_dir = fullfile(getenv('SINGULARITY_HOME'), 'AnalyticSignalHCPAging');
            tags = this.tags(sprintf("rsn%i-nlim%i", rsn, nlim));

            this.write_ciftis( ...
                Y_, ...
                sprintf('mean_Y_as_sub-all_ses-all_%s', tags), ...
                partitions=[], ...
                do_save_dynamic=false);
        end

        function ret = var_element(nlim, ele)
            arguments
                nlim = 4
                ele {mustBeTextScalar} = "T"
            end

            assert(any(contains(["T", "X", "Y", "Z", "angle", "unwrap"], ele)))

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

            tags_ = this.tags(sprintf('nlim%i', nlim));
            %tags_ = this.tags(sprintf('nlim%i', 725));  % DEBUGGING
            mu = cifti_read( ...
                fullfile(char(this.out_dir), ...
                sprintf('mean_%s_as_sub-all_ses-all_%s_avgt.dscalar.nii', char(ele), char(tags_))));

            mats = flip(asrow(mglob(fullfile(this.out_dir, "HCA*_MR/sub-*_ses-*-iFV-AnalyticSignalHCPAgingPar.mat"))));
            mats = mats(1:nlim);
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
                        this.(ele)(ld.this.bold_signal, ld.this.physio_signal));
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
            tags = this.tags(sprintf("nlim%i", nlim));

            this.write_ciftis( ...
                E_, ...
                sprintf("var_%s_as_sub-all_ses-all_%s", ele, tags), ...
                partitions=[], ...
                do_save_dynamic=false);
        end

        function ret = var_twistor_comparator(nlim)
            arguments
                nlim = 725;
            end

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

            tags_ = this.tags(sprintf('nlim%i', nlim));
            %tags_ = this.tags(sprintf('nlim%i', 725));  % DEBUGGING
            mu = cifti_read( ...
                fullfile(this.out_dir, sprintf('mean_comparator_as_sub-all_ses-all_%s_avgt.dscalar.nii', tags_)));

            mats = flip(asrow(mglob(fullfile(this.out_dir, "HCA*_MR/sub-*_ses-*-iFV-AnalyticSignalHCPAgingPar.mat"))));
            mats = mats(1:nlim);
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
            tags = this.tags(sprintf("nlim%i", nlim));

            this.write_ciftis( ...
                comparator_, ...
                sprintf('var_comparator_as_sub-all_ses-all_%s', tags));
        end

        function ret = var_twistor(nlim, rsn)
            arguments
                nlim = 4
                rsn = 7
            end

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

            tags_ = this.tags(sprintf('rsn%i-nlim%i', rsn, nlim));
            %tags_ = this.tags(sprintf('rsn%i-nlim%i', rsn, 725));  % DEBUGGING
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

            mats = flip(asrow(mglob(fullfile(this.out_dir, "HCA*_MR/sub-*_ses-*-iFV-AnalyticSignalHCPAgingPar.mat"))));
            mats = mats(1:nlim);
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
            tags = this.tags(sprintf("rsn%i-nlim%i", rsn, nlim));

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

        function ret = var_twistor_Y(nlim, rsn)
            arguments
                nlim = 4
                rsn = 7
            end

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

            tags_ = this.tags(sprintf('rsn%i-nlim%i', rsn, nlim));
            %tags_ = this.tags(sprintf('rsn%i-nlim%i', rsn, 725));  % DEBUGGING
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

            mats = flip(asrow(mglob(fullfile(this.out_dir, "HCA*_MR/sub-*_ses-*-iFV-AnalyticSignalHCPAgingPar.mat"))));
            mats = mats(1:nlim);
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
            tags = this.tags(sprintf("rsn%i-nlim%i", rsn, nlim));

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

        %% running {mean,var}_* on cluster

        function [j,c] = cluster_batch_stats(do_var, opts)
            %% for clusters running Matlab parallel server

            arguments
                do_var logical = false
                opts.nlim double = 4
                opts.only_comparator logical = false
                opts.include_comparator logical = false
            end
            if opts.only_comparator
                opts.include_comparator = true;
            end

            c = mlraut.CHPC3.propcluster_16gb_100h();
            disp(c.AdditionalProperties)

            %% T,X,Y,Z,angle,unwrap

            if ~opts.only_comparator
                if ~do_var
                    for rsn = 1:7
                        try
                            j = c.batch( ...
                                @mlraut.AnalyticSignalHCPAgingPar.mean_twistor, ...
                                1, ...
                                {opts.nlim, rsn}, ...
                                'CurrentFolder', '.', ...
                                'AutoAddClientPath', false);
                        catch ME
                            handwarning(ME)
                        end
                    end
                else
                    for rsn = 1:7
                        try
                            j = c.batch( ...
                                @mlraut.AnalyticSignalHCPAgingPar.var_twistor, ...
                                1, ...
                                {opts.nlim, rsn}, ...
                                'CurrentFolder', '.', ...
                                'AutoAddClientPath', false);
                        catch ME
                            handwarning(ME)
                        end
                    end
                end
            end

            %% comparator

            if opts.include_comparator
                if ~do_var
                    try
                        j = c.batch( ...
                            @mlraut.AnalyticSignalHCPAgingPar.mean_twistor_comparator, ...
                            1, ...
                            {opts.nlim}, ...
                            'CurrentFolder', '.', ...
                            'AutoAddClientPath', false);
                    catch ME
                        handwarning(ME)
                    end
                else
                    try
                        j = c.batch( ...
                            @mlraut.AnalyticSignalHCPAgingPar.var_twistor_comparator, ...
                            1, ...
                            {opts.nlim}, ...
                            'CurrentFolder', '.', ...
                            'AutoAddClientPath', false);
                    catch ME
                        handwarning(ME)
                    end
                end
            end
        end

        function [j,c] = cluster_batch_stats_2(do_var, opts)
            %% for clusters running Matlab parallel server

            arguments
                do_var logical = false
                opts.nlim double = 4
                opts.elements {mustBeText} = ["X", "Y"]  % ["T", "Z", "angle", "unwrap"]
            end

            c = mlraut.CHPC3.propcluster_tiny();
            disp(c.AdditionalProperties)

            %% X,Y; or T,Z,angle,unwrap

            if ~do_var
                for e = opts.elements
                    try
                        j = c.batch( ...
                            @mlraut.AnalyticSignalHCPAgingPar.mean_element, ...
                            1, ...
                            {opts.nlim, e}, ...
                            'CurrentFolder', '.', ...
                            'AutoAddClientPath', false);
                    catch ME
                        handwarning(ME)
                    end
                end
            else
                for e = opts.elements
                    try
                        j = c.batch( ...
                            @mlraut.AnalyticSignalHCPAgingPar.var_element, ...
                            1, ...
                            {opts.nlim, e}, ...
                            'CurrentFolder', '.', ...
                            'AutoAddClientPath', false);
                    catch ME
                        handwarning(ME)
                    end
                end
            end
        end

        function [j,c] = cluster_batch_stats_3(do_var, opts)
            %% for clusters running Matlab parallel server

            arguments
                do_var logical = false
                opts.nlim double = 4
                opts.elements {mustBeText} = 7  % rsn
            end

            c = mlraut.CHPC3.propcluster_tiny();
            disp(c.AdditionalProperties)

            %% T,Z,angle,unwrap

            if ~do_var
                for e = opts.elements
                    try
                        j = c.batch( ...
                            @mlraut.AnalyticSignalHCPAgingPar.mean_twistor_Y, ...
                            1, ...
                            {opts.nlim, e}, ...
                            'CurrentFolder', '.', ...
                            'AutoAddClientPath', false);
                    catch ME
                        handwarning(ME)
                    end
                end
            else
                for e = opts.elements
                    try
                        j = c.batch( ...
                            @mlraut.AnalyticSignalHCPAgingPar.var_twistor_Y, ...
                            1, ...
                            {opts.nlim, e}, ...
                            'CurrentFolder', '.', ...
                            'AutoAddClientPath', false);
                    catch ME
                        handwarning(ME)
                    end
                end
            end
        end

        %% running call on single server or cluster

        function server_call(cores, opts)
            %% for servers

            arguments
                cores {mustBeScalarOrEmpty} = 2
                opts.N_sub {mustBeScalarOrEmpty} = 725
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
                        do_save=true, ...
                        do_save_dynamic=false, ...
                        do_save_ciftis=false, ...
                        hp_thresh=0.01, ...
                        lp_thresh=0.05, ...
                        v_physio=50, ...
                        plot_range=1:250, ...
                        global_signal_regression=true, ...
                        tags="AnalyticSignalHCPAgingPar-parcall");
                    call(this);
                catch ME
                    handwarning(ME)
                end
            end
        end

        function [j,c] = cluster_batch_call(globbing_mat, opts)
            %% for clusters running Matlab parallel server

            arguments
                globbing_mat {mustBeFile} = ...
                    fullfile( ...
                    getenv('SINGULARITY_HOME'), 'HCPAging', 'HCPAgingRec', 'fmriresults01', ...
                    'mlraut_AnalyticSignalHCPAgingPar_globbing.mat')
                opts.sub_indices double = [1,3:100]
                opts.globbing_var = "globbed"
            end
            ld = load(globbing_mat);
            globbed = convertStringsToChars(ld.(opts.globbing_var));
            globbed = asrow(globbed);
            globbed = globbed(opts.sub_indices);

            c = mlraut.CHPC3.propcluster();
            disp(c.AdditionalProperties)
            for g = globbed
                try
                    j = c.batch( ...
                        @mlraut.AnalyticSignalHCPAgingPar.construct_and_call, ...
                        1, ...
                        {g(1)}, ...
                        'CurrentFolder', '.', ...
                        'AutoAddClientPath', false);
                catch ME
                    handwarning(ME)
                end
            end

            % j = c.batch(@mlraut.AnalyticSignalHCPAgingPar.construct_and_call, 1, {}, 'CurrentFolder', '.', 'AutoAddClientPath', false);
        end

        function duration = construct_and_call(subjects, opts)
            arguments
                subjects cell = {'HCA6002236_V1_MR'}
                opts.tasks cell = {'fMRI_CONCAT_ALL'}
                opts.tags {mustBeTextScalar} = "AnalyticSignalHCPAgingPar"
                opts.out_dir {mustBeFolder} = "/scratch/jjlee/Singularity/AnalyticSignalHCPAging"
            end
            
            mlraut.CHPC3.setenvs();
            setenv("VERBOSITY", "1");
            ensuredir(opts.out_dir);
            ensuredir(fullfile(opts.out_dir, subjects{1}));
            diary(fullfile(opts.out_dir, subjects{1}, "diary.log"));
            tic;
            disp("constructing mlraut.AnalyticSignalHCPAgingPar")
            this = mlraut.AnalyticSignalHCPAgingPar( ...
                subjects=subjects, ...
                tasks=opts.tasks, ...
                do_7T=false, ...
                do_plot_networks=false, ...
                do_resting=true, ...
                do_task=false, ...
                do_save=true, ...
                do_save_dynamic=false, ...
                do_save_ciftis=false, ...
                do_save_subset=false, ...
                force_band=false, ...
                global_signal_regression=true, ...
                hp_thresh=0.01, ...
                lp_thresh=0.05, ...
                out_dir=opts.out_dir, ...
                source_physio="iFV", ...
                tags=opts.tags, ...                
                v_physio=50);
            disp("calling this")
            call(this);
            duration = toc;
            fprintf("tic-toc duration: %s seconds", duration);
            diary("off");
        end

        function [j,c] = cluster_par_call(globbing_mat, opts)
            %% for clusters running Matlab parallel server

            arguments
                globbing_mat {mustBeFile} = ...
                    fullfile( ...
                    getenv('SINGULARITY_HOME'), 'HCPAging', 'HCPAgingRec', 'fmriresults01', ...
                    'mlraut_AnalyticSignalHCPAgingPar_globbing.mat')
                opts.sub_indices double = 101:200
                opts.N_max double = 1000
                opts.globbing_var = "globbed"
            end
            ld = load(globbing_mat);
            globbed = convertStringsToChars(ld.(opts.globbing_var));
            globbed = asrow(globbed);
            globbed = globbed(opts.sub_indices);
            N_sub = length(opts.sub_indices);

            c = mlraut.CHPC3.propcluster();
            disp(c.AdditionalProperties)
            try
                j = c.batch( ...
                    @mlraut.AnalyticSignalHCPAgingPar.par_construct_and_call, ...
                    1, ...
                    {globbed}, ...
                    'Pool', min(N_sub, opts.N_max), ...
                    'CurrentFolder', '.', ...
                    'AutoAddClientPath', false);
            catch ME
                handwarning(ME)
            end
        end

        function duration = par_construct_and_call(subjects, opts)
            arguments
                subjects cell = {'HCA6002236_V1_MR'}
                opts.tasks cell = {'fMRI_CONCAT_ALL'}
                opts.tags {mustBeTextScalar} = "AnalyticSignalHCPAgingPar"
                opts.out_dir {mustBeFolder} = "/scratch/jjlee/Singularity/AnalyticSignalHCPAging"
            end
            
            duration = [];

            parfor sidx = 1:length(subjects)            
                mlraut.CHPC3.setenvs();
                ensuredir(opts.out_dir); %#ok<PFBNS>
                setenv("VERBOSITY", "1");
                ensuredir(fullfile(opts.out_dir, subjects{sidx})); 
                diary(fullfile(opts.out_dir, subjects{sidx}, "diary.log")); 
                tic;
                disp("constructing mlraut.AnalyticSignalHCPAgingPar")
                this = mlraut.AnalyticSignalHCPAgingPar( ...
                    subjects=subjects(sidx), ...
                    tasks=opts.tasks, ...
                    do_7T=false, ...
                    do_plot_networks=false, ...
                    do_resting=true, ...
                    do_task=false, ...
                    do_save=true, ...
                    do_save_dynamic=false, ...
                    do_save_ciftis=false, ...
                    do_save_subset=false, ...
                    force_band=false, ...
                    global_signal_regression=true, ...
                    hp_thresh=0.01, ...
                    lp_thresh=0.05, ...
                    out_dir=opts.out_dir, ...
                    source_physio="iFV", ...
                    tags=opts.tags, ...
                    v_physio=50);
                disp("calling this")
                call(this);
                fprintf("tic-toc duration: %s seconds", duration);
                diary("off");
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
