classdef AnalyticSignalHCPAgingPar < handle & mlraut.AnalyticSignalHCPAging
    %% line1
    %  line2
    %  
    %  Created 07-Feb-2024 23:37:28 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    
    methods (Static)
        function this = median_twistor()
            this = mlraut.AnalyticSignalHCPAgingPar( ...
                subjects={'HCA9992517_V1_MR'}, ...
                do_7T=false, ...
                do_resting=true, ...
                do_task=false, ...
                do_save=false, ...
                do_save_dynamic=false, ...
                do_save_ciftis=false, ...
                hp_thresh=0.01, ...
                lp_thresh=0.1, ...
                global_signal_regression=true, ...
                tags="AnalyticSignalHCPAgingPar-median-twistor");

            mats = asrow(glob(fullfile(this.out_dir, 'HCA*_MR/sub-*_ses-*AnalyticSignalHCPAging*.mat')));
            n = length(mats);
            nx = this.num_nodes;
            nfail = 0;

            X_ = zeros(1, nx);
            Y_ = zeros(1, nx);
            Z_ = zeros(1, nx);
            T_ = zeros(1, nx);

            fprintf(stackstr() + "\n");
            for mat = mats
                try
                    tic
                    fprintf("loading %s\n", mat{1});
                    ld = load(mat{1});
                    this_subset = ld.this_subset;
                    try
                        ksi = this_subset.bold_signal;
                    catch ME
                        if strcmp(ME.identifier, 'MATLAB:nonExistentField')
                            ksi = this_subset.analytic_signal.*this_subset.physio_signal;  % overly normalized in this_subset
                        else
                            rethrow(ME)
                        end
                    end
                    eta = this_subset.physio_signal;
                    X = median(this.build_final_normalization((ksi.*conj(eta) + eta.*conj(ksi))/sqrt(2)), 1);
                    Y = median(this.build_final_normalization((ksi.*conj(eta) - eta.*conj(ksi))/sqrt(2i)), 1);
                    Z = median(this.build_final_normalization((ksi.*conj(ksi) - eta.*conj(eta))/sqrt(2)), 1);
                    T = median(this.build_final_normalization((ksi.*conj(ksi) + eta.*conj(eta))/sqrt(2)), 1);

                    X_ = X_ + X/n;
                    Y_ = Y_ + Y/n;
                    Z_ = Z_ + Z/n;
                    T_ = T_ + T/n;
                    toc
                catch ME
                    handwarning(ME)
                    nfail = nfail + 1;
                end
            end
            
            if nfail > 0
                assert(nfail < n)
                X_ = X_*(n/(n - nfail));
                Y_ = Y_*(n/(n - nfail));
                Z_ = Z_*(n/(n - nfail));
                T_ = T_*(n/(n - nfail));
            end

            this.out_dir = fullfile(getenv('HOME'), 'Singularity', 'AnalyticSignalHCPAging_zfs');
            ensuredir(this.out_dir);

            this.write_ciftis( ...
                X_, sprintf('X_as_sub-all_ses-all_%s', this.tags), ...
                do_final_normalization=false, ...
                do_save_dynamic=false);
            this.write_ciftis( ...
                Y_, sprintf('Y_as_sub-all_ses-all_%s', this.tags), ...
                do_final_normalization=false, ...
                do_save_dynamic=false);
            this.write_ciftis( ...
                Z_, sprintf('Z_as_sub-all_ses-all_%s', this.tags), ...
                do_final_normalization=false, ...
                do_save_dynamic=false);
            this.write_ciftis( ...
                T_, sprintf('T_as_sub-all_ses-all_%s', this.tags), ...
                do_final_normalization=false, ...
                do_save_dynamic=false);
            this.write_ciftis( ...
                T_+Z_, sprintf('T+Z_as_sub-all_ses-all_%s', this.tags), ...
                do_final_normalization=false, ...
                do_save_dynamic=false);
        end

        function parcall(cores, opts)
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
                        do_task=true, ...
                        do_save=true, ...
                        do_save_dynamic=false, ...
                        do_save_ciftis=false, ...
                        hp_thresh=0.01, ...
                        lp_thresh=0.1, ...
                        global_signal_regression=true, ...
                        tags="AnalyticSignalHCPAgingPar-parcall");
                    call(this);
                catch ME
                    handwarning(ME)
                end
            end
        end

        function [j,c] = parcluster(globbing_mat)
            arguments
                globbing_mat {mustBeFile} = fullfile( ...
                    getenv('SINGULARITY_HOME'), 'HCPAging', 'HCPAgingRec', 'fmriresults01', ...
                    'mlraut_AnalyticSignalHCPAgingPar_globbing.mat')
            end
            ld = load(globbing_mat);
            globbed = asrow(ld.globbed);
            % globbed = globbed(1:3);

            mlraut.CHPC3.propcluster();
            c = parcluster;
            disp(c.AdditionalProperties)
            for g = globbed
                try
                    j = c.batch( ...
                        @mlraut.AnalyticSignalHCPAgingPar.construct_and_call, ...
                        1, ...
                        {'subjects', g(1)}, ...
                        'CurrentFolder', '.', ...
                        'AutoAddClientPath', false);
                catch ME
                    handwarning(ME)
                end
            end
        end

        function duration = construct_and_call(opts)
            arguments
                opts.subjects cell
                opts.tasks cell = {'fMRI_CONCAT_ALL'}
                opts.tags {mustBeTextScalar} = "AnalyticSignalHCPAgingPar"
            end

            % setenv('TMPDIR', '/scratch/jjlee/tmp') % worker nodesk
            setenv('SINGULARITY_HOME', '/scratch/jjlee/Singularity')
            tic
            this = mlraut.AnalyticSignalHCPAgingPar( ...
                subjects=opts.subjects, ...
                tasks=opts.tasks, ...
                do_save=false, ...
                do_save_ciftis=true, ...
                tags=opts.tags);
            call(this);
            duration = toc;
        end
    end

    methods
        function this = AnalyticSignalHCPAgingPar(varargin)
            this = this@mlraut.AnalyticSignalHCPAging(varargin{:});
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
