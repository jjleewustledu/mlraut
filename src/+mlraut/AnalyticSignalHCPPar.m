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
                nmats {mustBeNumeric} = inf  % use finite for validity checks
            end

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
            if nmats > length(mats)
                nmats = length(mats);
            end
            mats = mats(1:nmats);
            nbin = this.cifti.Nbins;
            ngo = this.num_nodes;
            X_ = zeros(nbin, ngo);
            Y_ = zeros(nbin, ngo);
            Z_ = zeros(nbin, ngo);
            T_ = zeros(nbin, ngo);
            r_ = zeros(nbin, ngo);
            bold_ = zeros(nbin, ngo);
            plvs_ = zeros(nbin, ngo);
            errs = 0;

            % abbrev.
            bin = @this.bin_by_physio_angle;

            for mat = mats
                tic
                ld = load(mat{1});
                psi = ld.this_subset.bold_signal;
                phi = ld.this_subset.physio_signal;
                try
                    X = bin(this.X(psi, phi), phi);
                    Y = bin(this.Y(psi, phi), phi);
                    Z = bin(this.Z(psi, phi), phi);
                    T = bin(this.T(psi, phi), phi);
                    r = bin(this.connectivity(psi, phi), phi);
                    bold = bin(psi, phi);
                    plvs = bin(this.phase_locked_values(psi, phi), phi);

                    X_ = X_ + real(X/nmats);
                    Y_ = Y_ + real(Y/nmats);
                    Z_ = Z_ + real(Z/nmats);
                    T_ = T_ + real(T/nmats);
                    r_ = r_ + real(r/nmats);
                    bold_ = bold_ + real(bold/nmats);
                    plvs_ = plvs_ + real(plvs/nmats);
                catch ME
                    fprintf("while working with %s: ", mat{1});
                    handwarning(ME)
                    errs = errs + 1;
                end
                fprintf("time working with %s: ", mat{1});
                toc
            end
            fprintf("errs->%g\n", errs)

            dt = 2*pi/this.cifti_.Nbins;
            units_t = "RADIAN";
            this.cifti.write_cifti( ...
                X_, sprintf('X_as_sub-all_ses-all_%s', this.tags), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                Y_, sprintf('Y_as_sub-all_ses-all_%s', this.tags), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                Z_, sprintf('Z_as_sub-all_ses-all_%s', this.tags), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                T_, sprintf('T_as_sub-all_ses-all_%s', this.tags), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                r_, sprintf('comparator_as_sub-all_ses-all_%s', this.tags), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                bold_, sprintf('bold_as_sub-all_ses-all_%s', this.tags), dt=dt, units_t=units_t);
            this.cifti.write_cifti( ...
                plvs_, sprintf('plvs_as_sub-all_ses-all_%s', this.tags), dt=dt, units_t=units_t);
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
                opts.sub_indices double = 1:916  % total ~ 1:1113
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

            c = mlraut.CHPC3.propcluster(mempercpu='40gb', walltime='00:30:00');
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
    
    end

    methods
        function this = AnalyticSignalHCPPar(varargin)
            this = this@mlraut.AnalyticSignalHCP(varargin{:});
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
