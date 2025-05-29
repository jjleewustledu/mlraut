classdef AnalyticSignalHCPPar < handle & mlraut.AnalyticSignalHCP
    %% line1
    %  line2
    %  
    %  Created 13-Apr-2023 02:11:52 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 9.14.0.2206163 (R2023a) for MACI64.  Copyright 2023 John J. Lee.
    
    methods (Static)
        function ihcp = mean_bold()
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
                net_sigs2 = real(net_sigs2);
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

            ihcp.write_cifti(bins3, stackstr() + ".dtseries.nii");

            popd(pwd0);
        end

        function this = mean_twistor()
            this = mlraut.AnalyticSignalHCPPar( ...
                subjects={'995174'}, ...
                do_7T=false, ...
                do_resting=true, ...
                do_task=false, ...
                do_global_signal_regression=true, ...
                do_save=false, ...
                do_save_dynamic=false, ...
                do_save_ciftis=false, ...
                hp_thresh=0.01, ...
                lp_thresh=0.1, ...
                tags="AnalyticSignalHCPPar-mean-twistor");

            mats = asrow(glob(fullfile(this.out_dir, '*/sub-*_ses-*AnalyticSignalHCP*.mat')));
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
            %% 33 GiB memory needed per instance of this running on a single process

            arguments
                cores {mustBeScalarOrEmpty} = 8
                opts.N_sub {mustBeScalarOrEmpty} = 1113
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
            g = g(1:opts.N_sub);
            leng = length(g);
            %for idxg = 1:1
            %parfor (idxg = 1:2, 2)
            parfor (idxg = 1:leng, cores)
                try
                    % if isfolder(fullfile(out_dir, g{idxg}))
                    %     continue
                    % end
                    this = mlraut.AnalyticSignalHCPPar( ...
                        subjects=g(idxg), ...  
                        do_7T=false, ...
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
                        tags="AnalyticSignalHCPPar-parcall");
                    call(this);
                catch ME
                    handwarning(ME)
                end
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
