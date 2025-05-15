classdef FultzMulti < handle
    %% line1
    %  line2
    %  
    %  Created 06-May-2025 12:28:48 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 25.1.0.2897550 (R2025a) Prerelease Update 5 for MACA64.  Copyright 2025 John J. Lee.
    
    properties
        anat_weights
        asobj
        phi_multi  % may be iFV-brightest, other ventricles, brainstem, CE_on_T1w
        psi_multi
        regions = ["ctx", "cbm", "str", "thal"];  % regions to plot separately
        waves_dir = fullfile( ...
            getenv("HOME"), "MATLAB-Drive", "arousal-waves-main")
    end

    properties (Dependent)
        cbm_mask
        ctx_mask
        Fs
        hcp_signals
        out_dir
        source_physio
        str_mask
        thal_mask
        tr
    end

    methods  %% GET

        function g = get.cbm_mask(this)
            ld = load(fullfile(this.waves_dir, 'supporting_files', 'mask_cbm_HCP.mat'));
            g = ld.mask_cbm;
        end

        function g = get.ctx_mask(this)
            ld = load(fullfile(this.waves_dir, 'supporting_files', 'mask_ctx_HCP.mat'));
            g = ld.mask_ctx;
        end

        function g = get.Fs(this)
            g = this.asobj.Fs;
        end

        function g = get.hcp_signals(this)
            g = this.asobj.HCP_signals;
        end

        function g = get.out_dir(this)
            g = strrep(fileparts(this.asobj.out_dir), ...
                "/scratch/jjlee/Singularity", getenv("SINGULARITY_HOME"));
        end

        function g = get.source_physio(this)
            g = this.asobj.source_physio;
        end

        function g = get.str_mask(this)
            ld = load(fullfile(this.waves_dir, 'supporting_files', 'mask_str_HCP.mat'));
            g = ld.mask_str;
        end

        function g = get.thal_mask(this)
            ld = load(fullfile(this.waves_dir, 'supporting_files', 'mask_thal_HCP.mat'));
            g = ld.mask_thal;
        end
        
        function g = get.tr(this)
            g = this.asobj.tr;
        end
    end

    methods
        function this = FultzMulti(aso, opts)
            arguments
                aso = []
                opts.phi_multi = []
                opts.psi_multi = []
            end

            this.asobj = aso;
            this.phi_multi = opts.phi_multi;            
            this.psi_multi = opts.psi_multi;            

            % weights for use with AnalyticSignals.HCP_signals
            nets = mlraut.NetworkCollections(this.asobj, []);
            this.anat_weights = nets.anat_weights;
        end

        function this = add_phi_psi_from_aso(this)

            % preallocate large arrays
            Nt = size(this.asobj.bold_signal, 1);
            Nsub = 1;
            Nr = numel(this.regions);
            this.phi_multi = zeros(Nt, Nsub, Nr);
            this.psi_multi = zeros(Nt, Nsub, Nr);

            % load and assign phi, psi
            failures = 0;
            ig = 1;
            try
                for ir = 1:Nr
                    region = this.regions(ir);
                    signal = this.asobj.HCP_signals.(region);
                    this.phi_multi(:, ig, ir) = this.build_weighted(signal.phi, region);
                    this.psi_multi(:, ig, ir) = this.build_weighted(signal.psi, region);
                end
            catch ME
                fprintf("%s: fault in %s\n", ME.identifier, g);
                failures = failures + 1;
            end

            this.phi_multi = this.phi_multi(:,end-failures,:);
            this.psi_multi = this.psi_multi(:,end-failures,:);
        end

        function [phi,psi] = build_tseries_phi_psi_from_mat(this, opts)
            arguments
                this mlraut.FultzMulti
                opts.mat_pattern {mustBeTextScalar} = fullfile("**", "*" + this.source_physio + "-scaleiqr*.mat")
                opts.tag {mustBeTextScalar} = this.source_physio
            end
            if ~isemptytext(opts.tag)
                opts.tag = "_" + opts.tag;
            end

            pwd0 = pushd(this.out_dir);  % folder containing subject folders containing *.mat

            % glob *.mat
            globbed = mglob(opts.mat_pattern); 
            globbed = this.cull_globbed(globbed);  % there may be bold1, bold2, all, and "concat"
            
            % preallocate large arrays
            Nt = size(this.asobj.bold_signal, 1);
            Nsub = numel(globbed);
            Nr = numel(this.regions);
            this.phi_multi = zeros(Nt, Nsub, Nr);
            this.psi_multi = zeros(Nt, Nsub, Nr);

            % load and assign phi, psi
            failures = 0;
            ig = 1;
            for g = asrow(globbed)
                ld = load(g);
                try
                    for ir = 1:Nr
                        region = this.regions(ir);
                        signal = ld.this.HCP_signals.(region);
                        this.phi_multi(:, ig, ir) = this.build_weighted(signal.phi, region);
                        this.psi_multi(:, ig, ir) = this.build_weighted(signal.psi, region);
                    end
                    ig = ig + 1;
                catch ME
                    fprintf("%s: fault in %s\n", ME.identifier, g);
                    failures = failures + 1;
                end
            end

            this.phi_multi = this.phi_multi(:,end-failures,:);
            this.psi_multi = this.psi_multi(:,end-failures,:);

            % save to filesystem
            phi = this.phi_multi;
            save(fullfile(this.out_dir, stackstr() + opts.tag + "_phi.mat"), "phi", "-v7.3");
            psi = this.psi_multi;
            save(fullfile(this.out_dir, stackstr() + opts.tag + "_psi.mat"), "psi", "-v7.3");

            popd(pwd0);
        end

        function psi = build_weighted(this, psi, region)
            %  psi ~ Nt x 9 double complex
            %  region ~ text
            %  this.anat_weights ~ 4 x struct; foreach struct field ~ 1 x 9 double

            w = this.anat_weights;  % readability
            switch char(region)
                case 'ctx'
                    psi = psi * ascol(w.ctx);
                case 'cbm'
                    psi = psi * ascol(w.cbm);
                case 'str'
                    psi = psi * ascol(w.str);
                case 'thal'
                    psi(isnan(psi)) = 0;
                    psi = psi * ascol(w.thal);
                otherwise
                    error("mlraut:RuntimeError", stackstr())
            end
        end

        function g = cull_globbed(this, g)
            if isa(this.asobj, "mlraut.AnalyticSignalGBM")
                g = g(contains(g, "run-all"));
            end
        end

        function h = plot_coherencyc(this, opts)
            arguments
                this mlraut.FultzMulti
                opts.tseries {mustBeTextScalar} = "-dbold/dt"
                opts.physio_tag {mustBeTextScalar} = this.source_physio
                opts.physio_tag2 {mustBeTextScalar} = this.source_physio
            end

            % aufbau complex phi, psi ~ Nt x Nsub x Nr
            tseries_phi = this.phi_multi;
            tseries_psi = this.psi_multi;
            switch char(opts.tseries)
                case 'bold'
                    tseries_phi = real(tseries_phi);
                    tseries_psi = real(tseries_psi);
                case '-dbold/dt'
                    tseries_phi = real(tseries_phi(1:end-1, :, :));
                    bold = real(tseries_psi);
                    tseries_psi = -diff(bold, 1, 1);  % 1st deriv. of time
                    tseries_psi(tseries_psi < 0) = 0;  % see Fultz et al. 2019
                case 'X'
                    tseries_psi = this.asobj.X(tseries_psi, tseries_phi);
                case 'Y'
                    tseries_psi = this.asobj.Y(tseries_psi, tseries_phi);
                case 'Z'
                    tseries_psi = this.asobj.Z(tseries_psi, tseries_phi);
            end

            % prepare chronux params
            Nr = numel(this.regions);
            T = size(tseries_phi, 1) * this.tr;
            W = this.asobj.lp_thresh;  % bandwidth
            Ntapers = 30;
            p = floor(2 * T * W - Ntapers);
            pval = 0.05;
            % ...
            params.Fs = this.Fs;
            params.fpass = [0, W];
            params.trialave = 1;  % T/F
            params.tapers = [W, T, p];
            params.err = [2, pval];  % jack-knife resampling
            trial = coherencyc(tseries_psi(:,:,1), tseries_phi(:,:,1), params);
            Nt = size(trial, 1); % chronux calls its getfgrid()

            % call chronux per region
            coh = zeros(Nt, Nr);
            pha = zeros(Nt, Nr);
            S12 = zeros(Nt, Nr);
            S1 = zeros(Nt, Nr);
            S2 = zeros(Nt, Nr);
            f = zeros(1, Nt);  % chronux peculiarity
            conf_coh = zeros(1, Nr);
            pha_std = zeros(Nt, Nr);
            coh_err = zeros(2, Nt, Nr);
            for ir = 1:Nr
                [coh(:,ir),pha(:,ir),S12(:,ir),S1(:,ir),S2(:,ir),f,conf_coh(ir),pha_std(:,ir),coh_err(:,:,ir)] = ...
                    coherencyc(tseries_psi(:,:,ir), tseries_phi(:,:,ir), params);
            end
            f = ascol(f);

            % plotting config.
            [~,f_0p05_idx] = min(abs(f - 0.05));  % reference for power spectra
            pha = unwrap(pha);  % readable phases
            assert(numel(this.regions) == Nr)
            
            % https://www.mathworks.com/help/matlab/creating_plots/combine-multiple-plots.html            
            h = tiledlayout(5, 1);

            ax1 = nexttile;
            hold("on");
            for ir = 1:Nr
                plot(ax1, f, coh(:, ir));
            end            
            hold("off");
            ylabel(ax1, "coherence")
            legend(ax1, this.regions)

            ax2 = nexttile;
            hold("on");
            for ir = 1:Nr
                plot(ax2, f, pha(:, ir));
            end
            hold("off")
            ylabel(ax2, "phase")

            ax3 = nexttile;
            hold("on");
            for ir = 1:Nr
                plot(ax3, f, 10 * log10(S12(:, ir)/S12(f_0p05_idx, ir)));
            end
            hold("off")
            ylabel(ax3, "power (dB, 0.05 Hz ref.)")
            title(ax3, "cross-spectra")

            ax4 = nexttile;
            hold("on");
            for ir = 1:Nr
                plot(ax4, f, 10 * log10(S1(:, ir)/S1(f_0p05_idx, ir)));
            end
            hold("off")
            ylabel(ax4, "power (dB, 0.05 Hz ref.)")
            title(ax4, opts.tseries + " spectrum")

            ax5 = nexttile;
            S2 = S2(:, 1);
            plot(ax5, f, 10 * log10(S2/S2(f_0p05_idx)), LineWidth=1.25, Color="m");
            ylabel(ax5, "power (dB, 0.05 Hz ref.)")
            title(ax5, opts.physio_tag + " spectrum", Interpreter="none")
            legend(ax5, opts.physio_tag2)

            linkaxes([ax1,ax2,ax3,ax4,ax5], 'x');
            xlabel(h, "frequency (Hz)");
            xticklabels(ax1,{});
            xticklabels(ax2,{});
            xticklabels(ax3,{});
            xticklabels(ax4,{});
            h.TileSpacing = 'compact';
            fontsize(scale=1.2);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
