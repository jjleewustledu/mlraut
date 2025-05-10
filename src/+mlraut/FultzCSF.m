classdef FultzCSF < handle
    %% Examine CSF, CE tumor, control samples to be comparable to
    %  Fultz, et al., Science 2019.  https://www.science.org/doi/10.1126/science.aax5440
    %  
    %  Created 02-May-2025 23:14:47 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 24.2.0.2923080 (R2024b) Update 6 for MACA64.  Copyright 2025 John J. Lee.
    
    properties
        analytic_signal_obj
        time_series_window = 200  % sec
        waves_dir = fullfile( ...
            getenv("HOME"), "MATLAB-Drive", "arousal-waves-main")
    end

    properties (Dependent)
        cbm_mask
        ctx_mask
        Fs
        hcp_signals
        phi  % may be iFV-brightest, gray, CE_on_T1w
        psi
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
            g = this.analytic_signal_obj.Fs;
        end

        function g = get.hcp_signals(this)
            g = this.analytic_signal_obj.HCP_signals;
        end

        function g = get.phi(this)
            g = this.analytic_signal_obj.physio_signal(:, 1);
        end

        function g = get.psi(this)
            g = this.analytic_signal_obj.bold_signal;
        end

        function g = get.source_physio(this)
            g = this.analytic_signal_obj.source_physio;
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
            g = this.analytic_signal_obj.tr;
        end
    end

    methods
        function this = FultzCSF(aso)
            arguments
                aso = []
            end

            this.analytic_signal_obj = aso;
        end

        function psi = masked_psi(this, region)
            psi = mean(this.psi(:, this.(region + "_mask")), 2, 'omitnan');
        end

        function h = plot_coherencyc(this, opts)
            arguments
                this mlraut.FultzCSF
                opts.tseries {mustBeTextScalar} = "-dbold/dt"
            end

            regions = ["ctx", "cbm", "str", "thal"];
            lenr = length(regions);
            tseries_phi = repmat(this.phi, [1, lenr]);
            tseries_psi = zeros(size(tseries_phi, 1), lenr);
            switch char(opts.tseries)
                case 'bold'
                    tseries_phi = real(tseries_phi);
                    for ir = 1:lenr
                        tseries_psi(:,ir) = real(this.masked_psi(regions(ir)));
                    end
                case '-dbold/dt'
                    tseries_phi = real(tseries_phi(1:end-1, :));
                    tseries_psi = zeros(size(tseries_phi, 1), lenr);
                    for ir = 1:lenr
                        bold = real(this.masked_psi(regions(ir)));
                        tseries_psi(:, ir) = -diff(bold);
                    end
                    tseries_psi(tseries_psi < 0) = 0;
                case 'X'
                    for ir = 1:lenr
                        psi_ = this.masked_psi(regions(ir));
                        tseries_psi(:, ir) = this.analytic_signal_obj.X(psi_, this.phi);
                    end
                case 'Y'
                    for ir = 1:lenr
                        psi_ = this.masked_psi(regions(ir));
                        tseries_psi(:, ir) = this.analytic_signal_obj.Y(psi_, this.phi);
                    end
                case 'Z'
                    for ir = 1:lenr
                        psi_ = this.masked_psi(regions(ir));
                        tseries_psi(:, ir) = this.analytic_signal_obj.Z(psi_, this.phi);
                    end
            end

            T = size(tseries_phi, 1) * this.tr;
            W = 0.1;  % bandwidth
            Ntapers = 30;
            p = floor(2 * T * W - Ntapers);
            pval = 0.05;

            params.Fs = this.Fs;
            params.fpass = [0, W];
            params.trialave = 1;  % T/F
            params.tapers = [W, T, p];
            params.err = [2, pval];

            [coh,pha,S12,S1,S2,f,conf_coh,pha_std,coh_err] = coherencyc( ...
                tseries_psi, tseries_phi, params);
            f = ascol(f);
            [~,f_0p05_idx] = min(abs(f - 0.05));
            pha = unwrap(pha);
            
            % https://www.mathworks.com/help/matlab/creating_plots/combine-multiple-plots.html            
            h = tiledlayout(5, 1);

            ax1 = nexttile;
            plot(ax1, f, coh);
            ylabel(ax1, "coherence")
            legend(ax1, regions)

            ax2 = nexttile;
            plot(ax2, f, pha);
            ylabel(ax2, "phase")

            ax3 = nexttile;
            plot(ax3, f, 10 * log10(S12/S12(f_0p05_idx)));
            ylabel(ax3, "power (dB, 0.05 Hz ref.)")
            title(ax3, "cross-spectral")

            ax4 = nexttile;
            plot(ax4, f, 10 * log10(S1/S1(f_0p05_idx)));
            ylabel(ax4, "power (dB, 0.05 Hz ref.)")
            title(ax4, opts.tseries)

            ax5 = nexttile;
            S2 = S2(:, 1);
            plot(ax5, f, 10 * log10(S2/S2(f_0p05_idx)), LineWidth=1.25, Color="m");
            ylabel(ax5, "power (dB, 0.05 Hz ref.)")
            title(ax5, this.source_physio, Interpreter="none")

            linkaxes([ax1,ax2,ax3,ax4,ax5], 'x');
            xlabel(h, "frequency (Hz)");
            xticklabels(ax1,{});
            xticklabels(ax2,{});
            xticklabels(ax3,{});
            xticklabels(ax4,{});
            h.TileSpacing = 'compact';
            fontsize(scale=1.2);
        end

        function h = plot_xcorr(this)
            h = [];
        end

        function h = plot_mtspectrumc(this, opts)
            arguments
                this mlraut.FultzCSF
                opts.tseries {mustBeTextScalar} = "phi"
                opts.region {mustBeTextScalar} = "ctx"
            end

            switch char(opts.tseries)
                case {'phi', 'csf'}
                    tseries = real(this.phi);
                case 'bold'
                    tseries = real(this.masked_psi(opts.region));
                case '-dbold/dt'
                    bold = real(this.masked_psi(opts.region));
                    tseries = -diff(bold);
                    tseries(tseries < 0) = 0;
                case 'X'
                    psi_ = this.masked_psi(opts.region);
                    tseries = real(this.analytic_signal_obj.X(psi_, this.phi));
                case 'Y'
                    psi_ = this.masked_psi(opts.region);
                    tseries = real(this.analytic_signal_obj.Y(psi_, this.phi));
                case 'Z'
                    psi_ = this.masked_psi(opts.region);
                    tseries = real(this.analytic_signal_obj.Z(psi_, this.phi));
            end

            T = numel(tseries) * this.tr;
            W = 0.1;  % bandwidth
            p = 2 * T * W - 5;

            params.Fs = this.Fs;
            % params.fpass = [0, W];
            params.trialave = 1;  % true
            params.tapers = [W, T, p];
            params.err = [2, p];

            tseries = asrow(tseries);
            [S,f] = mtspectrumc(zscore(tseries), params);
            [~,f_0p05_idx] = min(abs(f - 0.05));
            h = plot(asrow(f), asrow(10 * log10(S/S(f_0p05_idx))));
            xlabel("frequency (Hz)")
            ylabel("power (dB against 0.05 Hz reference)")
            title(opts.tseries)
        end

        function h = plot_csf(this)
            if ~strcmp(this.source_physio, "iFV-brightest")
                h = [];
                return
            end

            Nt = round(this.time_series_window / this.tr);
            t1 = Nt * this.tr;
            t2 = (2 * Nt - 1) * this.tr;
            t = t1:this.tr:t2;

            h = plot(t, zscore(real(this.phi(Nt:2*Nt-1))));
            xlabel("time (s)");
            ylabel("BOLD signal (z)");
            title("Exemplar BOLD time-series for CSF")
            fontsize(scale=1.618);
        end

        function h = plot_bold(this)

            Nt = round(this.time_series_window / this.tr);
            t1 = Nt * this.tr;
            t2 = (2 * Nt - 1) * this.tr;
            t = t1:this.tr:t2;

            boldas = this.hcp_signals.ctx.psi ./ this.hcp_signals.ctx.phi;
            bold = asrow(mean(boldas(Nt:2*Nt-1, 1:7), 2));
            h = plot(t, zscore(real(bold)));
            xlabel("time (s)");
            ylabel("BOLD signal (z)");
            title("Exemplar BOLD time-series for cortex")
            % fontsize(scale=1.618);
        end

        function h = plot_neg_dbold_dt(this)

            Nt = round(this.time_series_window / this.tr);
            t1 = Nt * this.tr;
            t2 = (2 * Nt - 1) * this.tr;
            t = t1:this.tr:t2;

            boldas = this.hcp_signals.ctx.psi ./ this.hcp_signals.ctx.phi;
            bold = asrow(mean(boldas(Nt:2*Nt, 1:7), 2));
            neg_dbdt = -diff(real(bold));
            neg_dbdt(neg_dbdt < 0) = 0;
            h = plot(t, zscore(neg_dbdt));
            xlabel("time (s)");
            ylabel("-d/dt BOLD signal (z)");
            % fontsize(scale=1.618);
        end
    end

    properties (Access = private)
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
