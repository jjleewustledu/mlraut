classdef PlottingGBM < handle & mlraut.PlottingWavelets
    %% line1
    %  line2
    %  
    %  Created 23-Dec-2024 01:01:46 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 24.2.0.2806996 (R2024b) Update 3 for MACA64.  Copyright 2024 John J. Lee.
    
    properties (Dependent)
        gbm_list
    end

    methods %% GET
        function g = get.gbm_list(this)
            g = this.ias_.gbm_list;
        end
    end

    methods
        function this = PlottingGBM(varargin)
            this = this@mlraut.PlottingWavelets(varargin{:});
        end

        function plot_regions_gbm(this, funh, opts)
            arguments
                this mlraut.PlottingGBM
                funh function_handle
                opts.measure function_handle
            end

            % plot_regions@mlraut.PlottingWavelets(this, funh, measure=opts.measure);

            funh(measure=opts.measure, anatomy="gbm");
        end

        function [h1,h2] = plot_networks_gbm(this, opts)
            %  Args:
            %      this mlraut.Plotting           
            %      opts.measure function_handle = @this.ias_.X
            %      opts.anatomy {mustBeTextScalar} = 'ctx'
            %      opts.plot_range {mustBeInteger} = 150:300

            arguments
                this mlraut.PlottingGBM           
                opts.measure function_handle = @this.ias_.X
                opts.anatomy {mustBeTextScalar} = "gbm"
                opts.plot_range {mustBeInteger} = []
            end

            % plot_networks@mlraut.Plotting(this, measure=opts.measure, anatomy = opts.anatomy, plot_range=opts.plot_range);

            psi = this.HCP_signals.(lower(opts.anatomy)).psi;
            phi = this.HCP_signals.(lower(opts.anatomy)).phi;
            if isempty(opts.plot_range)
                opts.plot_range = this.plot_range;
            end
            secs_ = ascol(this.tr*((0:this.num_frames-1) + this.num_frames_to_trim));

            meas_label = char(opts.measure);
            meas_label = strrep(meas_label, "@(varargin)", "");
            meas_label = strrep(meas_label, "(varargin{:})", "");
            meas_label = strrep(meas_label, "this.", "");

            % get nice colours from colorbrewer
            % (https://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab)
            cb = cbrewer2('Set1', 3);
            cb = [cb, [1; 0.75; 0.5]];  % 80% brightness, 70% alpha for task+

            %% plot GBM regions

            h1 = figure;
            h1.Position = [0 0 2880 2880*0.618];
            hold on
            for k = 1:length(this.gbm_list)
                meas_ = real(opts.measure(psi(:, k), phi(:, k)));
                p = plot(secs_(opts.plot_range), meas_(opts.plot_range), LineWidth=7);
                p.Color=cb(k,:);
            end
            legend(this.gbm_list)
            xlabel('time/s')
            ylabel(sprintf('%s(%s signals)', meas_label, 'tumor'), Interpreter="none")
            title(sprintf('Tumor regions, %s, %s ', this.sub, this.source_physio), Interpreter="none")
            fontsize(scale=this.fontscale)
            hold off

            this.saveFigures(sprintf('%s_%s_', char(meas_label), opts.anatomy));
        end

        function [h1,h2] = plot_networks_dots_gbm(this, opts)
            %  Args:
            %      this mlraut.Plotting           
            %      opts.measure function_handle = @this.ias_.X
            %      opts.anatomy {mustBeTextScalar} = 'ctx'
            %      opts.plot_range {mustBeInteger} = 150:300

            arguments
                this mlraut.PlottingGBM           
                opts.measure function_handle = @this.ias_.angle
                opts.anatomy {mustBeTextScalar} = "gbm"
                opts.marker_size double = 100
                opts.plot_range {mustBeInteger} = []
            end

            % plot_networks_dots@mlraut.Plotting(this, measure=opts.measure, anatomy = opts.anatomy, marker_size=opts.marker_size, plot_range=opts.plot_range);

            ms = opts.marker_size;
            psi = this.HCP_signals.(lower(opts.anatomy)).psi;
            phi = this.HCP_signals.(lower(opts.anatomy)).phi;
            if isempty(opts.plot_range)
                opts.plot_range = this.plot_range;
            end
            secs_ = ascol(this.tr*((0:this.num_frames-1) + this.num_frames_to_trim));

            meas_label = char(opts.measure);
            meas_label = strrep(meas_label, "@(varargin)", "");
            meas_label = strrep(meas_label, "(varargin{:})", "");
            meas_label = strrep(meas_label, "this.", "");

            % get nice colours from colorbrewer
            % (https://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab)
            cb = cbrewer2('Set1', 3);

            %% plot GBM regions

            h1 = figure;
            h1.Position = [0 0 2880 2880*0.618];
            hold on
            for k = 1:length(this.gbm_list)
                meas_ = real(opts.measure(psi(:, k), phi(:, k)));
                scatter(secs_(opts.plot_range), meas_(opts.plot_range), ms, cb(k,:), 'filled', MarkerFaceAlpha=0.5, MarkerEdgeAlpha=0.5);
            end
            legend(this.gbm_list)
            xlabel('time/s')
            ylabel(sprintf('%s(%s signals)', meas_label, 'tumor'), Interpreter="none")
            title(sprintf('Tumor regions, %s, %s ', this.sub, this.source_physio), Interpreter="none")
            fontsize(scale=this.fontscale)
            hold off

            this.saveFigures(sprintf('%s_%s_', char(meas_label), opts.anatomy));
        end

        function h1 = plot_cmor_gbm(this, opts)
            %% Claude:  
            %  I need to detect features in time series data that have characteristic waveforms, but the starting times
            %  for emergence of features is unknown.  Are spectral methods best for identifying these features in new
            %  data?  How do traditional power spectra compare to wavelets or other spectral methods for identifying
            %  features identifiable as waveforms? 
            %
            %  When LIGO detected the first observed gravitational waves, what kind of  time-series analyses were used
            %  successfully?  Is the chirp of gravitational waves easy to find with continuous wavelet transformation?
            %
            %  In matlab, how should I look for power-law behavior in time-series using matlab's wavelet toolbox?
            %
            %  Args:
            %      this mlraut.Plotting           
            %      opts.measure function_handle = @this.ias_.X
            %      opts.anatomy {mustBeTextScalar} = 'ctx'
            %      opts.plot_range {mustBeInteger} = 150:300

            arguments
                this mlraut.Plotting           
                opts.measure function_handle = @this.ias_.X
                opts.anatomy {mustBeTextScalar} = 'gbm'
            end

            % plot_cmor@mlraut.PlottingWavelets(this, measure=opts.measure, anatomy = opts.anatomy);

            psi = this.HCP_signals.(lower(opts.anatomy)).psi;
            phi = this.HCP_signals.(lower(opts.anatomy)).phi;

            meas_label = char(opts.measure);
            meas_label = strrep(meas_label, "@(varargin)", "");
            meas_label = strrep(meas_label, "(varargin{:})", "");
            meas_label = strrep(meas_label, "this.", "");

            % get nice colours from colorbrewer
            % (https://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab)
            cb = flip(cbrewer2('Spectral', 7), 1);  % task+ is blue, task- is red
            cb = [0.8*cb, [0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.9]];  % 80% brightness, 70% alpha for task+

            %% plot Yeo's 7 RSNs

            h1 = figure;
            h1.Position = [0 0 1600 1600*0.618];
            hold on
            % CE, edema, WT
            for k = 1:length(this.gbm_list)
                meas_ = real(opts.measure(psi(:, k), phi(:, k)));
                
                try
                    % Perform continuous wavelet transform
                    % Using complex Morlet wavelet (cmor) for good frequency localization
                    [cwt_coefficients,frequencies]  = cwt(meas_, this.Fs);

                    % Calculate wavelet power spectrum
                    power_spectrum = abs(cwt_coefficients).^2;

                    % Calculate average power at each scale
                    mean_power = mean(power_spectrum, 2);

                    % Plot in log-log scale to identify power-law behavior
                    if k < 7
                        p = loglog(frequencies, mean_power, LineWidth=5);
                    else
                        p = loglog(frequencies, mean_power, LineWidth=7);
                    end
                    p.Color=cb(k,:);
                    grid on;

                    %disp(this)

                    % Open wavelet analysis tool
                    %wvtool(meas_, fs);
                catch ME
                    handwarning(ME)
                end
            end
            xlim([0, this.ias_.lp_thresh]);

            legend(this.gbm_list)
            xlabel('Frequency (Hz)')
            ylabel(sprintf('Wavelet power (%s, %s)', meas_label, this.expand_abbr(opts.anatomy)), Interpreter="none")
            title(sprintf('Yeo RSNs, %s, %s ', this.sub, this.source_physio), Interpreter="none")
            fontsize(scale=this.fontscale)
            hold off

            this.saveFigures(sprintf('%s_%s_', char(meas_label), opts.anatomy));
        end

        function h1 = plot_cwt_gbm(this, opts)
            arguments
                this mlraut.Plotting
                opts.measure = []  % unused, maintaining interface
                opts.anatomy {mustBeTextScalar} = 'gbm'
            end

            % plot_cwt@mlraut.PlottingWavelets(this, measure=opts.measure, anatomy = opts.anatomy);

            psi = this.HCP_signals.(lower(opts.anatomy)).psi;
            phi = this.HCP_signals.(lower(opts.anatomy)).phi;
            phi = ones(size(psi)).*mean(phi, 2);
            zeta = psi ./ phi;

            % CE, edema, WT
            for k = 1:length(this.gbm_list)
                try
                    h1 = figure;
                    cwt(zeta(:, k), this.Fs);
                    fontsize(scale=1.5)
                    sgtitle(sprintf('Magnitude scalogram (%s, %s)', this.expand_abbr(opts.anatomy), this.gbm_list(k)))
                catch ME
                    handwarning(ME)
                end
            end

            this.saveFigures(sprintf('%s_%s_', 'cwt', opts.anatomy));
        end

        function h1 = plot_wcoherence_gbm(this, opts)
            arguments
                this mlraut.Plotting
                opts.measure = []  % unused, maintaining interface
                opts.anatomy {mustBeTextScalar} = 'gbm'
            end

            % plot_wcoherence@mlraut.PlottingWavelets(this, measure=opts.measure, anatomy = opts.anatomy);

            psi = this.HCP_signals.(lower(opts.anatomy)).psi;
            phi = this.HCP_signals.(lower(opts.anatomy)).phi;

            % CE, edema, WT
            for k = 1:length(this.gbm_list)
                try
                    h1 = figure;
                    x = real(phi(:, k));
                    y = real(psi(:, k));
                    wcoherence(x, y);
                    title(sprintf('Wavelet coherence (%s, %s)', this.expand_abbr(opts.anatomy), this.gbm_list(k)));
                    fontsize(scale=1.5)
                catch ME
                    handwarning(ME)
                end
            end

            this.saveFigures(sprintf('%s_%s_', 'wcoherence', opts.anatomy));
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
