classdef PlottingWavelets < handle & mlraut.Plotting
    %% line1
    %  line2
    %  
    %  Created 23-Dec-2024 21:13:42 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 24.2.0.2806996 (R2024b) Update 3 for MACA64.  Copyright 2024 John J. Lee.
    
    methods
        function this = PlottingWavelets(varargin)
            this = this@mlraut.Plotting(varargin{:});
        end

        function h1 = plot_cmor(this, opts)
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
                opts.anatomy {mustBeTextScalar} = 'ctx'
            end

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
            % RSNs 1-6 are extrinsic; 7 is DMN
            for k = 1:7
                meas_ = real(opts.measure(psi(:, k), phi(:, k)));
                
                try
                    % Perform continuous wavelet transform
                    % Using complex Morlet wavelet (cmor) for good frequency localization
                    [cwt_coefficients,frequencies] = cwt(meas_, this.Fs);

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
            xlim([this.ias_.hp_thresh, this.ias_.lp_thresh]);

            legend(this.rsn_list(1:7))
            xlabel('Frequency (Hz)')
            ylabel(sprintf('Wavelet power (%s, %s)', meas_label, this.expand_abbr(opts.anatomy)), Interpreter="none")
            title(sprintf('Yeo RSNs, %s, %s ', this.sub, this.source_physio), Interpreter="none")
            fontsize(scale=this.fontscale)
            hold off

            this.saveFigures(sprintf('%s_%s_', char(meas_label), opts.anatomy));
        end

        function h1 = plot_cwt(this, opts)
            arguments
                this mlraut.Plotting
                opts.measure = []  % unused, maintaining interface
                opts.anatomy {mustBeTextScalar} = 'ctx'
            end

            psi = this.HCP_signals.(lower(opts.anatomy)).psi;
            phi = this.HCP_signals.(lower(opts.anatomy)).phi;
            zeta = psi ./ phi;

            % RSNs 1-6 are extrinsic; 7 is DMN
            for k = 1:7
                try
                    h1 = figure;
                    cwt(zeta(:, k), this.Fs);
                    fontsize(scale=1.5)
                    sgtitle(sprintf('Magnitude scalogram (%s, %s)', this.expand_abbr(opts.anatomy), this.rsn_list{k}));
                catch ME
                    handwarning(ME)
                end
            end

            this.saveFigures(sprintf('%s_%s_', 'cwt', opts.anatomy));
        end

        function h1 = plot_wcoherence(this, opts)
            arguments
                this mlraut.Plotting
                opts.measure = []  % unused, maintaining interface
                opts.anatomy {mustBeTextScalar} = 'ctx'
            end

            psi = this.HCP_signals.(lower(opts.anatomy)).psi;
            phi = this.HCP_signals.(lower(opts.anatomy)).phi;

            % RSNs 1-6 are extrinsic; 7 is DMN
            for k = 1:7
                try
                    h1 = figure;
                    x = real(phi(:, k));
                    y = real(psi(:, k));
                    wcoherence(x, y);
                    title(sprintf('Wavelet coherence (%s, %s)', this.expand_abbr(opts.anatomy), this.rsn_list{k}));
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
