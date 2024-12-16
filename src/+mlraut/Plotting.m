classdef Plotting < handle & mlsystem.IHandle
    %% line1
    %  line2
    %  
    %  Created 10-Feb-2024 13:59:48 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    

    properties
        do_plot_emd
        plot_range
    end

    properties (Dependent)
        anatomy_list
        rsn_list
        do_save
        Fs
        HCP_signals
        global_signal
        hp_thresh
        lp_thresh
        num_frames
        out_dir
        physio_signal
        sub
        subjects
        tags
        task
        tasks
        tr
    end

    methods %% GET
        function g = get.anatomy_list(this)
            g = this.ias_.anatomy_list;
        end
        function g = get.rsn_list(this)
            g = this.ias_.rsn_list;
        end
        function g = get.do_save(this)
            g = this.ias_.do_save;
        end
        function g = get.Fs(this)
            g = this.ias_.Fs;
        end
        function g = get.global_signal(this)
            g = this.ias_.global_signal;
        end
        function g = get.HCP_signals(this)
            g = this.ias_.HCP_signals;
        end
        function g = get.hp_thresh(this)
            g = this.ias_.hp_thresh;
        end
        function g = get.lp_thresh(this)
            g = this.ias_.lp_thresh;
        end
        function g = get.num_frames(this)
            g = this.ias_.num_frames;
        end
        function g = get.out_dir(this)
            g = this.ias_.out_dir;
        end
        function g = get.physio_signal(this)
            g = this.ias_.physio_signal;
        end
        function g = get.sub(this)
            g = this.ias_.current_subject;
        end
        function g = get.subjects(this)
            g = this.ias_.subjects;
        end
        function g = get.tags(this)
            g = this.ias_.tags;
        end
        function g = get.task(this)
            g = this.ias_.current_task;
        end
        function g = get.tasks(this)
            g = this.ias_.tasks;
        end
        function g = get.tr(this)
            g = this.ias_.tr;
        end
    end

    methods
        function plot_emd(this, opts)
            %  Args:
            %      this mlraut.Plotting           
            %      opts.measure function_handle = @angle
            %      opts.anatomy {mustBeTextScalar} = 'ctx'

            arguments
                this mlraut.Plotting    
                opts.freq_limits {mustBeNumeric} = []
                opts.measure function_handle = @real
                opts.anatomy {mustBeTextScalar} = 'ctx'
            end
            assert(contains(opts.anatomy, {'cbm', 'ctx', 'str', 'thal'}))
            signals = this.HCP_signals.(lower(opts.anatomy));
            if isreal(signals)
                return
            end
            if isempty(opts.freq_limits)
                Flim = [this.hp_thresh this.lp_thresh];
            else
                Flim = opts.freq_limits;
            end

            % emd
            if ~this.do_plot_emd % which automates many plots
                for k = 1:length(this.rsn_list)
                    emd(opts.measure(signals(:,k)), SiftMaxIterations=256)
                    set(gcf, Position=[0 0 2880*0.618 2880*0.618]);
                    title(sprintf('EMD %s, showing 3 IMFs, RSN %s\n', char(opts.measure), this.rsn_list{k}), FontSize=14)
                end
            end

            % hht
            h2 = figure;
            h2.Position = [0 0 2880 2880*0.618];
            tiledlayout(3,3);
            for k = 1:length(this.rsn_list)
                nexttile
                hht(emd(opts.measure(signals(:,k))), this.Fs, FrequencyLimits=Flim); % MaxNumIMF=5
                title(sprintf('Hilbert Spectrum, opts.measure, %s', char(opts.measure), this.rsn_list{k}))
            end

            % fsst
            h3 = figure;
            h3.Position = [0 0 2880 2880*0.618];
            tiledlayout(3,3);
            for k = 1:length(this.rsn_list)
                nexttile
                fsst(opts.measure(signals(:,k)), this.Fs, 'yaxis')
                title(sprintf('Fourier synchrosqueezed transform, %s, %s', char(opts.measure), this.rsn_list{k}))
            end

            this.saveFigures(sprintf('%s_%s', char(opts.measure), opts.anatomy));
        end
        function h1 = plot_global_physio(this, opts)
            %  Args:
            %      this mlraut.Plotting   
            %      opts.measure function_handle = @angle    

            arguments
                this mlraut.Plotting
                opts.measure function_handle = @abs
            end
            secs_ = this.tr * (0:this.num_frames-1);

            h1 = figure;
            h1.Position = [0 0 2880 2880*0.618];
            hold on
            plot(secs_, opts.measure(this.global_signal));
            plot(secs_, opts.measure(this.physio_signal), '--', LineWidth=2);
            legend({'global', 'arousal'}, FontSize=18)
            xlabel('time/s', FontSize=24)
            ylabel(char(opts.measure), FontSize=24, Interpreter="none")
            title(sprintf('Global Signal, Arousal Signal, sub-%s, %s ', this.sub, this.task), FontSize=24, Interpreter="none")
            hold off

            if this.do_save
                this.saveFigures(char(opts.measure));
            end
        end
        function plot_regions(this, funh, opts)
            arguments
                this mlraut.Plotting
                funh function_handle
                opts.measure function_handle
                opts.anatomy {mustBeTextScalar} = 'ctx'
            end

            for anat = this.anatomy_list
                opts.anatomy = anat{1};
                funh(measure=opts.measure, anatomy=opts.anatomy);
            end
        end
        function [h1,h3] = plot_networks(this, opts)
            %  Args:
            %      this mlraut.Plotting           
            %      opts.measure function_handle = @this.ias_.X
            %      opts.anatomy {mustBeTextScalar} = 'ctx'
            %      opts.plot_range {mustBeInteger} = 150:300

            arguments
                this mlraut.Plotting           
                opts.measure function_handle = @this.ias_.X
                opts.anatomy {mustBeTextScalar} = 'ctx'
                opts.plot_range {mustBeInteger} = []
            end
            assert(contains(opts.anatomy, {'cbm', 'ctx', 'str', 'thal'}))
            ksi = this.HCP_signals.(lower(opts.anatomy)).ksi;
            eta = this.HCP_signals.(lower(opts.anatomy)).eta;
            if isempty(opts.plot_range)
                opts.plot_range = this.plot_range;
            end
            secs_ = this.tr * (0:this.num_frames-1);

            % plot Yeo's 7 RSNs
            h1 = figure;
            h1.Position = [0 0 2880 2880*0.618];
            hold on
            % RSNs 1-5 are extrinsic
            for k = 1:5
                meas_ = real(opts.measure(ksi(:, k), eta(:, k)));
                plot(secs_(opts.plot_range), meas_(opts.plot_range));
            end
            % RSNs 6-7 are instrinsic
            meas6_ = real(opts.measure(ksi(:, 6), eta(:, 6)));
            plot(secs_(opts.plot_range), meas6_(opts.plot_range), '--', LineWidth=2); % frontoparietal
            meas7_ = real(opts.measure(ksi(:, 7), eta(:, 7)));
            plot(secs_(opts.plot_range), meas7_(opts.plot_range), '--', LineWidth=2); % default mode
            legend(this.rsn_list(1:7), FontSize=18)
            xlabel('time/s', FontSize=24)
            ylabel(sprintf('%s(%s_signals)', char(opts.measure), opts.anatomy), FontSize=24, Interpreter="none")
            title(sprintf('Yeo RSNs, sub-%s, %s ', this.sub, this.task), FontSize=24, Interpreter="none")
            hold off

            % plot task+, task- RSNs
            h3 = figure;
            h3.Position = [0 0 2880 2880*0.618];
            hold on
            % 8 is extrinsic (task+)
            meas8_ = real(opts.measure(ksi(:, 8), eta(:, 8)));
            plot(secs_(opts.plot_range), meas8_(opts.plot_range));
            % 9 is instrinsic (task-)
            meas9_ = real(opts.measure(ksi(:, 9), eta(:, 9)));
            plot(secs_(opts.plot_range), meas9_(opts.plot_range), '--', LineWidth=2);
            legend(this.rsn_list(8:9), FontSize=18)
            xlabel('time/s', FontSize=24)
            ylabel(sprintf('%s(%s_signals)', char(opts.measure), opts.anatomy), FontSize=24, Interpreter="none")
            title(sprintf('Task +/-, sub-%s, %s ', this.sub, this.task), FontSize=24, Interpreter="none")
            hold off

            this.saveFigures(sprintf('%s_%s_', char(opts.measure), opts.anatomy));
        end
        function [h,h1,h2] = plot_radar(this, opts)
            %  Args:
            %      this mlraut.Plotting           
            %      opts.anatomy {mustBeTextScalar} = 'ctx'

            arguments
                this mlraut.Plotting
                opts.measure function_handle = @this.identity
                opts.anatomy {mustBeText} = 'ctx'
            end
            assert(contains(opts.anatomy, {'cbm', 'ctx', 'str', 'thal'}))
            signals = this.HCP_signals.(lower(opts.anatomy));
            if isreal(signals)
                return
            end

            % plot "radar" of RSNs
            h = figure;
            h.Position = [0 0 2880*0.618 2880*0.618];
            hold on
            for k = 1:3
                plot(signals(:, k), '.', MarkerSize=4);
            end
            plot(signals(:, 6), '.', MarkerSize=8)
            plot(signals(:, 7), '.', MarkerSize=8)
            legend([this.rsn_list(1:3) this.rsn_list(6:7)], FontSize=18)
            xlabel(sprintf('real(%s_signals)', opts.anatomy), FontSize=24, Interpreter="none")
            ylabel(sprintf('imag(%s_signals)', opts.anatomy), FontSize=24, Interpreter="none")
            hold off

            % plot "radar" of task+ and task-
            h1 = figure;
            h1.Position = [0 0 2880*0.618 2880*0.618];
            hold on
            plot(signals(:, 8), '.', MarkerSize=8)
            plot(signals(:, 9), '.', MarkerSize=8)
            %plot(this.global_signal, '-.')
            %plot(this.physio_signal, '-.')
            legend({'task+', 'task-'}, FontSize=18)
            xlabel(sprintf('real(%s_signals)', opts.anatomy), FontSize=24, Interpreter="none")
            ylabel(sprintf('imag(%s_signals)', opts.anatomy), FontSize=24, Interpreter="none")
            hold off

            % plot "radar" of global signal, arousal signal
            h2 = figure;
            h2.Position = [0 0 2880*0.618 2880*0.618];
            hold on
            plot(this.global_signal, '.', MarkerSize=8)
            plot(this.physio_signal, '.', MarkerSize=8)
            legend({'global', 'arousal'}, FontSize=18)
            xlabel('real(signal)', FontSize=24, Interpreter="none")
            ylabel('imag(signal)', FontSize=24, Interpreter="none")
            hold off

            this.saveFigures(sprintf('%s_%s', char(opts.measure), opts.anatomy));
        end
        function [h,h1] = plot_timeseries_qc(this, tseries, opts)
            arguments
                this mlraut.Plotting
                tseries double {mustBeNonempty}
                opts.ylabel {mustBeTextScalar} = "time-series (arbitrary)"
                opts.Fs double = []
                opts.tr double = []
                opts.L double = []
            end
            if isempty(opts.Fs); opts.Fs = this.Fs; end
            if isempty(opts.tr); opts.tr = this.tr; end
            if isempty(opts.L); opts.L = this.num_frames; end
            times = ascol((0:opts.L-1)*opts.tr);
            tseries_centered = mean(tseries, 2) - mean(tseries, 'all');

            % plot times -> tseries
            h = figure;
            if isreal(tseries_centered)
                plot(times, tseries_centered);
            else
                plot(times, real(tseries_centered), times, imag(tseries_centered))
                legend(["real", "imag"])
            end
            xlabel("times (s)");
            ylabel(opts.ylabel);
            title(sprintf("%s: %s: %s", stackstr(3), stackstr(2), opts.ylabel), Interpreter="none");

            % plot times -> Fourier(tseries)
            h1 = figure;
            P2 = abs(fft(tseries_centered))/opts.L;
            P1 = P2(1:opts.L/2+1);
            P1(2:end-1) = 2*P1(2:end-1);
            logP1 = log(P1);
            freqs = opts.Fs*(0:(opts.L/2))/opts.L;
            plot(freqs(2:end), logP1(2:end));            
            xlabel("frequency (Hz)");
            ylabel("log "+opts.ylabel);
            title(sprintf("%s: %s: log %s", stackstr(3), stackstr(2), opts.ylabel), Interpreter="none");
        end
        function saveFigures(this, label, varargin)
            label = strrep(label, "@(varargin)", "");
            label = strrep(label, "(varargin{/})", "");
            saveFigures(this.out_dir, ...
                closeFigure=true, ...
                prefix=sprintf('%s_%s%s_%s_%s', stackstr(3, use_dashes=true), label, this.tags, this.sub, this.task));
        end

        function this = Plotting(ias, opts)
            %%
            %  Args:
            %      opts.do_plot_emd logical = false
            %      opts.plot_range {mustBeInteger} = 1:572 (1:158 for GBM)

            arguments
                ias mlraut.AnalyticSignal {mustBeNonempty}
                opts.do_plot_emd logical = false
                opts.plot_range {mustBeInteger} = 1:572
            end

            this.ias_ = ias;
            this.do_plot_emd = opts.do_plot_emd;
            this.plot_range = opts.plot_range;
        end
    end

    %% PROTECTED

    properties (Access = protected)
        ias_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
