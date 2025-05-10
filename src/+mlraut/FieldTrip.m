classdef FieldTrip
    %% line1
    %  line2
    %  
    %  Created 17-Sep-2024 11:13:11 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 24.1.0.2689473 (R2024a) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    

    properties
        Fs = 508.6275  % from *_rmegpreproc.mat:  data.fsample
        hp_thresh = 0.01
        lp_thresh = 0.1
        sub
    end

    properties (Dependent)
        datadir
        fieldtripdir
        headmodel_file
        indices_inside_iFV_file
        inside_3p75mm_ifc
        neighbours_file
        rmegpreproc_files  % string array
        sourcemodel2d_file
        sourcemodel3d_file
        megdir
    end

    methods %% GET
        function g = get.datadir(this)
            if contains(hostname, "twistor")
                g = fullfile("/Volumes", "PrecunealSSD2", "HCP", "AWS", "hcp-openaccess", "HCP_1200");
                return
            end
            if contains(hostname, "MATTA")
                g = fullfile("D:", "HCP_1200");
                return
            end
        end
        function g = get.fieldtripdir(this)
            g = fullfile(getenv("HOME"), "MATLAB-Drive", "fieldtrip-20240731");
        end
        function g = get.headmodel_file(this)
            g = fullfile(this.datadir, this.sub, "MEG", "anatomy", this.sub + "_MEG_anatomy_headmodel.mat");
        end
        function g = get.indices_inside_iFV_file(this)
            g = fullfile(this.datadir, this.sub, "MEG", "indices_inside_iFV.mat");
        end
        function g = get.inside_3p75mm_ifc(this)
            g = mlfourd.ImagingFormatContext2( ...
                fullfile(this.megdir, "inside_3p75mm.nii.gz"));
        end
        function g = get.neighbours_file(this)
            g = fullfile(this.fieldtripdir, "template", "neighbours", "bti248_neighb.mat");
        end
        function g = get.rmegpreproc_files(this)
            g = fullfile(this.datadir, this.sub, "MEG", "Restin", "rmegpreproc", this.sub + "_MEG_*-Restin_rmegpreproc.mat");
            g = asrow(mglob(g));
        end
        function g = get.megdir(this)
            g = fullfile(this.datadir, this.sub, "MEG");
        end
        function g = get.sourcemodel2d_file(this)
            g = fullfile(this.datadir, this.sub, "MEG", "anatomy", this.sub + "_MEG_anatomy_sourcemodel_2d.mat");
        end
        function g = get.sourcemodel3d_file(this)
            g = fullfile(this.datadir, this.sub, "MEG", "anatomy", this.sub + "_MEG_anatomy_sourcemodel_3d4mm.mat");
        end
    end

    methods
        function this = FieldTrip(opts)
            %% FIELDTRIP 
            %  Args:
            %      arg1 (its_class): Description of arg1.
            
            arguments
                opts.sub {mustBeTextScalar} = "175237"
            end

            addpath(fullfile(getenv("HOME"), "MATLAB-Drive", "fieldtrip-20240731"))
            ft_defaults;
            
            this.sub = opts.sub;
        end

        function dat1 = build_band_passed(this, dat)
            %% Implements butter:  web(fullfile(docroot, 'signal/ref/butter.html?browser=F1help#bucsfmj')) .
            %  See also web(fullfile(docroot, 'signal/ug/practical-introduction-to-digital-filtering.html')) .
            %  Returns:
            %      dat1 same num. type as dat

            if all(dat == 0)
                dat1 = dat;
                return
            end
            if isempty(this.lp_thresh) && isempty(this.hp_thresh)
                dat1 = dat;
                return
            end
            [z,p,k] = butter(2, [this.hp_thresh, this.lp_thresh - eps('single')]/(this.Fs/2)); % digital Wn in [0, 1]
            [sos,g] = zp2sos(z, p, k);
            dat1 = filtfilt(sos, g, double(dat));
            if isa(dat, 'single')
                dat1 = single(dat1);
            end
            if isa(dat, 'double')
                dat1 = double(dat1);
            end
        end

        function psi = build_centered(~, psi)
            assert(~isempty(psi))
            if all(psi == 0)
                return
            end
            psi = psi - median(psi, 'all', 'omitnan');
        end

        function psi = build_centered_and_rescaled(this, psi, pos)
            %% Mimics z-score of |psi(t,x)> using median and mad.
            %  To avoid round-off errors, rescale before centering.

            psi = this.build_rescaled(psi, pos);
            psi = this.build_centered(psi);
        end

        function psi = build_rescaled(~, psi, pos)
            %% MEG timeseries drop off significantly in centrally located sources; "central bias",
            %  so rescale concentric 2-sphere of sources separately.
            %  To ignore location biases, use pos = [].

            assert(~isempty(psi))
            assert(isempty(pos) || size(psi, 2) == size(pos, 1))  % one pos for each source in psi
            if all(psi == 0)
                return
            end

            if isempty(pos)
                d = mad(abs(psi), 1, 'all');  % median abs. dev.
                psi = psi./d;
                return
            end

            % index 2-spheres; N.B.:  z \in (-5.51, 8.66) cm
            layers = 0:ceil(max(pos, [], "all"));
            for ilayer = 2:length(layers)
                radius = sqrt(pos(:,1).^2 + pos(:,2).^2 + pos(:,3).^2);
                select = layers(ilayer-1) < radius & radius <= layers(ilayer);  % selected 2-sphere
                if sum(select) > 0
                    psi_selected = psi(:, select);
                    d = mad(abs(psi_selected), 1, 'all');  % median abs. dev.
                    psi_selected = psi_selected./d;
                    psi(:, select) = psi_selected;
                end
            end     
        end

        function source = call(this, source, opts)
            %% CALL
            %
            % >> source
            %
            %   struct with fields:
            %
            %         time: [0 0.0020 0.0039 0.0059 0.0079 0.0098 0.0118 0.0137 0.0157 0.0177 0.0196 0.0216 0.0236 0.0255 0.0275 0.0295 0.0314 0.0334 0.0354 0.0373 0.0393 0.0412 0.0432 0.0452 … ] (1×132340 double)
            %          cfg: [1×1 struct]
            %     coordsys: 'bti'
            %          dim: [38 48 41]
            %       inside: [74784×1 logical]
            %          pos: [74784×3 double]
            %         unit: 'cm'
            %       method: 'average'
            %          avg: [1×1 struct]
            %
            % >> source.avg
            %
            %   struct with fields:
            %
            %          mom: {74784×1 cell}
            %          pow: [74784×132340 double]
            %     noisecov: {1×74784 cell}
            arguments
                this mlraut.FieldTrip
                source struct {mustBeNonempty}
                opts.use_cache logical = false
            end

            signal = source.avg.pow';  % N_times x N_pos
            signal = sqrt(signal);  % but signal ~ 20 log10(amplitude)?

            %% use N_times x N_positions

            N_1min = 509*60;  % sampling rate * sec
            N_times = N_1min; % N_1min or length(source.time);
            % N_pos = sum(source.inside);
            % pos_inside = source.pos(source.inside, :);  % z \in (-5.51, 8.66) cm

            field_fqfn = fullfile(this.megdir, "field_"+stackstr()+".mat");
            if opts.use_cache && isfile(field_fqfn)
                ld1 = load(field_fqfn);
                field = ld1.field;
                clear("ld1");
            else
                field = hilbert( ...
                    this.build_centered_and_rescaled( ...
                    signal(:, source.inside), []));
                save(field_fqfn, "field");
            end

            % physio with independent normalizations
            ld_iFV = load(this.indices_inside_iFV_file);  % 32 x 1
            physio = hilbert( ...
                mean( ...
                this.build_centered_and_rescaled( ...
                signal(:, ld_iFV.indices_inside_iFV, :), []), 2));
            save(fullfile(this.megdir, "physio_"+stackstr()+".mat"), "physio");

            analytic_signal_fqfn = fullfile(this.megdir, "analytic_signal_"+stackstr()+".mat");
            if opts.use_cache && isfile(analytic_signal_fqfn)
                ld2 = load(analytic_signal_fqfn);
                analytic_signal = ld2.analytic_signal;
                clear("ld2");
            else
                analytic_signal = conj(physio) .* field;
                save(analytic_signal_fqfn, "analytic_signal");
            end

            % save abs_as
            ifc = copy(this.inside_3p75mm_ifc);
            mat = zeros(N_times, prod(source.dim));
            mat(:, source.inside') = abs(analytic_signal(1:N_times, :));
            ifc.img = reshape(mat', [source.dim N_times]);
            ifc.filepath = this.megdir;
            ifc.fileprefix = "abs_as_frame1to"+N_times+"_"+stackstr();
            ifc.save();

            % save angle_as
            ifc = copy(this.inside_3p75mm_ifc);
            mat = zeros(N_times, prod(source.dim));
            mat(:, source.inside') = angle(analytic_signal(1:N_times, :));
            ifc.img = reshape(mat', [source.dim N_times]);
            ifc.filepath = this.megdir;
            ifc.fileprefix = "angle_as_frame1to"+N_times+"_"+stackstr();
            ifc.save();

            % save abs_as_avgt
            ifc = copy(this.inside_3p75mm_ifc);
            mat = zeros(1, prod(source.dim));
            mat(:, source.inside') = mean(abs(analytic_signal(1:N_times, :)), 1);  % time average
            ifc.img = reshape(mat', [source.dim]);
            ifc.filepath = this.megdir;
            ifc.fileprefix = "abs_as_avgt_"+stackstr();
            ifc.save();

            % save angle_as_avgt
            ifc = copy(this.inside_3p75mm_ifc);
            mat = zeros(1, prod(source.dim));
            mat(:, source.inside') = mean(angle(analytic_signal(1:N_times, :)), 1);  % time average
            ifc.img = reshape(mat', [source.dim]);
            ifc.filepath = this.megdir;
            ifc.fileprefix = "angle_as_avgt_"+stackstr();
            ifc.save();
        end

        function cdata = cell_rmegpreproc_data(this, opts)
            %% Resting MEG timeseries data in canonical FieldTrip structure, following HCP sensor censoring, denoising
            %  as described in https://www.humanconnectome.org/storage/app/media/documentation/meg1/MEG1_Release_Reference_Manual.pdf 
            %
            %  Args:
            %      this mlraut.FieldTrip
            %      opts.select_cell double = 1  % empty or < 1 collapses all cells into singleton
            % 
            %  data time ~ 1-cell containing double of size 1 x (130*1018); units ~ sec
            %  data.trial ~ 1-cell containing double of size ~ 245 x (130*1018) ~ Nsensors x Ntimes; units of
            %  magnetization density?
            % 
            %  Returns cell array, len ~ 3

            arguments
                this mlraut.FieldTrip
                opts.select_cell double = []  % empty or < 1 collapses all cells into singleton
            end
            
            cdata = cell(1, length(this.rmegpreproc_files));
            ci = 0;
            for globbed = this.rmegpreproc_files
                ci = ci + 1;
                ld = load(globbed);

                % select trialinfo, trial, time
                ld.data.trialinfo = ld.data.trialinfo(1);

                % select data.trial, data.time
                if isempty(opts.select_cell) || opts.select_cell < 1
                    ld.data.trial = {cell2mat(ld.data.trial)};
                    for ti = 2:length(ld.data.time)  % len ~ 130
                        T = ld.data.time{ti-1}(end);
                        ld.data.time{ti} = T + ld.data.time{ti};
                    end
                    ld.data.time = {cell2mat(ld.data.time)};
                else                    
                    ld.data.trial = ld.data.trial(opts.select_cell);
                    ld.data.time = ld.data.time(opts.select_cell);
                end

                cdata{ci} = ld.data;
            end
        end

        function data = selected_rmegpreproc_data(this, opts)

            arguments
                this mlraut.FieldTrip
                opts.run double = 1  % 1:3 for resting MEG, usually
            end

            cdata = this.cell_rmegpreproc_data();
            data = cdata{opts.run};
        end

        function source = sourcemodel(this, opts)
            arguments
                this mlraut.FieldTrip
                opts.kind {mustBeTextScalar} = "sourcemodel_3d4mm"
            end

            if contains(opts.kind, "2d")
                ld = load(this.sourcemodel2d_file);
                source = ld.sourcemodel2d;
                return
            end
            if contains(opts.kind, "3d")
                ld = load(this.sourcemodel3d_file);
                source = ld.sourcemodel3d;
                return
            end
        end

        %% See also
        %  https://www.fieldtriptoolbox.org/tutorial/visual_artifact_rejection/
        %  https://www.fieldtriptoolbox.org/tutorial/beamformer_lcmv/
        %  https://www.fieldtriptoolbox.org/tutorial/minimumnormestimate/#inverse-solution
        %  https://www.fieldtriptoolbox.org/workshop/paris2019/handson_sourceanalysis/#preprocessing
        %  https://www.fieldtriptoolbox.org/tutorial/networkanalysis/
        %  https://www.fieldtriptoolbox.org/assets/img/template/layout/4d248_helmet.mat.png

        function ft_databrowser(this)
            data = this.selected_rmegpreproc_data;
            cfg = [];
            ft_databrowser(cfg, data);
        end

        function tlck = ft_timelockanalysis(this, opts)
            %% build covariance matrix

            arguments
                this mlraut.FieldTrip
                opts.data struct {mustBeNonempty} = this.selected_rmegpreproc_data()
                opts.source_method {mustBeTextScalar} = "mne"
                opts.do_visualize logical = false
            end

            switch opts.source_method
                case "mne"
                    cfg = [];
                    cfg.covariance = 'yes';
                    tlck = ft_timelockanalysis(cfg, opts.data);
                case "lcmv"
                    cfg = [];
                    cfg.covariance = 'yes';
                    cfg.channel = 'A';
                    %cfg.vartrllength = 2;
                    cfg.covariancewindow = 'all';
                    cfg.lcmv.keepfilter = 'yes';
                    tlck = ft_timelockanalysis(cfg, opts.data);
                otherwise
                    error("mlraut:FieldTrip", stackstr())
            end


            if opts.do_visualize

                %% https://www.fieldtriptoolbox.org/tutorial/beamformer_lcmv/

                % MEG sensor covariance matrix
                figure;
                imagesc(tlck.cov);

                % simple timeseries plot
                figure;
                plot(tlck.time, tlck.avg(1, :));

                % view the results
                cfg        = [];
                cfg.layout = '4D248_helmet.mat';
                cfg.xlim   = [0.045 0.050];
                ft_topoplotER(cfg, tlck);

                % calculate planar gradients
                cfg                 = [];
                cfg.feedback        = 'yes';
                cfg.method          = 'template';
                cfg.template        = 'bti248_neighb.mat';
                cfg.neighbours      = ft_prepare_neighbours(cfg, tlck);

                cfg.planarmethod    = 'sincos';
                timelock_planar     = ft_megplanar(cfg, tlck);

                % combine planar gradients
                cfg                 = [];
                timelock_planarcomb = ft_combineplanar(cfg, timelock_planar);

                % view the results
                cfg        = [];
                cfg.layout = '4D248_helmet.mat';
                cfg.xlim   = [0.045 0.050];
                ft_topoplotER(cfg, timelock_planarcomb);
            end
        end

        function [leadfield,headmodel,tlck] = ft_prepare_leadfield(this, opts)
            %% Foward solution:  build leadfield using headmodel

            arguments
                this mlraut.FieldTrip
                opts.data struct {mustBeNonempty} = this.selected_rmegpreproc_data()
                opts.sourcemodel_kind {mustBeTextScalar} = "sourcemodel_3d4mm"  % "2d" or "3d" are minimal specifications
                opts.source_method {mustBeTextScalar} = "lcmv"
            end

            % build timelocked data with covariance matrix
            tlck = this.ft_timelockanalysis(data=opts.data, source_method=opts.source_method);

            cfg = [];
            cfg.grad = tlck.grad;  % gradiometer distances
            cfg.channel = tlck.label;  % selected channels
            cfg.grid = this.sourcemodel(kind=opts.sourcemodel_kind);  % source points
            ld = load(this.headmodel_file);
            headmodel = ld.headmodel;
            cfg.headmodel = headmodel;  % volume conduction headmodel
            cfg.singleshell.batchsize = 5000;  % speed-up
            leadfield = ft_prepare_leadfield(cfg);
        end

        function source = ft_sourceanalysis(this, opts)
            %% Inverse solution:  requires functional data, forward solution (leadfield), volume conduction model (headmodel),
            %  noise-covariance (from ft_timelockanalysis())

            arguments
                this mlraut.FieldTrip
                opts.leadfield struct = []
                opts.headmodel struct = []
                opts.timelock struct = [] 
                opts.method {mustBeTextScalar} = "lcmv"
            end
            if isempty(opts.leadfield) || isempty(opts.headmodel) || isempty(opts.timelock)
                 [opts.leadfield,opts.headmodel,opts.timelock] = this.ft_prepare_leadfield(source_method=opts.method);
            end

            %% create spatial filters for timeseries

            switch lower(opts.method)
                case "lcmv"

                    %% https://www.fieldtriptoolbox.org/tutorial/virtual_sensors/

                    cfg = [];
                    cfg.method = 'lcmv';
                    cfg.headmodel = opts.headmodel;
                    cfg.sourcemodel = opts.leadfield;
                    cfg.lcmv.keepfilter = 'yes';
                    cfg.lcmv.fixedori = 'yes';  % project on axis of most variance using SVD
                    source = ft_sourceanalysis(cfg, opts.timelock);

                    save(fullfile(this.megdir, stackstr()), "source");
                    
                case "mne"
                    cfg = [];
                    cfg.method = 'mne';
                    cfg.grid = opts.leadfield;
                    cfg.headmodel = opts.headmodel;
                    cfg.mne.prewhiten = 'yes';
                    cfg.mne.lambda = 3;
                    cfg.mne.scalesourcecov = 'yes';
                    source = ft_sourceanalysis(cfg, opts.timelock);

                    save(fullfile(this.megdir, stackstr()), "source");

                    % complex_signal = this.complex_signal(source);
                    % save(fullfile(this.megdir, "complex_signal.mat"), "complex_signal")

                    % m = source.avg.pow(:,450); % plotting the result at the 450th time-point that is
                    % ft_plot_mesh(source, 'vertexcolor', m);
                    % view([180 0]); h = light; set(h, 'position', [0 1 0.2]); lighting gouraud; material dull
                    
                case "dspm"
                case "sloretta"
                case "eloretta"
                otherwise
                    error("mlraut:ValueError", stackstr())
            end
        end

        function ft_sourceplot(this)
        end

        function spectral_analysis(this)
            %% https://www.fieldtriptoolbox.org/tutorial/networkanalysis/#spectral-analysis

            ld = load(this.rmegpreproc_files{1});
            dataica = ld.data;

            %% compute the power spectrum
            cfg              = [];
            cfg.output       = 'pow';
            cfg.method       = 'mtmfft';
            cfg.taper        = 'dpss';
            cfg.tapsmofrq    = 1;
            cfg.keeptrials   = 'no';
            datapow          = ft_freqanalysis(cfg, dataica);

            %% compute the planar transformation, this is not really necessary, but instructive anyhow
            ld = load(this.neighbours_file); % loads the neighbourhood structure for the channels

            dataicatmp      = dataica;
            dataicatmp.grad = dataica.grad;

            cfg               = [];
            cfg.neighbours    = ld.neighbours;
            cfg.planarmethod  = 'sincos';
            planar            = ft_megplanar(cfg, dataicatmp);
            clear dataicatmp;

            %% compute the power spectrum
            cfg              = [];
            cfg.output       = 'pow';
            cfg.method       = 'mtmfft';
            cfg.taper        = 'dpss';
            cfg.tapsmofrq    = 1;
            cfg.keeptrials   = 'no';
            datapow_planar   = ft_freqanalysis(cfg, planar);

            %% plot the topography and the spectrum
            figure;
            tiledlayout(2,2)

            cfg        = [];
            cfg.layout = '4D248_helmet.mat';
            % cfg.xlim   = [9 11];
            nexttile; ft_topoplotER(cfg, datapow);
            nexttile; ft_topoplotER(cfg, ft_combineplanar([], datapow_planar));

            cfg         = [];
            cfg.channel = {'A1', 'A12', 'A28'};
            nexttile; ft_singleplotER(cfg, datapow);
        end

        function neural_activity_index(this)
            %% https://www.fieldtriptoolbox.org/tutorial/networkanalysis/#source-reconstruction

            %% compute sensor level Fourier spectra, to be used for cross-spectral density computation.
            cfg            = [];
            cfg.method     = 'mtmfft';
            cfg.output     = 'fourier';
            cfg.keeptrials = 'yes';
            cfg.tapsmofrq  = 1;
            cfg.foi        = 10;
            freq           = ft_freqanalysis(cfg, dataica);

            %% do the source reconstruction
            cfg                   = [];
            cfg.frequency         = freq.freq;
            cfg.method            = 'pcc';
            cfg.grid              = lf;
            cfg.headmodel         = hdm;
            cfg.keeptrials        = 'yes';
            cfg.pcc.lambda        = '10%';
            cfg.pcc.projectnoise  = 'yes';
            cfg.pcc.fixedori      = 'yes';
            source = ft_sourceanalysis(cfg, freq);
            source = ft_sourcedescriptives([], source); % to get the neural-activity-index

            %% plot the neural activity index (power/noise)
            cfg               = [];
            cfg.method        = 'surface';
            cfg.funparameter  = 'nai';
            cfg.maskparameter = cfg.funparameter;
            cfg.funcolorlim   = [0.0 8];
            cfg.opacitylim    = [3 8];
            cfg.opacitymap    = 'rampup';
            cfg.funcolormap   = 'jet';
            cfg.colorbar      = 'no';
            ft_sourceplot(cfg, source);
            view([-90 30]);
            light;
        end
    end

    methods (Static)
        function source = trim_frames(source, N)
            arguments
                source struct {mustBeNonempty}
                N double {mustBeInteger}
            end
            assert(N < length(source.time))

            source.time = source.time(1:N);
            source.avg.pow = source.avg.pow(:, 1:N);
        end
    end

    %% PRIVATE

    methods (Access = private)
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
