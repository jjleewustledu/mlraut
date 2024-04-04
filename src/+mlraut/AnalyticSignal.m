classdef AnalyticSignal < handle & mlraut.HCP
    %% Implements Raut, et al.  Global waves synchronize the brain's functional systems with fluctuating arousal.
    %  https://www.science.org/doi/10.1126/sciadv.abf2709
    %  Extends physio_phase_mapping.m.
    %  Emphasizes consistency of analytic signal generated by Hilbert transform.
    %  
    %  Created 29-Nov-2022 13:47:07 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 9.13.0.2105380 (R2022b) Update 2 for MACI64.  Copyright 2022 John J. Lee.
    

    properties
        do_only_resting
        do_only_task
        do_plot_emd
        do_plot_global_physio
        do_plot_networks
        do_plot_radar
        do_save
        do_save_ciftis
        do_save_dynamic

        force_band  % force bandpass to [0.01 0.1] Hz
        final_normalization
        roi
        source_physio
    end

    properties (Dependent)
        global_signal
        global_signal_regression  % logical
        hp_thresh  % lower bound, Ryan ~ 0.01 Hz -> units of 1/frame_duration
        json
        lp_thresh  % higher bound, Ryan ~ 0.05 Hz -> units of 1/frame_duration
        no_physio
        num_nets
        num_sub
        num_tasks
        scale_to_hcp  % adjust norm by time of scanning
        tags  % for filenames

        analytic_signal
        bold_signal
        HCP_signals
        physio_signal
    end

    methods %% GET, SET
        function g = get.global_signal(this)
            g = this.global_signal_;
        end
        function g = get.global_signal_regression(this)
            g = this.global_signal_regression_;
        end
        function g = get.hp_thresh(this)
            N = this.num_frames - 2*this.num_frames_to_trim;
            Nyquist = 2/N; % Nyquist limited
            if this.force_band
                g = max(0.01*this.tr, Nyquist);
                return
            end
            if ~isempty(this.hp_thresh_)
                g = max(this.hp_thresh_, Nyquist);
                return
            end
            g = Nyquist; 
        end
        function g = get.json(this)
            g = this.cohort_data_.json;
        end
        function     set.json(this, s)
            this.cohort_data_.json = s;
        end
        function g = get.lp_thresh(this)
            Nyquist = 1/2; % Nyquist limited
            if this.force_band
                g = min(0.1*this.tr, Nyquist);
                return
            end
            if ~isempty(this.lp_thresh_)
                g = min(this.lp_thresh_, Nyquist);
                return
            end
            g = Nyquist;
        end        
        function g = get.no_physio(this)
            g = ~stricmp(this.source_physio, "none") && ~stricmp(this.source_physio, "nophy");
        end
        function g = get.num_nets(this)
            g = length(mlraut.NetworkData.NETWORKS_YEO_NAMES);
        end
        function g = get.num_sub(this)
            g = numel(this.subjects);
        end
        function g = get.num_tasks(this)
            g = numel(this.tasks);
        end        
        function g = get.scale_to_hcp(this)
            if ~isempty(this.scale_to_hcp_)
                g = this.scale_to_hcp_;
                return
            end

            this.scale_to_hcp_ = min(1192, this.max_frames)*0.72/(this.num_frames*this.tr);
            g = this.scale_to_hcp_;
        end
        function g = get.tags(this)
            if isempty(this.final_normalization) && strcmp(this.source_physio, "iFV") && ...
                    ~isempty(this.hp_thresh) && ~isempty(this.lp_thresh)
                % provide legacy compatibility
                g = "";
                return
            end

            g = "proc";
            if ~isempty(this.lp_thresh)
                g = g + "-lp" + strrep(num2str(this.lp_thresh), ".", "p");
            end
            if ~isempty(this.hp_thresh) 
                g = g + "-hp" + strrep(num2str(this.hp_thresh), ".", "p");
            end
            if ~isemptytext(this.final_normalization)
                g = g + "-" + this.final_normalization;
            end
            if ~isemptytext(this.source_physio)
                g = g + "-" + this.source_physio;
            end
            if isfinite(this.max_frames)
                g = g + "-maxframes" + num2str(this.max_frames);
            end
            if ~isempty(this.tags_) 
                g = g + "-" + this.tags_;
            end
        end

        function g = get.analytic_signal(this)
            g = this.analytic_signal_;
        end
        function g = get.bold_signal(this)
            g = this.bold_signal_;
        end
        function g = get.HCP_signals(this)
            g = this.HCP_signals_;
        end
        function g = get.physio_signal(this)
            g = this.physio_signal_;
        end
    end
    
    methods
        function a = angle(~, as)
            a = unwrap(angle(as)); % [-pi pi] -> [-inf inf]
            a = mod(a, 2*pi); % [-inf inf] -> [0 2*pi]
        end
        function psis = average_network_signals(this, psi)
            cbm = mlraut.CerebellarData(this, psi);
            this.HCP_signals_.cbm = cbm.build_Yeo_signals();
            ctx = mlraut.CorticalData(this, psi);
            this.HCP_signals_.ctx = ctx.build_Yeo_signals();
            str = mlraut.StriatalData(this, psi);
            this.HCP_signals_.str = str.build_Yeo_signals();
            thal = mlraut.ThalamicData(this, psi);
            this.HCP_signals_.thal = thal.build_Yeo_signals();

            psis = this.HCP_signals;
        end
        function this = call(this, opts)
            %% CALL all subjects

            arguments
                this mlraut.AnalyticSignal
                opts.do_qc logical = false
            end

            % exclude subjects
            this.subjects = this.subjects(~contains(this.subjects, '_7T'));
            %this.subjects = this.subjects(~contains(this.subjects, 'sub-'));

            out_dir_ = this.out_dir;
            for s = 1:this.num_sub
                try
                    this.current_subject = this.subjects{s};
                    if ~contains(out_dir_, this.current_subject)
                        proposed_dir = fullfile(out_dir_, this.current_subject);
                        this.out_dir = proposed_dir;
                        ensuredir(proposed_dir);
                    end
                    if opts.do_qc
                        this.call_subject_qc(s);
                    else
                        this.call_subject(s);
                    end
                catch ME
                    handexcept(ME)
                end
            end
        end
        function this = call_subject(this, s)
            arguments
                this mlraut.AnalyticSignal
                s double
            end

            this.malloc();
            for t = 1:this.num_tasks     
                try
                    this.current_task = this.tasks{t};

                    % BOLD
                    try
                        bold_ = this.task_dtseries();
                        bold_ = this.build_global_signal_regressed(bold_);
                        bold_ = hilbert(bold_);
                    catch ME
                        disp([this.current_subject ' ' this.current_task ' BOLD missing or defective:']);
                        handwarning(ME)
                        continue
                    end

                    % Physio
                    try
                        physio_ = this.task_physio();
                        physio_ = this.build_global_signal_regressed(physio_);
                        physio_ = hilbert(physio_);
                    catch ME
                        disp([this.current_subject ' ' this.current_task ' physio missing or defective:']);
                        handwarning(ME)
                        continue
                    end

                    % Store BOLD, Physio, and Analytic signals
                    this.bold_signal_ = this.build_band_passed(this.build_centered_and_rescaled(bold_));
                    this.bold_signal_ = this.build_final_normalization(this.bold_signal_);
                    if ~all(physio_ == 0)
                        this.physio_signal_ = this.build_band_passed(this.build_centered_and_rescaled(physio_));
                        this.physio_signal_ = this.build_final_normalization(this.physio_signal_);
                        % <psi_p|BOLD_operator|psi_p> ~ <psi_p|psi_b>, not unitary
                        this.analytic_signal_ = this.build_band_passed( ...
                            this.build_centered_and_rescaled(conj(physio_)) .* ...
                            this.build_centered_and_rescaled(bold_));
                        this.analytic_signal_ = this.build_final_normalization(this.analytic_signal_);
                    else
                        this.physio_signal_ = physio_;
                        this.analytic_signal_ = this.bold_signal_;
                    end

                    % Averages for networks
                    this.average_network_signals(this.analytic_signal_);

                    % Store reduced analytic signal, real(), imag(), abs(), angle()
                    if this.do_save
                        % Store reduced analytic signals for all s, t
                        % grid of data from s, t may be assessed with stats
                        save(this, s, t);
                    end
                    if this.do_save_ciftis
                        this.write_ciftis( ...
                            abs(this.analytic_signal_), ...
                            sprintf('abs_as_sub-%s_ses-%s_%s', this.subjects{s}, this.tasks{t}, this.tags), ...
                            do_save_dynamic=this.do_save_dynamic);
                        this.write_ciftis( ...
                            angle(this.analytic_signal_), ...
                            sprintf('angle_as_sub-%s_ses-%s_%s', this.subjects{s}, this.tasks{t}, this.tags), ...
                            do_save_dynamic=this.do_save_dynamic);

                        % analytic_signal_ - bold_signal_
                        diff_ = this.analytic_signal_ - this.bold_signal_;
                        this.write_ciftis( ...
                            abs(diff_), ...
                            sprintf('abs_diff_sub-%s_ses-%s_%s', this.subjects{s}, this.tasks{t}, this.tags), ...
                            do_save_dynamic=this.do_save_dynamic);
                        this.write_ciftis( ...
                            angle(diff_), ...
                            sprintf('angle_diff_sub-%s_ses-%s_%s', this.subjects{s}, this.tasks{t}, this.tags), ...
                            do_save_dynamic=this.do_save_dynamic);
                    end

                    % do plot
                    if this.do_plot_global_physio
                        this.plot_global_physio(measure=@this.unwrap);
                        this.plot_global_physio(measure=@angle);
                        this.plot_global_physio(measure=@abs);
                        this.plot_global_physio(measure=@real);
                    end
                    if this.do_plot_networks
                        this.plot_regions(@this.plot_networks, measure=@this.unwrap);
                        this.plot_regions(@this.plot_networks, measure=@angle);
                        this.plot_regions(@this.plot_networks, measure=@abs);
                        this.plot_regions(@this.plot_networks, measure=@real);
                    end
                    if this.do_plot_radar
                        this.plot_regions(@this.plot_radar, measure=@this.identity);
                    end
                    if this.do_plot_emd
                        this.plot_regions(@this.plot_emd, measure=@this.unwrap);
                        this.plot_regions(@this.plot_emd, measure=@angle);
                        this.plot_regions(@this.plot_emd, measure=@abs);
                    end
                catch ME
                    handwarning(ME)
                end
            end
        end
        function this = call_subject_qc(this, s, t)
            arguments
                this mlraut.AnalyticSignal
                s double = 1
                t double = 1
            end

            this.malloc();            
            this.current_task = this.tasks{t};
            tic

            % BOLD
            bold_ = this.task_dtseries(); 
            bold_ = this.build_global_signal_regressed(bold_);
            this.write_ciftis(bold_, "bold-subgs-"+this.global_signal_regression_);
             
            % Physio
            physio_ = this.task_physio();
            physio_ = this.build_global_signal_regressed(physio_);
            this.write_nii(physio_, this.source_physio+"-subgs-"+this.global_signal_regression_);

            % BOLD for Analytic signal
            psi = this.build_centered(bold_);
            this.write_ciftis(psi, "centered-bold-subgs-"+this.global_signal_regression_);
            psi = this.build_rescaled(psi);
            this.write_ciftis(psi, "rescaled-centered-bold-subgs-"+this.global_signal_regression_);
            psi = this.build_band_passed(psi);
            this.write_ciftis(psi, "lp"+this.lp_thresh+"-hp"+this.hp_thresh+"-rescaled-centered-bold-subgs-"+this.global_signal_regression_);
            psi = hilbert(psi);
            this.write_ciftis(abs(psi), "abs-psi");
            this.write_ciftis(angle(psi), "angle-psi");

            % Physio for Analytic signal
            phi = this.build_centered(physio_);
            this.write_nii(phi, "centered-"+this.source_physio+"-subgs-"+this.global_signal_regression_);
            phi = this.build_rescaled(phi);
            this.write_nii(phi, "rescaled-centered-"+this.source_physio+"-subgs-"+this.global_signal_regression_);
            phi = this.build_band_passed(phi);
            this.write_nii(phi, "lp"+this.lp_thresh+"-hp"+this.hp_thresh+"-rescaled-centered-"+this.source_physio+"-subgs-"+this.global_signal_regression_);
            phi = hilbert(phi);
            this.write_nii(abs(phi), "abs-phi");
            this.write_nii(angle(phi), "angle-phi");

            % Analytic signal aufbau
            as = conj(phi) .* psi;
            this.write_ciftis(abs(as), "abs-as");
            this.write_ciftis(angle(as), "angle-as");
            as = this.build_final_normalization(as); 
            if ~isempyttext(this.final_normalization)
                this.write_ciftis(abs(as), "abs-"+this.final_normalization+"-as");
                this.write_ciftis(angle(as), "angle-"+this.final_normalization+"-as");
            end

            this.save(s, t);

            toc
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
            [z,p,k] = butter(2, [this.hp_thresh/2, 2*this.lp_thresh - eps('single')]); % digital Wn in [0, 1]
            [sos,g] = zp2sos(z, p, k);
            dat1 = filtfilt(sos, g, double(dat));
            if isa(dat, 'single')
                dat1 = single(dat1);
            end
            if isa(dat, 'double')
                dat1 = double(dat1);
            end
        end
        function psi = build_centered(this, psi)
            assert(~isempty(psi))
            if all(psi == 0)
                return
            end
            psi = psi - median(psi, 'all', 'omitnan');
        end
        function psi = build_centered_and_rescaled(this, psi)
            %% Mimics z-score of |psi(t,x)> using median and mad.

            psi = this.build_centered(psi);
            psi = this.build_rescaled(psi);
        end
        function as = build_final_normalization(this, as)
            %% provides final normalization by max(abs()) for more interpretable visualization

            switch convertStringsToChars(this.final_normalization)
                case 'normt'
                    % allowing fluctuations in xyz
                    as = as ./ max(abs(as), 1);
                case 'normxyz'
                    % allowing fluctuations in t
                    as = as ./ max(abs(as), 2);
                case 'normxyzt'
                    % allowing fluctuations in xyz & t
                    as = as / max(abs(as), [], "all");
                otherwise
                    return
            end
        end
        function [gs,beta] = build_global_signal_for(this, sig)
            %% global_signal := median(sig, "space"), then formatted for greyordinates or 4D voxels

            % task niigz in 4D, reshaped to 2D 
            niigz = this.task_niigz();
            sz = size(niigz);
            assert(sz(4) == this.num_frames)
            niigz = reshape(niigz, [sz(1)*sz(2)*sz(3), sz(4)]);

            % mask ~ 3D task signal reference; reshaped to 2D
            msk = this.task_signal_mask();
            msk = reshape(msk, [sz(1)*sz(2)*sz(3), 1]);

            % img ~ task niigz masked; then median
            img_g = double(niigz);
            img_g = img_g(logical(msk), :);
            img_g = median(img_g, 1);

            % format for greyordinates
            this.global_signal_ = ascol(img_g);
            this.global_signal_beta_ = 1;

            if isnumeric(sig)
                assert(size(sig, 1) == this.num_frames)
                gs = this.global_signal_;  % col
                beta = this.global_signal_beta_;  % row
                return
            end
            if isa(sig, "mlfourd.ImagingContext2")
                sz = size(sig);
                assert(sz(4) == this.num_frames)
                sig = reshape(sig, [sz(1)*sz(2)*sz(3), sz(4)]);

                msk = this.task_signal_mask();
                msk = reshape(msk, [sz(1)*sz(2)*sz(3), 1]);

                img_g = double(sig);
                img_g(logical(msk), :) = repmat(asrow(this.global_signal_), [sum(logical(msk)), 1]);
                img_g = reshape(img_g, [sz(1), sz(2), sz(3), sz(4)]);
  
                img_b = double(sig);
                img_b(logical(msk), :) = this.global_signal_beta_ * ones(sum(logical(msk), sz(4)));
                img_b = reshape(img_b, [sz(1), sz(2), sz(3), sz(4)]);

                gs = sig.selectImagingTool(img=img_g);
                gs.fileprefix = "global-signal";
                beta = sig.selectImagingTool(img=img_b);
                beta.fileprefix = "global-signal-beta";
                return
            end
            error("mlraut:TypeError", stackstr())
        end
        function psi = build_global_signal_regressed(this, psi)
            if ~this.global_signal_regression
                return
            end
            if all(psi == 0)
                return
            end

            psi = psi - this.build_global_signal_for(psi);
        end
        function psi = build_rescaled(this, psi)
            assert(~isempty(psi))
            if all(psi == 0)
                return
            end

            d = mad(abs(psi), 1, 'all');  % median abs. dev.
            psi = psi./d;
        end
        function obj = identity(~, obj)
        end
        function this = malloc(this)

            % accumulate for statistics on serialized AnalyticSignal
            this.bold_signal_ = complex(nan(this.num_frames, this.num_nodes));  % largest
            this.physio_signal_ = complex(nan(this.num_frames, 1));   
            this.analytic_signal_ = complex(nan(this.num_frames, this.num_nodes));  % largest

            this.HCP_signals_.cbm = complex(nan(this.num_frames,this.num_nets));
            this.HCP_signals_.ctx = complex(nan(this.num_frames,this.num_nets));
            this.HCP_signals_.str = complex(nan(this.num_frames,this.num_nets));
            this.HCP_signals_.thal = complex(nan(this.num_frames,this.num_nets));
        end
        function b = omit_late_frames(this, b)
            %% Keep frames 1:this.max_frames, following use of trim_frames() to remove this.num_frames_to_trim
            %  from start and end of frames, for purposes of omitting brain/cognitive responses to start and conclusion 
            %  of the scanning session.            

            if isnumeric(b)
                bound = min(this.max_frames, size(b, 1));
                b = b(1:bound, :);
                return
            end
            if isa(b, "mlfourd.ImagingContext2")
                bound = min(this.max_frames, size(b, 4));
                img = double(b);
                img = img(:,:,:,1:bound);
                b.selectImagingTool(img=img);
                j = b.json_metadata;
                j.timesMid = j.timesMid(1:bound);
                b.addJsonMetadata(j);
                return
            end
            error("mlraut:TypeError", stackstr())
        end

        %% PLOTTING

        function plot_emd(this, varargin)
            this.plotting_.plot_emd(varargin{:});
        end
        function h1 = plot_global_physio(this, varargin)
            this.plotting_.plot_global_physio(varargin{:});
        end
        function plot_regions(this, varargin)
            this.plotting_.plot_regions(varargin{:});
        end
        function [h1,h3] = plot_networks(this, varargin)
            [h1,h3] = this.plotting_.plot_networks(varargin{:});
        end
        function [h,h1,h2] = plot_radar(this, varargin)
            [h,h1,h2] = this.plotting_.plot_radar(varargin{:});
        end
        function [h,h1] = plot_timeseries_qc(this, varargin)
            [h,h1] = this.plotting_.plot_timeseries_qc(varargin{:});
        end
        
        %%

        function save(this, s, t)

            % reduce size of saved            
            this_subset.do_only_resting = this.do_only_resting;
            this_subset.do_only_task = this.do_only_task;
            this_subset.do_save = this.do_save;
            this_subset.do_save_ciftis = this.do_save_ciftis;
            this_subset.do_save_dynamic = this.do_save_dynamic;
            this_subset.force_band = this.force_band;
            this_subset.final_normalization = this.final_normalization;
            this_subset.roi = this.roi;
            this_subset.source_physio = this.source_physio;
            this_subset.global_signal = this.global_signal;
            this_subset.global_signal_regression = this.global_signal_regression;
            this_subset.hp_thresh = this.hp_thresh;
            this_subset.lp_thresh = this.lp_thresh;
            this_subset.num_nets = this.num_nets;
            this_subset.num_sub = this.num_sub;
            this_subset.num_tasks = this.num_tasks;
            this_subset.scale_to_hcp = this.scale_to_hcp;
            this_subset.tags = this.tags;
            this_subset.analytic_signal = this.analytic_signal;
            % this_subset.bold_signal = this.bold_signal;
            this_subset.HCP_signals = this.HCP_signals;
            this_subset.physio_signal = this.physio_signal;
            this_subset.max_frames = this.max_frames;
            this_subset.current_subject = this.current_subject;
            this_subset.current_task = this.current_task;
            this_subset.subjects = this.subjects;
            this_subset.tasks = this.tasks;
            % this_subset.bold_data = this.bold_data;
            % this_subset.cohort_data = this.cohort_data;
            % this_subset.cifti_last = this.cifti_last;
            this_subset.Fs = this.Fs;
            this_subset.num_frames = this.num_frames;
            this_subset.num_frames_ori = this.num_frames_ori;
            this_subset.num_frames_to_trim = this.num_frames_to_trim;
            this_subset.num_nodes = this.num_nodes;
            this_subset.out_dir = this.out_dir;
            this_subset.root_dir = this.root_dir;
            this_subset.task_dir = this.task_dir;
            this_subset.task_dtseries_fqfn = this.task_dtseries_fqfn;
            this_subset.task_niigz_fqfn = this.task_niigz_fqfn;
            this_subset.task_signal_reference_fqfn = this.task_signal_reference_fqfn;
            this_subset.t1w_fqfn = this.t1w_fqfn;
            this_subset.tr = this.tr;
            this_subset.waves_dir = this.waves_dir;
            this_subset.wmparc_fqfn = this.wmparc_fqfn;
            this_subset.workbench_dir = this.workbench_dir;

            the_tags_ = this.tags;
            the_out_dir_ = this.out_dir;
            
            try
                save(fullfile(the_out_dir_, ...
                    sprintf("sub-%s_ses-%s_%s.mat", this.subjects{s}, strrep(this.tasks{t}, "_", "-"), the_tags_)), ...
                    'this_subset');
            catch ME
                handwarning(ME)
            end
        end
        function mat = task_dtseries(this, varargin)
            %  Args:
            %      this mlraut.AnalyticSignal
            %      sub {mustBeTextScalar} = this.current_subject
            %      task {mustBeTextScalar} = this.current_task
            %  Returns:
            %      mat (numeric):  time x grayordinate from BOLDData

            mat = task_dtseries@mlraut.HCP(this, varargin{:});
            mat = this.trim_frames(mat);
            mat = this.omit_late_frames(mat);
            mat = single(mat);
        end
        function ic = task_niigz(this)
            ic = task_niigz@mlraut.HCP(this);
            ic = this.trim_frames(ic);
            ic = this.omit_late_frames(ic);
            ic.ensureSingle();
            ic.fileprefix = stackstr(use_dashes=true);
        end
        function physio = task_physio(this, opts)
            %  Returns:
            %      physio numeric Nt x Nx
            %  Throws:
            %      mlraut:ValueError if this.source_physio not supported

            arguments
                this mlraut.AnalyticSignal
                opts.roi = this.roi
                opts.flipLR logical = false
                opts.source_physio {mustBeText} = this.source_physio
            end

            bold = this.task_niigz();
            switch convertStringsToChars(this.source_physio)
                case 'RV'
                    RV = mlraut.PhysioRV(this, bold);
                    physio = RV.call();
                case 'HRV'
                    HRV = mlraut.PhysioHRV(this, bold);
                    physio = HRV.call();
                case 'iFV'
                    iFV = mlraut.IFourthVentricle(this, bold);
                    physio = iFV.call();
                case 'ROI'
                    pROI = mlraut.PhysioRoi(this, bold, ...
                        from_imaging_context=opts.roi, flipLR=opts.flipLR);
                    physio = pROI.call();
                case {'no-physio', 'nophys', 'none'}
                    physio = zeros(size(bold, ndims(bold)), 1);
                otherwise
                    error("mlraut:ValueError", stackstr())
            end
            physio = single(physio);
            assert(~isempty(physio))
        end
        function ic = task_signal_mask(this)
            % if ~isempty(this.task_signal_mask_)
            %     ic = this.task_signal_mask_;
            %     return
            % end

            ic = this.task_signal_reference();
            ic = ic.binarized();
            
            ic1 = mlfourd.ImagingContext2(this.wmparc_fqfn);
            ic1 = ic1.binarized();
            ic1 = ic1.blurred(6);
            ic1 = ic1.thresh(0.1);
            ic1 = ic1.binarized();
            if ~isempty(getenv("DEBUG"))
                ic1.view_qc(ic);
            end

            % mask should not have greater coverage than blurred binarized wmparc
            if dipsum(ic) > dipsum(ic1)
                ic = ic1;
            end
            ic.ensureSingle();
            this.task_signal_mask_ = ic;
        end
        function ic = task_signal_reference(this)
            ic = task_signal_reference@mlraut.HCP(this);
            ic.ensureSingle();
        end
        function tseries = trim_frames(this, tseries)
            nt = this.num_frames_to_trim + 1;
            if isnumeric(tseries)
                tseries = tseries(nt:end-nt+1,:);
                return
            end
            if isa(tseries, "mlfourd.ImagingContext2")
                img = double(tseries);
                img = img(:,:,:,nt:end-nt+1);
                tseries.selectImagingTool(img=img);
                j = tseries.json_metadata;
                j.timesMid = j.timesMid(nt:end-nt+1);
                tseries.addJsonMetadata(j);
                return
            end
            error("mlraut:TypeError", stackstr())
        end
        function u = unwrap(~, psi)
            u = unwrap(angle(psi));
        end
        function write_ciftis(this, varargin)
            this.cifti_.write_ciftis(varargin{:});
        end
        function write_nii(this, varargin)
            this.cifti_.write_nii(varargin{:});
        end

        function this = AnalyticSignal(opts)
            %% ANALYTICSIGNAL 
            %  Args:
            %      opts.do_only_resting logical = true
            %      opts.do_only_tast logical = false
            %      opts.do_plot_global_physio logical = true
            %      opts.do_plot_networks logical = true
            %      opts.do_plot_radar logical = true
            %      opts.do_save logical = true: save fully populated this to mlraut_AnalyticSignal.mat
            %      opts.do_save_ciftis logical = true: save ciftis of {abs,angle} of analytic_signal.
            %      opts.do_save_dynamic logical = false; save large dynamic dtseries
            %      opts.final_normalization {mustBeTextScalar} = 'normxyzt': also: 'normt' | 'normxyz' | ''
            %      opts.force_band logical = tru e: force bandpass to [0.01 0.1] Hz
            %      opts.hp_thresh {mustBeScalarOrEmpty} : default := 0.009*0.72, Dworetsky; support ~ 2/this.num_frames ~ 0.0019, compared to Ryan's 0.01.
            %                                             nan =: 2/(this.num_frames - this.num_frames_to_trim).
            %      opts.lp_thresh {mustBeScalarOrEmpty} : default := 0.08*0.72, Dworetsky; support ~ 1/(2*this.tr), compared to Ryan's 0.05.
            %                                             nan =: 1/2
            %      opts.max_frames {mustBeScalarOrEmpty} = nan: try 158 for assessing GBM rsfMRI
            %      opts.out_dir {mustBeFolder} = pwd
            %      opts.roi = []:  e.g. fqfn;
            %                      ImagingContext2 for "WT_on_T1w", "CE_on_T1w", for files found in sub-*/MNINonLinear; 
            %                      double row_vec for mlraut.PhysioRoi(from_wmparc_indices=row_vec)
            %                      used with this.task_physio()
            %      opts.scale_to_hcp {mustBeScalar,mustBePositive} = 1: scaling factor
            %      opts.source_physio {mustBeTextScalar} = 'iFV'
            %      opts.global_signal_regression logical = true
            %      opts.subjects cell {mustBeText} = {}
            %      opts.tags {mustBeTextScalar} = ""
            %      opts.tasks cell {mustBeText} = {}
            
            arguments
                opts.do_only_resting logical = true
                opts.do_only_task logical = false
                opts.do_plot_emd logical = false
                opts.do_plot_global_physio logical = false
                opts.do_plot_networks logical = false
                opts.do_plot_radar logical = false
                opts.do_save logical = false
                opts.do_save_ciftis logical = false
                opts.do_save_dynamic logical = false
                opts.final_normalization {mustBeTextScalar} = "normxyzt"
                opts.force_band logical = true
                opts.global_signal_regression logical = true
                opts.hp_thresh {mustBeScalarOrEmpty} = []
                opts.lp_thresh {mustBeScalarOrEmpty} = []
                opts.max_frames double = Inf
                opts.out_dir {mustBeTextScalar} = ""
                opts.plot_range double = 1:572
                opts.roi = []
                opts.scale_to_hcp double {mustBePositive} = 1
                opts.source_physio {mustBeTextScalar} = "iFV"
                opts.subjects = {}
                opts.tags {mustBeTextScalar} = ""
                opts.tasks = {}
            end

            this = this@mlraut.HCP(max_frames=opts.max_frames, subjects=opts.subjects, tasks=opts.tasks)

            addpath(genpath(fullfile(this.waves_dir, 'Dependencies', '-end')));
            addpath(genpath(fullfile(this.waves_dir, 'supporting_files', '')));
            this.cifti_ = mlraut.Cifti(this);
            this.plotting_ = mlraut.Plotting(this, plot_range=opts.plot_range);

            this.do_only_resting = opts.do_only_resting;
            this.do_only_task = opts.do_only_task;
            this.do_plot_emd = opts.do_plot_emd;
            this.do_plot_global_physio = opts.do_plot_global_physio;
            this.do_plot_networks = opts.do_plot_networks;
            this.do_plot_radar = opts.do_plot_radar;
            this.do_save = opts.do_save;
            this.do_save_ciftis = opts.do_save_ciftis;
            this.do_save_dynamic = opts.do_save_dynamic;

            this.force_band = opts.force_band;
            this.global_signal_regression_ = opts.global_signal_regression;
            this.hp_thresh_ = opts.hp_thresh;
            this.lp_thresh_ = opts.lp_thresh;
            this.max_frames = opts.max_frames;
            this.final_normalization = opts.final_normalization;
            this.cohort_data_.out_dir = opts.out_dir;
            this.build_roi(opts.roi);
            this.scale_to_hcp_ = opts.scale_to_hcp;
            this.source_physio = opts.source_physio;
            this.tags_ = opts.tags;

            this.debugging_ = struct();
        end
    end

    %% PROTECTED

    properties (Access = protected)
        cifti_
        debugging_
        hp_thresh_
        lp_thresh_
        plotting_
        scale_to_hcp_
        tags_
        task_signal_mask_

        analytic_signal_
        bold_signal_
        global_signal_
        global_signal_beta_
        HCP_signals_
        physio_signal_
        global_signal_regression_
    end

    methods (Access = protected)
        function that = copyElement(this)
            that = copyElement@matlab.mixin.Copyable(this);
            if ~isempty(this.cifti_)
                that.cifti_ = copy(this.cifti_); end
            if ~isempty(this.plotting_)
                that.plotting_ = copy(this.plotting_); end
            if ~isempty(this.task_signal_mask_)
                that.task_signal_mask_ = copy(this.task_signal_mask_); end
        end
    end

    %% PRIVATE

    methods (Access = private)
        function this = build_roi(this, roi)
            arguments
                this mlraut.AnalyticSignal
                roi = []
            end
            if isempty(roi)
                this.roi = [];
                return
            end
            
            this.source_physio = "ROI";
            if istext(roi) && isfile(roi)
                this.roi = mlfourd.ImagingContext2(roi);
                return
            end
            if isnumeric(roi) && ismatrix(roi)
                pr = mlraut.PhysioRoi(this, this.task_niigz, from_wmparc_indices=roi);
                this.roi = pr.roi_mask;
                return
            end
            if isnumeric(roi) && ~ismatrix(roi)
                this.roi = mlfourd.ImagingContext2(roi);
                return
            end
            error("mlraut:ValueError", stackstr())
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
