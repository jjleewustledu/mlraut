classdef Lee2024 < handle
    %% line1
    %
    % [~,min_idx] = min(t.OS_Diagnosis)
    % min_idx =
    %    115
    %
    % [~,max_idx] = max(t.OS_Diagnosis)
    % max_idx =
    %     41
    %
    % [t.I3CRID(115), t.I3CRID(41)]
    % ans =
    %   1×2 string array
    %     "I3CR1156"    "I3CR0433"
    %  
    % [t.OS_Diagnosis(115), t.OS_Diagnosis(41)]
    % ans =
    %           18        2696
    %
    %  Created 10-Feb-2024 13:43:26 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    
    properties
        analytic_signal_hcpaging_dir
        analytic_signal_gbm_dir
        gbm_datashare_dir
        matlabout_dir
        hcpaging_connex
        hcpaging_angle
        hcpaging_T
        hcpaging_X
        hcpaging_Y
        hcpaging_Z
    end

    methods
        function this = Lee2024(varargin)
            this.analytic_signal_hcpaging_dir = ...
                fullfile(getenv("SINGULARITY_HOME"), "AnalyticSignalHCPAging");
            this.analytic_signal_gbm_dir = ...
                fullfile(getenv("SINGULARITY_HOME"), "AnalyticSignalGBM");
            this.gbm_datashare_dir = ...
                fullfile(this.analytic_signal_gbm_dir, "GBM_datashare");
            this.matlabout_dir = ...
                fullfile(this.analytic_signal_gbm_dir, "analytic_signal", "matlabout");

            this.hcpaging_connex = fullfile( ...
                this.analytic_signal_hcpaging_dir, ...
                "mean_comparator_as_sub-all_ses-all_proc-v50-iFV-brightest-mean-twistor-nlim725_avgt.dscalar.nii");
            this.hcpaging_angle = fullfile(  ...
                this.analytic_signal_hcpaging_dir, ...
                "mean_angle_as_sub-all_ses-all_proc-v50-iFV-brightest-mean-element-nlim725_avgt.dscalar.nii");
            this.hcpaging_T = fullfile( ...
                this.analytic_signal_hcpaging_dir, ...
                "mean_T_as_sub-all_ses-all_proc-v50-iFV-brightest-mean-element-nlim725_avgt.dscalar.nii");
            this.hcpaging_X = fullfile( ...
                this.analytic_signal_hcpaging_dir, ...
                "mean_X_as_sub-all_ses-all_proc-v50-iFV-brightest-mean-twistor-rsn7-nlim725_avgt.dscalar.nii");
            this.hcpaging_Y = fullfile( ...
                this.analytic_signal_hcpaging_dir, ...
                "mean_Y_as_sub-all_ses-all_proc-v50-iFV-brightest-mean-twistor-rsn7-nlim725_avgt.dscalar.nii");
            this.hcpaging_Z = fullfile( ...
                this.analytic_signal_hcpaging_dir, ...
                "mean_Z_as_sub-all_ses-all_proc-v50-iFV-brightest-mean-element-nlim725_avgt.dscalar.nii");
        end

        function build_mean_for_gbm(this, opts)
            %  Lee2024_build_mean_for_gbm: Nerr->0, Nsubs->203, for iFV-brightest
            %  Lee2024_build_mean_for_gbm: Nerr->46, Nsubs->203, for CE
            %  Lee2024_build_mean_for_gbm: Nerr->47, Nsubs->203, for edema

            arguments
                this mlraut.Lee2024
                opts.physio {mustBeTextScalar} = "iFV-brightest"
            end

            data = {};
            Nsubs = numel(this.unique_subs);
            Nerr = 0;
            for s = this.unique_subs
                try
                    [datum,ld] = this.gather_signals(s, weight=1, build_residuals=false, physio=opts.physio);
                    datum.sub = s;
                    data = [data, {datum}]; %#ok<AGROW>
                catch ME
                    Nerr = Nerr + 1;
                    fprintf("%s:%s for %s\n", stackstr(), ME.identifier, s)
                    handwarning(ME)
                end
            end

            % re-weight using Nerr
            if Nerr > 0
                fprintf("%s: Nerr->%g, Nsubs->%g\n", stackstr(), Nerr, Nsubs)
            end

            save(fullfile(this.matlabout_dir, stackstr()+".mat"), "data");
            this.write_mean_ciftis(ld, data, "connex", physio=opts.physio);
            this.write_mean_ciftis(ld, data, "angle", physio=opts.physio);
            this.write_mean_ciftis(ld, data, "T", physio=opts.physio);
            this.write_mean_ciftis(ld, data, "X", physio=opts.physio);
            this.write_mean_ciftis(ld, data, "Y", physio=opts.physio);
            this.write_mean_ciftis(ld, data, "Z", physio=opts.physio);
        end

        function build_var_for_gbm(this, opts)
            arguments
                this mlraut.Lee2024
                opts.physio {mustBeTextScalar} = "iFV-brightest"
            end

            build_mean_for_gbm_mat = fullfile(this.matlabout_dir, "Lee2024_build_mean_for_gbm.mat");
            assert(isfile(build_mean_for_gbm_mat))
            ld_ = load(build_mean_for_gbm_mat, "data");
            data = ld_.data;
            ld = load( ...
                fullfile(this.matlabout_dir, ...
                "sub-I3CR0000", ...
                "sub-RT089_ses-1-task-rest-run-01-desc-preproc_proc-v50-"+opts.physio+"-AnalyticSignalGBMPar.mat"));

            this.write_var_ciftis(ld, data, "connex", physio=opts.physio);
            this.write_var_ciftis(ld, data, "angle", physio=opts.physio);
            this.write_var_ciftis(ld, data, "T", physio=opts.physio);
            this.write_var_ciftis(ld, data, "X", physio=opts.physio);
            this.write_var_ciftis(ld, data, "Y", physio=opts.physio);
            this.write_var_ciftis(ld, data, "Z", physio=opts.physio);
        end

        function build_uthr_dscalar(this, fn, alpha)
            arguments
                this mlraut.Lee2024 %#ok<INUSA>
                fn {mustBeFile}
                alpha {mustBeScalarOrEmpty}
            end
            fn1 = extractBefore(fn, "Lee2024");
            alpha_str = strrep(num2str(alpha), ".", "p");
            fn1 = sprintf("%s%s-alpha%s.dscalar.nii", fn1, stackstr(use_dashes=true), alpha_str);
            if alpha > 1
                alpha = 0.01*alpha;
            end

            cifti = cifti_read(convertStringsToChars(fn));
            c = cifti.cdata';
            uthr = prctile(c, 100*(1 - alpha), method="exact");
            c1 = c > uthr;
            cifti.cdata = c1';
            cifti_write(cifti, convertStringsToChars(fn1))
        end

        function concat_ciftify_runs(this, results_path)
        end

        function ifc = concat_nii(~, varargin)
            ifc = mlraut.BOLDData.concat_nii(varargin{:});
        end

        function c = concat_dtseries(~, varargin)
            c = mlraut.BOLDData.concat_dtseries(varargin{:});
        end

        function that = concat_frames_and_save(this, srcdir, opts)
            %% assumes run-01, run-02, run-03, run-all may exist

            arguments
                this mlraut.Lee2024  %#ok<INUSA>
                srcdir {mustBeFolder} = pwd
                opts.physios {mustBeText} = ["CE_on_T1w", "iFV-brightest"]
                opts.do_save logical = true
                opts.do_save_ciftis logical = false
                opts.do_save_dynamic logical = false
            end

            pwd0 = pushd(srcdir);

            for phys = opts.physios
                g = mglob(sprintf("sub-*-run-0*%s*.mat", phys));
                g = g(~contains(g, "-concat"));
                if isempty(g)
                    continue
                end

                % single file -> rename to run-all
                if isscalar(g)
                    copyfile(g, regexprep(g, 'run-\d+', 'run-all'));
                    continue
                end

                % multiple files -> invoke AnalyticSignalHCP.concat_frames()
                ld1 = load(g(1));  % adjust first run
                that = ld1.this;
                template_cifti_ = that.template_cifti;
                if isscalar(that.tasks)
                    that.tasks = ensureCell(regexprep(that.tasks, 'run-\d+', 'run-all'));
                    that.current_task = that.tasks{1};
                else
                    that.tasks = ensureCell(regexprep(that.tasks{1}, 'run-\d+', 'run-all'));
                    that.current_task = that.tasks{1};
                end
                for gidx = 2:length(g)  % concat subsequent runs into that
                    ld_ = load(g(gidx));
                    that.concat_frames(ld_.this);
                end
                
                that.template_cifti = template_cifti_;
                that.out_dir = srcdir;
                that.do_save = opts.do_save;
                that.do_save_ciftis = opts.do_save_ciftis;
                that.do_save_dynamic = opts.do_save_dynamic;
                if any([opts.do_save, opts.do_save_ciftis, opts.do_save_dynamic])
                    that.meta_save();
                end
            end

            popd(pwd0)
        end

        function tag = i3cr(this)
            subs_tag = ascol(extractAfter(this.subs, "-"));
            tag = subs_tag;
            t = this.table_found_RT_I3CR;

            % Assuming:
            % subs_tag = your string array with mixed tags
            % t = table with columns RT and I3CR (where RT contains the kind RT tags and I3CR contains replacement values)

            % Find positions where subs_tag matches any tag in t.RT
            [is_kind_RT, loc] = ismember(subs_tag, t.RT);

            % Replace only the kind RT tags with their corresponding I3CR values
            % loc(is_kind_RT) gives the indices in t where matches were found
            tag(is_kind_RT) = t.I3CR(loc(is_kind_RT));            
        end

        function [ts1,ts2,Fs] = match_time_series(this, this1, this2)
            arguments
                this mlraut.Lee2024
                this1 mlraut.AnalyticSignalHCP
                this2 mlraut.AnalyticSignalHCP
            end

            physio1 = mean(this1.physio_signal, 2);
            physio2 = mean(this2.physio_signal, 2);
            tr = max(this1.tr, this2.tr);
            Fs = min(this1.Fs, this2.Fs);
            Nt1 = size(physio1, 1);
            Nt2 = size(physio2, 1);
            Nt = min(Nt1, Nt2);
            t = 0:tr:tr*(Nt - 1);

            % immediately return if time series are already matched
            if this1.tr == this2.tr && Nt1 == Nt2
                ts1 = physio1;
                ts2 = physio2;
                Fs = this1.Fs;
                return
            end

            % interpolate to worst case tr and Nt
            t1 = 0:this1.tr:this1.tr*(Nt1 - 1);
            ts1 = interp1(t1, physio1, t);

            t2 = 0:this2.tr:this2.tr*(Nt2 - 1);
            ts2 = interp1(t2, physio2, t);
        end

        function psi = measure_uthr_avgxyzt(~, mat_fn, uthr_fn, meas)
            if ~isfile(mat_fn) || ~contains(mat_fn, "concat")
                psi = nan;
                return
            end
            ld = load(mat_fn);
            psi = ld.this.(meas)(ld.this.bold_signal, ld.this.physio_signal);  % Nt x Nx
            psi = squeeze(psi);

            if ~isfile(uthr_fn) || ~contains(uthr_fn, "uthr")
                psi = nan;
                return
            end
            cifti = cifti_read(convertStringsToChars(uthr_fn));
            select = logical(cifti.cdata');  % 1 x Nx
            if isvector(psi)
                psi = asrow(psi);
                psi = mean(psi(select), "all");
            else
                psi = mean(psi(:, select), "all");
            end
        end

        function h = plot_table(this, t, opts)
            arguments
                this mlraut.Lee2024
                t = []
                opts.var2 {mustBeTextScalar}
                opts.var1 {mustBeTextScalar} = "OS_Diagnosis"
                opts.f function_handle = @identity
                opts.measure {mustBeTextScalar} = ""
            end
            if isempty(t)
                ld = load(fullfile(this.matlabout_dir, "table_gbm_ifv.mat"));
                t = ld.t;
            end

            h = figure;
            var2 = opts.f(t.(opts.var2));
            scatter( ...
                t.(opts.var1), var2, 100, var2, 'filled', ...
                MarkerFaceAlpha=0.6, MarkerEdgeAlpha=0.6, MarkerEdgeColor=[0.1 0.1 0.1] ...
            )

            if any(var2 < 0)
                colormap(by_bc_bl)
                caxis([-max(var2), max(var2)])
            else
                colormap(magma)
                caxis([0, max(var2)])
            end
            xlabel("overall survival (days)")
            if isequal(opts.f, @abs)
                ylabel("|" + opts.measure + "(\psi_{\alpha=0.01}, \phi_{iFV})|")
            else
                ylabel(opts.measure + "(\psi_{\alpha=0.01}, \phi_{iFV})")
            end
            fontsize(scale=1.618)
        end

        function save_ciftis(this, asobj, opts)
            arguments
                this mlraut.Lee2024 %#ok<INUSA>
                asobj mlraut.AnalyticSignalHCP  % minimally must understand cifti
                opts.out_dir {mustBeFolder} = pwd
                opts.do_save_dynamic logical = false
                opts.do_save_regularized logical = true
                opts.template_cifti {mustBeFile} = cifti_read( ...
                    fullfile( ...
                    getenv("SINGULARITY_HOME"), ...
                    "AnalyticSignalGBM", "analytic_signal", "tmp", ...
                    "connectivity_sub-I3CR1488_ses-1_task-rest_run-all_desc-preproc_proc-v50-scaleiqr-iFV-brightest.dscalar.nii")) 
            end

            %% KLUDGE
            asobj.template_cifti = opts.template_cifti;
            
            asobj.out_dir = opts.out_dir;
            asobj.do_save = false;
            asobj.do_save_ciftis = true;
            asobj.do_save_dynamic = opts.do_save_dynamic;
            asobj.do_save_regularized = opts.do_save_regularized;
            meta_save(asobj);
        end

        function m = similarity_physios(this, p1, p2, opts)
            arguments
                this mlraut.Lee2024 %#ok<INUSA>
                p1 {mustBeTextScalar}
                p2 {mustBeTextScalar}
                opts.measure {mustBeTextScalar} = "connectivity"
            end

            if contains(p1, "*") && contains(p2, "*")
                g1 = mglob(p1);
                g2 = mglob(p2);
                Ng = min(length(g1), length(g2));
                p1 = [];
                p2 = [];
                for ig = 1:Ng
                    ld1 = load(g1(ig));
                    ld2 = load(g2(ig));
                    p1 = [p1; mean(ld1.this.physio_signal, 2)]; %#ok<AGROW>
                    p2 = [p2; mean(ld2.this.physio_signal, 2)]; %#ok<AGROW>
                end
                that = ld2.this;
            end

            if istext(p1) && istext(p2) && isfile(p1) && isfile(p2)
                ld1 = load(p1);
                ld2 = load(p2);
                p1 = mean(ld1.this.physio_signal, 2);
                p2 = mean(ld2.this.physio_signal, 2);
                that = ld2.this;
            end

            assert(isa(that, "mlraut.AnalyticSignal"))
            m = that.(opts.measure)(p1, p2);
        end

        function s = spectra(this, mat1, mat2, opts)
            arguments
                this mlraut.Lee2024
                mat1 {mustBeFile}
                mat2 {mustBeFile}
                opts.do_plot logical = true
                opts.Nf {mustBeScalarOrEmpty} = 256
            end

            ld1 = load(mat1);
            ld2 = load(mat2);
            [physio1,physio2,Fs] = this.match_time_series(ld1.this, ld2.this);            

            range_f = linspace(0.01, 0.1, opts.Nf);
            [Pxy, Pf] = cpsd(physio1, physio2, [], [], range_f, Fs);

            % Extract magnitude and phase
            magnitude = abs(Pxy);
            phase = angle(Pxy);

            % Plot
            if opts.do_plot
                figure;
                subplot(2,1,1);
                plot(Pf, magnitude);
                % xlim([.01, .05]);
                grid on;
                xlabel('Frequency (Hz)');
                ylabel('Magnitude');
                title('Cross-Spectral Magnitude');
                subplot(2,1,2);
                plot(Pf, phase);
                % xlim([.01, .05]);
                grid on;
                xlabel('Frequency (Hz)');
                ylabel('Phase (radians)');
                title('Cross-Spectral Phase');
            end

            [Cxy, Cf] = mscohere(physio1, physio2, [], [], range_f, Fs);

            % Plot coherence
            if opts.do_plot
                figure;
                plot(Cf, Cxy);
                % xlim([.01, .05]);
                grid on;
                xlabel('Frequency (Hz)');
                ylabel('Magnitude-Squared Coherence');
                title('Coherence Spectrum');
            end

            % return data struct
            s.Pxy = Pxy;
            s.Pf = Pf;
            s.Cxy = Cxy;
            s.Cf = Cf;
        end

        function spectra_by_os_dx(this)
            linew = 1;  % linewidth
            alpha = 0.5;  % alpha of line
            t = this.table_gbm_ce;  % N = 149
            t = sortrows(t, "OS_Diagnosis", "descend");  % survivors first
            os_dx = t.OS_Diagnosis;
            colors = flip(magma(max(t.OS_Diagnosis)));  % brighter first

            pwd0 = pushd(this.matlabout_dir);
            Nrows = size(t, 1);
            valid = true(1, Nrows);
            for idx = 1:Nrows
                i3cr = t.I3CRID(idx);
                try
                    g_ifv = mglob(sprintf("sub-%s/sub-%s*-iFV-brightest-*-concat.mat", i3cr, i3cr));
                    assert(~isempty(g_ifv));
                    g_ce = mglob(sprintf("sub-%s/sub-%s*-CE-*-concat.mat", i3cr, i3cr));
                    assert(~isempty(g_ce));
                    s_ = this.spectra(g_ifv(1), g_ce(1), do_plot=false); %#ok<AGROW>
                    s(idx).Pxy = s_.Pxy;
                    s(idx).Pf = s_.Pf;
                    s(idx).Cxy = s_.Cxy;
                    s(idx).Cf = s_.Cf;
                    s(idx).color = colors(os_dx(idx), :); %#ok<AGROW>
                catch ME
                    valid(idx) = false;
                    handwarning(ME)
                end
            end
            popd(pwd0);
            s = s(valid);
            Ns = sum(valid);

            figure;

            subplot(3,1,1);
            hold on;
            for idx = 1:Ns
                plot(s(idx).Pf, abs(s(idx).Pxy), Color=[s(idx).color, alpha], LineWidth=linew);
            end
            hold off;
            grid on;
            xlabel('Frequency (Hz)');
            ylabel('Magnitude');
            title('Cross-Spectral Magnitude');

            subplot(3,1,2);
            hold on;
            for idx = 1:Ns
                plot(s(idx).Pf, angle(s(idx).Pxy), Color=[s(idx).color, alpha], LineWidth=linew);
            end
            hold off;
            grid on;
            xlabel('Frequency (Hz)');
            ylabel('Phase (radians)');
            title('Cross-Spectral Phase');

            subplot(3,1, 3);
            hold on;
            for idx = 1:Ns
                plot(s(idx).Cf, s(idx).Cxy, Color=[s(idx).color, alpha], LineWidth=linew);
            end
            hold off;
            grid on;
            xlabel('Frequency (Hz)');
            ylabel('Magnitude-Squared Coherence');
            title('Coherence Spectrum');
        end

        function s = subs(this)
            s = asrow(mglob(fullfile(this.matlabout_dir, "sub-*")));
        end

        function t = table_found_RT_I3CR(this)
            ld = load(fullfile(this.gbm_datashare_dir, "found_RT_I3CR.mat"));
            t = ld.found_RT_I3CR;
        end

        function t = table_gbm_ce(this)
            if ~isempty(this.table_gbm_ce_)
                t = this.table_gbm_ce_;
                return
            end

            fqfn = fullfile(this.gbm_datashare_dir, "table_gbm_ce.mat");
            if isfile(fqfn)
                ld = load(fqfn);
                t = ld.t;
                return
            end

            % filter sub with physio CE
            t = this.table_gbm_datashare();
            filter = ismember(t.I3CRID, extractAfter(mybasename(this.unique_subs_with_ce), "-"));
            t = t(filter, :);
            Nrows = size(t, 1);

            % new variables for table
            angle = nan(Nrows, 1);
            correlation = nan(Nrows, 1);
            T = nan(Nrows, 1);
            X = nan(Nrows, 1);
            Y = nan(Nrows, 1);
            Z = nan(Nrows, 1);

            % build new variables for table
            for row = 1:Nrows
                try
                    pwd_ = pushd(fullfile(this.matlabout_dir, "sub-"+t.I3CRID(row)));
                    angle(row) = mean(this.similarity_physios( ...
                        "sub-*-CE*.mat", "sub-*-iFV-brightest*.mat", ...
                        measure="angle" ...
                    ));
                    correlation(row) = this.similarity_physios( ...
                        "sub-*-CE*.mat", "sub-*-iFV-brightest*.mat", ...
                        measure="connectivity" ...
                    );
                    T(row) = mean(this.similarity_physios( ...
                        "sub-*-CE*.mat", "sub-*-iFV-brightest*.mat", ...
                        measure="T" ...
                    ));
                    X(row) = mean(this.similarity_physios( ...
                        "sub-*-CE*.mat", "sub-*-iFV-brightest*.mat", ...
                        measure="X" ...
                    ));
                    Y(row) = mean(this.similarity_physios( ...
                        "sub-*-CE*.mat", "sub-*-iFV-brightest*.mat", ...
                        measure="Y" ...
                    ));
                    Z(row) = mean(this.similarity_physios( ...
                        "sub-*-CE*.mat", "sub-*-iFV-brightest*.mat", ...
                        measure="Z" ...
                    ));
                    popd(pwd_);
                catch ME
                    handwarning(ME);
                end
            end
            t.angle = angle;
            t.correlation = correlation;
            t.T = T;
            t.X = X;
            t.Y = Y;
            t.Z = Z;

            % cache table
            this.table_gbm_ce_ = t;
        end

        function t = table_gbm_ifv(this)
            if ~isempty(this.table_gbm_ifv_)
                t = this.table_gbm_ifv_;
                return
            end

            % filter sub with physio CE
            t = this.table_gbm_datashare();
            filter = ismember(t.I3CRID, extractAfter(mybasename(this.unique_subs), "-"));
            t = t(filter, :);
            Nrows = size(t, 1);

            % new variables for table
            angle = nan(Nrows, 1);
            correlation = nan(Nrows, 1);
            T = nan(Nrows, 1);
            X = nan(Nrows, 1);
            Y = nan(Nrows, 1);
            Z = nan(Nrows, 1);

            angle_uthr_fn = fullfile(this.matlabout_dir, ...
                "mean_angle_as_sub-all_ses-all_proc-v50-iFV-brightest-Lee2024-build-uthr-dscalar-alpha0p01.dscalar.nii");
            connex_uthr_fn = fullfile(this.matlabout_dir, ...
                "mean_connex_as_sub-all_ses-all_proc-v50-iFV-brightest-Lee2024-build-uthr-dscalar-alpha0p01.dscalar.nii");
            T_uthr_fn = fullfile(this.matlabout_dir, ...
                "mean_T_as_sub-all_ses-all_proc-v50-iFV-brightest-Lee2024-build-uthr-dscalar-alpha0p01.dscalar.nii");
            X_uthr_fn = fullfile(this.matlabout_dir, ...
                "mean_X_as_sub-all_ses-all_proc-v50-iFV-brightest-Lee2024-build-uthr-dscalar-alpha0p01.dscalar.nii");
            Y_uthr_fn = fullfile(this.matlabout_dir, ...
                "mean_Y_as_sub-all_ses-all_proc-v50-iFV-brightest-Lee2024-build-uthr-dscalar-alpha0p01.dscalar.nii");
            Z_uthr_fn = fullfile(this.matlabout_dir, ...
                "mean_Z_as_sub-all_ses-all_proc-v50-iFV-brightest-Lee2024-build-uthr-dscalar-alpha0p01.dscalar.nii");

            % build new variables for table
            for row = 1:Nrows
                try
                    fprintf("working in %s ...\n", t.I3CRID(row));
                    pwd_ = pushd(fullfile(this.matlabout_dir, "sub-"+t.I3CRID(row)));
                    globbed = mglob("sub-*iFV-brightest*concat.mat");
                    mat_fn = globbed(1);                    
                    angle(row) = this.measure_uthr_avgxyzt( ...
                        mat_fn, angle_uthr_fn, ...
                        "angle" ...
                    );
                    correlation(row) = this.measure_uthr_avgxyzt( ...
                        mat_fn, connex_uthr_fn, ...
                        "connectivity" ...
                    );
                    T(row) = this.measure_uthr_avgxyzt( ...
                        mat_fn, T_uthr_fn, ...
                        "T" ...
                    );
                    X(row) = this.measure_uthr_avgxyzt( ...
                        mat_fn, X_uthr_fn, ...
                        "X" ...
                    );
                    Y(row) = this.measure_uthr_avgxyzt( ...
                        mat_fn, Y_uthr_fn, ...
                        "Y" ...
                    );
                    Z(row) = this.measure_uthr_avgxyzt( ...
                        mat_fn, Z_uthr_fn, ...
                        "Z" ...
                    );
                    popd(pwd_);
                catch ME
                    handwarning(ME);
                end
            end
            t.angle = angle;
            t.correlation = correlation;
            t.T = T;
            t.X = X;
            t.Y = Y;
            t.Z = Z;

            % cache table
            this.table_gbm_ifv_ = t;
        end

        function t = table_gbm_datashare(this)
            ld = load(fullfile(this.gbm_datashare_dir, "GBMClinicalDatabasesorted.mat"));
            t = ld.GBMClinicalDatabasesorted;

            % rename RT089 -> I3CR0000
            row_rt089 = contains(t.I3CRID, "RT089");
            t.I3CRID(row_rt089) = "I3CR0000";
            t = sortrows(t, "I3CRID");
        end

        function s = unique_subs(this)
            % returns ~ "/Users/jjlee/Singularity/AnalyticSignalGBM/analytic_signal/matlabout/sub-I3CR0000"
            % mybasename(unique_subs()) ~ "sub-I3CR0000"

            s = asrow(mglob(fullfile(this.matlabout_dir, "sub-*", "*.mat")));
            s = fileparts(s);
            s = unique(s);
        end

        function s = unique_subs_with_ce(this)
            % returns ~ "/Users/jjlee/Singularity/AnalyticSignalGBM/analytic_signal/matlabout/sub-I3CR0000"
            % mybasename(unique_subs()) ~ "sub-I3CR0000"

            s = asrow(mglob(fullfile(this.matlabout_dir, "sub-*", "*CE*.mat")));
            s = fileparts(s);
            s = unique(s);
        end

        function view_physio_roi(~, sub)
            %% start in ciftify folder
            
            globbed = mglob( ...
                fullfile(sub, "MNINonLinear", "Results", "ses-1_task-rest_run-0*", "ses*.nii.gz"));
            globbed = globbed(~contains(globbed, "_avgt"));
            for g = globbed
                fprintf("%s: viewing %s\n", stackstr(), g);
                bold = mlfourd.ImagingContext2(g);
                wmp = mlsurfer.Wmparc( ...
                    fullfile(sub, "MNINonLinear", "wmparc.nii.gz"));
                roi = wmp.select_4th_ventricle + wmp.select_3rd_ventricle;
                bold.view_qc([roi, wmp.T1], cache_memory=false, options='');
            end
        end

        function write_mean_ciftis(this, ld, data, field, opts)
            arguments
                this mlraut.Lee2024
                ld struct
                data cell
                field {mustBeTextScalar}
                opts.physio {mustBeTextScalar} = "iFV-brightest"
                opts.rsn {mustBeScalarOrEmpty} = -1
            end
            valid_rsn = 1 <= opts.rsn && opts.rsn <= 9;
            rsn_tag = "rsn" + opts.rsn;

            Nd = numel(data);
            cdata = [];
            for id = 1:Nd                
                cdata = [cdata; data{id}.(field)]; %#ok<AGROW>
            end
            cdata = mean(cdata, 1);  % average over subs

            caller = strrep(stackstr(3, use_dashes=true), "Lee2024-", "");
            if valid_rsn && (contains(field, "X") || contains(field, "Y"))
                fn = fullfile( ...
                    this.matlabout_dir, ...
                    "mean_"+field+"_as_sub-all_ses-all_proc-v50-"+opts.physio+"-Lee2024-"+caller+"-"+rsn_tag+".dscalar.nii");
            else
                fn = fullfile( ...
                    this.matlabout_dir, ...
                    "mean_"+field+"_as_sub-all_ses-all_proc-v50-"+opts.physio+"-Lee2024-"+caller+".dscalar.nii");
            end
            ld.this.write_cifti(cdata, fn);
        end

        function write_metric_stats(this, metric_lbl)
            arguments
                this mlraut.Lee2024 %#ok<INUSA>
                metric_lbl {mustBeTextScalar} = "abs_as_"
            end

            out_dir = '/home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalHCPAging';
            %out_dir = '/home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalHCP';
            cd(out_dir);

            g = glob(fullfile(out_dir, 'HCA*', sprintf('%s*_avgt.dscalar.nii', metric_lbl)));
            %g = glob(fullfile(out_dir, '*', 'AnalyticSignalPar_call_subject_as_proc-norm_xyzt-ROI_*.mat'));
            leng = length(g);
            metric_mu = zeros(91282, 1);
            metric_sig2 = zeros(91282, 1);

            for gidx = 1:leng
                try
                    % tic
                    the_cifti = cifti_read(g{gidx});
                    metric_mu = metric_mu + the_cifti.cdata/leng; % mean
                    % toc % 0.8 - 1.6 sec
                catch ME
                    handwarning(ME)
                end
            end
            save(metric_lbl+"_mu.mat", "metric_mu");
            the_cifti.cdata = metric_mu;
            cifti_write(the_cifti, convertStringsToChars(metric_lbl+"_mu.dscalar.nii"));
        
            for gidx = 1:leng
                try
                    the_cifti = cifti_read(g{gidx});
                    metric_sig2 = metric_sig2 + (the_cifti.cdata - metric_mu).^2/leng; % var
                catch ME
                    handwarning(ME)
                end
            end
            metric_sig = metric_sig2.^(0.5);
            save(metric_lbl+"_sig.mat", "metric_sig");
            the_cifti.cdata = metric_sig;
            cifti_write(the_cifti, convertStringsToChars(metric_lbl+"_sig.dscalar.nii"));

            % coeff. of var.
            cov = the_cifti;
            cov.cdata  = metric_sig ./ metric_mu;
            cifti_write(cov, convertStringsToChars(metric_lbl+"_cov.dscalar.nii"));

            % snr
            snr = the_cifti;
            snr.cdata  = metric_mu ./ metric_sig;
            cifti_write(snr, convertStringsToChars(metric_lbl+"_snr.dscalar.nii"));
        end
            
        function write_var_ciftis(this, ld, data, field, opts)
            arguments
                this mlraut.Lee2024
                ld struct
                data cell
                field {mustBeTextScalar}
                opts.physio {mustBeTextScalar} = "iFV-brightest"
                opts.rsn {mustBeScalarOrEmpty} = -1
            end
            valid_rsn = 1 <= opts.rsn && opts.rsn <= 9;
            rsn_tag = "rsn" + opts.rsn;

            Nd = numel(data);
            cdata = [];
            for id = 1:Nd                
                cdata = [cdata; data{id}.(field)]; %#ok<AGROW>
            end
            cdata = var(cdata, 0, 1);  % norm by Nsubs - 1; var over subs

            caller = strrep(stackstr(3, use_dashes=true), "Lee2024-", "");
            if valid_rsn && (contains(field, "X") || contains(field, "Y"))
                fn = fullfile( ...
                    this.matlabout_dir, ...
                    "var_"+field+"_as_sub-all_ses-all_proc-v50-"+opts.physio+"-Lee2024-"+caller+"-"+rsn_tag+".dscalar.nii");
            else
                fn = fullfile( ...
                    this.matlabout_dir, ...
                    "var_"+field+"_as_sub-all_ses-all_proc-v50-"+opts.physio+"-Lee2024-"+caller+".dscalar.nii");
            end
            ld.this.write_cifti(cdata, fn);
        end
    end

    %% PROTECTED
    
    properties (Access = protected)
        table_gbm_ce_
        table_gbm_ifv_
    end

    methods (Access = protected)
        function [datum,ld] = gather_signals(this, subdir, opts)
            %% excludes *-concat.mat

            arguments
                this mlraut.Lee2024
                subdir {mustBeFolder}
                opts.physio {mustBeTextScalar} = "iFV-brightest"
                opts.rsn {mustBeInteger} = -1
                opts.weight {mustBeScalarOrEmpty} = 1
                opts.build_residuals logical = false
            end
            use_rsn = 1 <= opts.rsn && opts.rsn <= 9;

            m = mglob(fullfile(subdir, "*run-all*"+opts.physio+"*.mat"));
            m = m(~contains(m, "-concat.mat"));
            if isempty(m)
                error("mlraut:RuntimeError", stackstr(use_dashes=true) + "-data-missing")
            end
            psi = [];
            phi = [];
            ctx_psi = [];
            ctx_phi = [];
            for m1 = asrow(m)
                ld = load(m1);
                psi = [psi; ld.this.bold_signal]; %#ok<AGROW>
                phi = [phi; ld.this.physio_signal]; %#ok<AGROW>
                if use_rsn
                    ctx_psi = [ctx_psi; ld.this.HCP_signals.ctx.psi]; %#ok<AGROW>
                    ctx_phi = [ctx_phi; ld.this.HCP_signals.ctx.phi]; %#ok<AGROW>
                end
            end
            Nt = ld.this.num_frames;
            twist = mlraut.Twistors(ld.this);
            connex = asrow(ld.this.connectivity(psi, mean(phi, 2)));
            angle = mean(twist.angle(psi, phi), 1);  % average over t

            % cortical X(psi, phi) >= 0, region ~ opts.rsn ~ DMN, biased but informative
            if use_rsn
                angle_rsn = twist.angle(ctx_psi(:, opts.rsn), ctx_phi(:, opts.rsn));
                t_interesting = cos(angle_rsn) > 0;
                u_interesting = sin(angle_rsn) > 0;
            else
                t_interesting = true(Nt, 1);
                u_interesting = true(Nt, 1);
            end

            % mean statistic
            T = mean(twist.T(psi, phi), 1);
            X = twist.X(psi, phi);
            X = mean(X(t_interesting, :), 1);
            Y = twist.Y(psi, phi);
            Y = mean(Y(u_interesting, :), 1);
            Z = mean(twist.Z(psi, phi), 1);

            if opts.build_residuals

                % reference from 725 HCP Aging subs
                c = cifti_read(this.hcpaging_connex);
                connex_ref = c.cdata';
                c = cifti_read(this.hcpaging_angle);
                angle_ref = c.cdata';
                c = cifti_read(this.hcpaging_T);
                T_ref = c.cdata';
                c = cifti_read(this.hcpaging_X);
                X_ref = c.cdata';
                c = cifti_read(this.hcpaging_Y);
                Y_ref = c.cdata';
                c = cifti_read(this.hcpaging_Z);
                Z_ref = c.cdata';

                datum = struct( ...
                    "connex", opts.weight * (connex - connex_ref), ...
                    "angle", opts.weight * (angle - angle_ref), ...
                    "T", opts.weight * (T - T_ref), ...
                    "X", opts.weight * (X - X_ref), ...
                    "Y", opts.weight * (Y - Y_ref), ...
                    "Z", opts.weight * (Z - Z_ref));
            else
                datum = struct( ...
                    "connex", opts.weight * connex , ...
                    "angle", opts.weight * angle, ...
                    "T", opts.weight * T, ...
                    "X", opts.weight * X, ...
                    "Y", opts.weight * Y, ...
                    "Z", opts.weight * Z);
            end
        end

        function data = reweight_data(~, data, weight)
            for id = 1:numel(data)
                datum = data{id};
                datum.connex = weight * datum.connex;
                datum.angle = weight * datum.angle;
                datum.T = weight * datum.T;
                datum.X = weight * datum.X;
                datum.Y = weight * datum.Y;
                datum.Z = weight * datum.Z;
                data{id} = datum;
            end
        end
    end

    methods (Static, Access = private)
        function obj = identity(obj)
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
