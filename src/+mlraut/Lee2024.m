classdef Lee2024 < handle
    %% line1
    %  line2
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
                "mean_comparator_as_sub-all_ses-all_proc-v50-iFV-mean-twistor-nlim725_avgt.dscalar.nii");
            this.hcpaging_angle = fullfile(  ...
                this.analytic_signal_hcpaging_dir, ...
                "mean_angle_as_sub-all_ses-all_proc-v50-iFV-mean-element-nlim725_avgt.dscalar.nii");
            this.hcpaging_T = fullfile( ...
                this.analytic_signal_hcpaging_dir, ...
                "mean_T_as_sub-all_ses-all_proc-v50-iFV-mean-element-nlim725_avgt.dscalar.nii");
            this.hcpaging_X = fullfile( ...
                this.analytic_signal_hcpaging_dir, ...
                "mean_X_as_sub-all_ses-all_proc-v50-iFV-mean-twistor-rsn7-nlim725_avgt.dscalar.nii");
            this.hcpaging_Y = fullfile( ...
                this.analytic_signal_hcpaging_dir, ...
                "mean_Y_as_sub-all_ses-all_proc-v50-iFV-mean-twistor-rsn7-nlim725_avgt.dscalar.nii");
            this.hcpaging_Z = fullfile( ...
                this.analytic_signal_hcpaging_dir, ...
                "mean_Z_as_sub-all_ses-all_proc-v50-iFV-mean-element-nlim725_avgt.dscalar.nii");
        end

        function build_mean_for_gbm(this, opts)
            arguments
                this mlraut.Lee2024
                opts.physio {mustBeTextScalar} = "iFV"
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

        function s = subs(this)
            s = asrow(mglob(fullfile(this.matlabout_dir, "sub-*")));
        end

        function t = table_found_RT_I3CR(this)
            ld = load(fullfile(this.gbm_datashare_dir, "found_RT_I3CR.mat"));
            t = ld.found_RT_I3CR;
        end

        function t = table_gbm_datashare(this)
            ld = load(fullfile(this.gbm_datashare_dir, "GBMClinicalDatabasesorted.mat"));
            t = ld.GBMClinicalDatabasesorted;
        end

        function s = unique_subs(this)
            s = asrow(mglob(fullfile(this.matlabout_dir, "sub-*", "*.mat")));
            s = fileparts(s);
            s = unique(s);
        end

        function write_mean_ciftis(this, ld, data, field, opts)
            arguments
                this mlraut.Lee2024
                ld struct
                data cell
                field {mustBeTextScalar}
                opts.physio {mustBeTextScalar} = "iFV"
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
    end

    methods (Access = protected)
        function [datum,ld] = gather_signals(this, subdir, opts)
            arguments
                this mlraut.Lee2024
                subdir {mustBeFolder}
                opts.physio {mustBeTextScalar} = "iFV"
                opts.rsn {mustBeInteger} = -1
                opts.weight {mustBeScalarOrEmpty} = 1
                opts.build_residuals logical = true
            end
            use_rsn = 1 <= opts.rsn && opts.rsn <= 9;

            m = mglob(fullfile(subdir, "*"+opts.physio+"*.mat"));
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
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
