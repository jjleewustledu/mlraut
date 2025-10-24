classdef BrainSmashAdapter < handle
    %% line1
    %  line2
    %  
    %  Created 21-Oct-2025 10:13:51 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 24.2.0.2923080 (R2024b) Update 6 for MACA64.  Copyright 2025 John J. Lee.
    
    methods
        function this = BrainSmashAdapter(varargin)
        end
    end

    methods (Static)
        function [j,c] = cluster_sample_subs_tasks(globbing_mat, opts)
            arguments
                globbing_mat {mustBeText} = ...
                    fullfile( ...
                    getenv('SINGULARITY_HOME'), 'AnalyticSignalHCP', 'mlraut_AnalyticSignalHCPPar_globbing.mat')
                opts.globbing_var = "globbed"
                opts.new_physio {mustBeTextScalar} = "iFV"
                opts.anatomy {mustBeTextScalar} = "ctx"
                opts.measure {mustBeTextScalar} = "phase_locked_values"
                opts.measure_name {mustBeTextScalar} = "plvs"
                opts.out_dir {mustBeTextScalar} = "/scratch/jjlee/Singularity/AnalyticSignalHCP"
                opts.real logical = true  % real() | imag()
            end

            % additional internal params
            Ncol = 16;  % Ncol ~ num jobs; Nrow ~ ceil(1113/16) = 70 subs per job
            account_name = 'joshua_shimony';
            mempercpu = '128gb';
            walltime = '24:00:00';

            % `globbing_mat` -> `subs`
            if isscalar(globbing_mat) && endsWith(globbing_mat, ".mat") && ~isemptytext(opts.globbing_var)
                ld = load(globbing_mat);
                subs = ld.(opts.globbing_var);
            else
                subs = globbing_mat;
            end
            subs = convertCharsToStrings(subs);
            subs = asrow(subs);
            subs_nontrivial = subs;

            % pad and reshape globbed
            Nrow = ceil(numel(subs)/Ncol);
            padding = repmat("", [1, Ncol*Nrow - numel(subs)]);
            subs = [subs, padding];
            subs = reshape(subs, Nrow, Ncol);
            fprintf("%s:subs:\n", stackstr());
            fprintf("\t%s ... %s (len=%i)\n", subs_nontrivial(1), subs_nontrivial(end), length(subs_nontrivial));
            fprintf("\tstring array size: %s\n", mat2str(size(subs)));

            %% contact cluster slurm

            warning('off', 'MATLAB:legacy:batchSyntax');
            warning('off', 'parallel:convenience:BatchFunctionNestedCellArray');
            warning('off', 'MATLAB:TooManyInputs');

            c = mlraut.CHPC3.propcluster(account_name, mempercpu=mempercpu, walltime=walltime);
            subs = convertStringsToChars(subs);
            for col_idx = 1:Ncol
                try
                    j = c.batch( ...
                        @mlraut.BrainSmashAdapter.sample_subs_tasks, ...
                        1, ...
                        {subs(:, col_idx), ...
                         'globbing_var', '', ...
                         'new_physio', opts.new_physio, ...
                         'anatomy', opts.anatomy, ...
                         'measure', opts.measure, ...
                         'measure_name', opts.measure_name, ...
                         'col_idx', col_idx, ...
                         'out_dir', opts.out_dir, ...
                         'real', opts.real}, ...
                        'CurrentFolder', '/scratch/jjlee/Singularity/AnalyticSignalHCP', ...
                        'AutoAddClientPath', false);
                catch ME
                    handwarning(ME)
                end
            end

            warning('on', 'MATLAB:legacy:batchSyntax');
            warning('on', 'parallel:convenience:BatchFunctionNestedCellArray');
            warning('on', 'MATLAB:TooManyInputs');
        end

        function measure_data = sample_subs_tasks(globbing_mat, opts)
            arguments
                globbing_mat {mustBeText} = ...
                    fullfile( ...
                    getenv('SINGULARITY_HOME'), 'AnalyticSignalHCP', 'mlraut_AnalyticSignalHCPPar_globbing.mat')
                opts.globbing_var = "globbed"
                opts.new_physio {mustBeTextScalar} = "RV-std"
                opts.anatomy {mustBeTextScalar} = "ctx"
                opts.measure {mustBeTextScalar} = "phase_locked_values"
                opts.measure_name {mustBeTextScalar} = "plvs"
                opts.col_idx {mustBeScalarOrEmpty} = nan
                opts.out_dir {mustBeFolder} = "/home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalHCP"
                opts.stat = []
                opts.real logical = true  % real() | imag()
                opts.do_save logical = true
            end

            % `globbing_mat` -> `subs`
            if isscalar(globbing_mat) && endsWith(globbing_mat, ".mat") && ~isemptytext(opts.globbing_var)
                ld = load(globbing_mat);
                subs = ld.(opts.globbing_var);
            else
                subs = globbing_mat;
            end
            subs = convertCharsToStrings(subs);
            subs = asrow(subs);

            % construct `as` which supplies utilities from AnalyticSignalHCP
            warning("off", "MATLAB:class:LoadDefinitionUpdated");
            mat_fqfn = [];
            idx = 0;
            while isempty(mat_fqfn) && idx < length(subs)
                idx = idx + 1;
                patt = fullfile( ...
                    opts.out_dir, ...
                    subs(idx), ...
                    sprintf("sub-%s_ses-*_proc-*ASHCPPar*.mat", subs(idx)));
                mat_fqfn = mglob(patt);
            end
            as = mlraut.AnalyticSignalHCPPar.load(mat_fqfn(1), class="mlraut.AnalyticSignalHCPPar");
            as.out_dir = opts.out_dir;
            warning("on", "MATLAB:class:LoadDefinitionUpdated");
            
            % sample subjects and tasks per subject (typically 2-4)
            samples = [];
            for sub = asrow(subs)
                if isemptytext(sub)
                    continue
                end
                mat_fqfns = mglob( ...
                    fullfile(as.out_dir, sub, sprintf("sub-%s_ses-*ASHCPPar*.mat", sub)));
                if any(contains(mat_fqfns, "-all"))
                    mat_fqfns = mat_fqfns(~contains(mat_fqfns, "-all"));  % don't double count
                end
                if isemptytext(mat_fqfns)
                    continue
                end

                img_ = [];
                nerror = 0;
                for mat_fqfn = mat_fqfns
                    tic
                    try
                        ld = load(mat_fqfn);
                        psi = ld.this_subset.bold_signal;
                        if isemptytext(opts.new_physio)
                            phi = ld.this_subset.physio_signal;
                        elseif strcmpi(opts.new_physio, ld.this_subset.source_physio)
                            phi = ld.this_subset.physio_signal;
                        elseif strcmpi(opts.new_physio, 'vis')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 1);
                        elseif strcmpi(opts.new_physio, 'sms')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 2);
                        elseif strcmpi(opts.new_physio, 'dan')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 3);
                        elseif strcmpi(opts.new_physio, 'van')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 4);
                        elseif strcmpi(opts.new_physio, 'lim')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 5);
                        elseif strcmpi(opts.new_physio, 'fpn')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 6);
                        elseif strcmpi(opts.new_physio, 'dmn')
                            phi = ld.this_subset.HCP_signals.(opts.anatomy).psi(:, 7);
                        else
                            re_phi = ld.this_subset.physio_supplementary(opts.new_physio);
                            assert(~isempty(re_phi))
                            assert(all(isfinite(re_phi)))
                            phi = hilbert(re_phi);
                        end

                        img_ = [img_; as.(opts.measure)(psi, phi)]; %#ok<AGROW>
                    catch ME
                        fprintf("while working with %s: ", mat_fqfn);
                        nerror = nerror + 1;
                        handwarning(ME)
                    end
                    fprintf("time working with %s: ", mat_fqfn);
                    toc
                end

                img_ = mean(img_, 1, "omitnan");  % 1 x `ngo`
                samples = [samples; img_];  %#ok<AGROW> % accum N_{subs} x `ngo`
            end

            % build parts of filenames
            ntag = size(samples, 1);
            if ~strcmp(opts.anatomy, "ctx")
                ptags = strrep(as.tags, as.source_physio, sprintf("%s-%s", opts.anatomy, opts.new_physio));
            else
                ptags = strrep(as.tags, as.source_physio, opts.new_physio);
            end

            % assemble struct of measure
            measure_data = struct( ...
                "measure", opts.measure, ...
                "measure_name", opts.measure_name, ...
                "as", as, "subs", subs, "samples", samples, "ntag", ntag, "ptags", ptags, "nerror", nerror);

            % save intermediates to out_fqfn
            if opts.do_save
                out_fqfn = fullfile( ...
                    opts.out_dir, ...
                    sprintf('%s_%s_as_sub-%i_ses-%i_%s_par%i.mat', ...
                    stackstr(use_dashes=false), opts.measure_name, ntag, ntag, ptags, opts.col_idx));
                save(out_fqfn, "measure_data", "-v7.3");
            end

            % apply stat function handle, then save out_fqfn as cifti, then return            
            if isa(opts.stat, "function_handle")
                stat_samples_ = opts.stat(samples, 1);
                if opts.real
                    stat_samples_ = real(stat_samples_);
                    prefix = "re";
                else
                    stat_samples_ = imag(stat_samples_);
                    prefix = "im";
                end
                measure_data_ = measure_data;  % expensive copy; set opts.stat <- [] if unnecessary
                measure_data_.samples = stat_samples_;
                mlraut.BrainSmashAdapter.write_cifti( ...
                    measure_data_, ...
                    prefix=prefix+func2str(opts.stat), ...
                    out_dir=opts.out_dir);
            end
        end

        function out_fqfn = write_cifti(measure_data, opts)
            arguments
                measure_data struct {mustBeNonempty}
                opts.prefix {mustBeTextScalar} = ""
                opts.out_dir {mustBeFolder} = "/home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalHCP"
            end

            measure_name = measure_data.measure_name;
            ntag = measure_data.ntag;  % abbrevs
            ptags = measure_data.ptags;
            out_fqfn =  fullfile(opts.out_dir, ...
                sprintf('%s%s_as_sub-%i_ses-%i_%s', ...
                opts.prefix, measure_name, ntag, ntag, ptags));
            measure_data.as.cifti.write_cifti(measure_data.samples, out_fqfn);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
