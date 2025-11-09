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
        function [j,c] = cluster_test()

            % N.B.:  Matlab requires ~ 1101020 kB
            c = mlraut.CHPC3.propcluster_free("joshua_shimony", mempercpu="10gb", walltime="00:10:00");
            try
                j = c.batch( ...
                    @mlraut.CHPC3.parallel_example, ...
                    1, ...
                    {8}, ...
                    'CurrentFolder', '/scratch/jjlee/tmp', ...
                    'AutoAddClientPath', false);
            catch ME
                handwarning(ME)
            end
        end

        function [j,c] = cluster_sample_subs_tasks(globbing_mat, opts)
            arguments
                globbing_mat {mustBeText} = ...
                    fullfile( ...
                    getenv('SINGULARITY_HOME'), 'AnalyticSignalHCP', 'mlraut_AnalyticSignalHCPPar_verified_globbed.mat')
                opts.globbing_var = "verified_globbed"
                opts.new_physio {mustBeTextScalar} = "iFV"
                opts.anatomy {mustBeTextScalar} = "ctx"
                opts.measure {mustBeTextScalar} = "phase_locked_values"  % "Z_sup_0_tuned"
                opts.measure_name {mustBeTextScalar} = "plvs"  % "Zsup0tuned"
                opts.out_dir {mustBeTextScalar} = "/scratch/jjlee/Singularity/AnalyticSignalHCP/cache"
                opts.real logical = true  % real() | imag()
                opts.funch function_handle = @mlraut.BrainSmashAdapter.sample_subs_tasks_scrambled
                % opts.funch function_handle = @mlraut.BrainSmashAdapter.sample_subs_tasks
                opts.mempercpu char {mustBeTextScalar} = '32gb'  % '32gb'
                opts.walltime char {mustBeTextScalar} = '24:00:00'
                %opts.walltime char {mustBeTextScalar} = '8:00:00'
                opts.starting_index_col {mustBeInteger} = 1
                opts.Ncol {mustBeInteger} = 86
                opts.partition_free logical = false
            end

            % additional internal params
            Ncol = opts.Ncol;  % Ncol ~ num jobs; Nrow ~ ceil(1097/16) = 69 subs per job
            account_name = 'joshua_shimony';

            % `globbing_mat` -> `subs`
            if isscalar(string(globbing_mat)) && endsWith(globbing_mat, ".mat") && ~isemptytext(opts.globbing_var)
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
            if Ncol > 1
                padding = repmat("", [1, Ncol*Nrow - numel(subs)]);
                subs = [subs, padding];                
                subs = reshape(subs, Nrow, Ncol);
            end
            fprintf("%s:subs:\n", stackstr());
            fprintf("\t%s ... %s (len=%i)\n", subs_nontrivial(1), subs_nontrivial(end), length(subs_nontrivial));
            fprintf("\tstring array size: %s\n", mat2str(size(subs)));

            %% contact cluster slurm

            warning('off', 'MATLAB:legacy:batchSyntax');
            warning('off', 'parallel:convenience:BatchFunctionNestedCellArray');
            warning('off', 'MATLAB:TooManyInputs');

            if opts.partition_free
                opts.walltime = '3:00:00';
                c = mlraut.CHPC3.propcluster_free(account_name, mempercpu=opts.mempercpu, walltime=opts.walltime);
            else
                c = mlraut.CHPC3.propcluster(account_name, mempercpu=opts.mempercpu, walltime=opts.walltime);
            end
            subs = ensureCell(convertStringsToChars(subs));            
            for col_idx = opts.starting_index_col:opts.starting_index_col+Ncol-1
                try
                    j = c.batch( ...
                        opts.funch, ...
                        1, ...
                        {subs(:, col_idx-opts.starting_index_col+1), ...
                         'globbing_var', '', ...
                         'new_physio', opts.new_physio, ...
                         'anatomy', opts.anatomy, ...
                         'measure', opts.measure, ...
                         'measure_name', opts.measure_name, ...
                         'col_idx', col_idx, ...
                         'out_dir', opts.out_dir, ...
                         'real', opts.real}, ...
                        'CurrentFolder', opts.out_dir, ...
                        'AutoAddClientPath', false);
                catch ME
                    handwarning(ME)
                end
            end

            warning('on', 'MATLAB:legacy:batchSyntax');
            warning('on', 'parallel:convenience:BatchFunctionNestedCellArray');
            warning('on', 'MATLAB:TooManyInputs');
        end

        function globbing_mat = sample_subs_tasks(globbing_mat, opts)
            %% collects opts.measure for all tasks for all sub
            arguments
                globbing_mat {mustBeText} = ...
                    fullfile( ...
                    getenv('SINGULARITY_HOME'), 'AnalyticSignalHCP', 'mlraut_AnalyticSignalHCPPar_globbing.mat')
                opts.globbing_var = "globbed"
                opts.new_physio {mustBeTextScalar} = "RV-std"
                opts.anatomy {mustBeTextScalar} = "ctx"
                opts.measure {mustBeTextScalar} = "Z_sup_2_tuned"
                opts.measure_name {mustBeTextScalar} = "Zsup2tuned"
                opts.col_idx {mustBeScalarOrEmpty} = nan
                opts.out_dir {mustBeFolder} = "/home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalHCP"
                opts.stat = []  % e.g., @identity, @mean; [] skips writing cifti
                opts.real logical = true  % real() | imag()
                opts.do_save logical = true
            end

            % `globbing_mat` -> `subs`
            if isscalar(string(globbing_mat)) && endsWith(globbing_mat, ".mat") && ~isemptytext(opts.globbing_var)
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
                    fullfile(opts.out_dir, sub, sprintf("sub-%s_ses-*ASHCPPar*.mat", sub)));
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

                        img_ = [img_; as.(opts.measure)(psi, phi)]; %#ok<AGROW> % accum N_{angles} x `ngo`
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
                "as", [], "subs", subs, "samples", samples, "ntag", ntag, "ptags", ptags, "nerror", nerror);

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
                try
                    stat_samples_ = opts.stat(samples, 1);
                catch
                    try
                        stat_samples_ = opts.stat(samples, [], 1);
                    catch
                        stat_samples_ = opts.stat(samples);
                    end
                end
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

        function globbing_mat = sample_subs_tasks_physio(globbing_mat, opts)
            %% collects physio time-seres for all tasks for all sub

            arguments
                globbing_mat {mustBeText} = ...
                    fullfile( ...
                    getenv('SINGULARITY_HOME'), 'AnalyticSignalHCP', 'mlraut_AnalyticSignalHCPPar_globbing.mat')
                opts.globbing_var = "globbed"
                opts.new_physio {mustBeTextScalar} = "RV-std"
                opts.anatomy {mustBeTextScalar} = "ctx"
                opts.measure {mustBeTextScalar} = ""
                opts.measure_name {mustBeTextScalar} = ""
                opts.col_idx {mustBeScalarOrEmpty} = nan
                opts.out_dir {mustBeFolder} = "/home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalHCP"
                opts.stat = []
                opts.real logical = true  % real() | imag()
                opts.do_save logical = true
            end

            % `globbing_mat` -> `subs`
            if isscalar(string(globbing_mat)) && endsWith(globbing_mat, ".mat") && ~isemptytext(opts.globbing_var)
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
                    fullfile(opts.out_dir, sub, sprintf("sub-%s_ses-*ASHCPPar*.mat", sub)));
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
                            phi = re_phi;
                        end

                        img_ = [img_, real(phi)]; %#ok<AGROW>  accum N_t x N_{tasks}
                    catch ME
                        fprintf("while working with %s: ", mat_fqfn);
                        nerror = nerror + 1;
                        handwarning(ME)
                    end
                    fprintf("time working with %s: ", mat_fqfn);
                    toc
                end

                samples = [samples, img_];  %#ok<AGROW> % accum N_t x (N_{tasks} + N_{subs})
            end

            % build parts of filenames
            ntag = size(samples, 2);
            if ~strcmp(opts.anatomy, "ctx")
                ptags = strrep(as.tags, as.source_physio, sprintf("%s-%s", opts.anatomy, opts.new_physio));
            else
                ptags = strrep(as.tags, as.source_physio, opts.new_physio);
            end

            % assemble struct of measure
            measure_data = struct( ...
                "new_physio", opts.new_physio, ...
                "subs", subs, "samples", samples, "ntag", ntag, "ptags", ptags, "nerror", nerror);

            % save intermediates to out_fqfn
            if opts.do_save
                out_fqfn = fullfile( ...
                    opts.out_dir, ...
                    sprintf('%s_as_sub-%i_ses-%i_%s_par%i.mat', ...
                    stackstr(use_dashes=false), ntag, ntag, ptags, opts.col_idx));
                save(out_fqfn, "measure_data", "-v7.3");
            end
        end
        
        function globbing_mat = sample_subs_tasks_scrambled(globbing_mat, opts)
            %% collects opts.measure for all tasks for all sub
            arguments
                globbing_mat {mustBeText} = "200917"
                opts.globbing_var = "verified_globbed_251to500"
                opts.new_physio {mustBeTextScalar} = "RV-std"
                opts.anatomy {mustBeTextScalar} = "ctx"
                opts.measure {mustBeTextScalar} = "phase_locked_values"
                opts.measure_name {mustBeTextScalar} = "plvs"
                opts.col_idx {mustBeScalarOrEmpty} = nan
                opts.out_dir {mustBeFolder} = "/home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalHCP/cache"
                opts.sub_folders {mustBeTextScalar} = fullfile("BrainSmashAdapter", "sample_subs_tasks_scrambled")
                opts.stat = []  % e.g., "cifti", @identity, @mean; [] skips writing cifti
                opts.real logical = true  % real() | imag()
                opts.num_scrambles {mustBeInteger} = 500
            end

            % `globbing_mat` -> `subs`
            if isscalar(string(globbing_mat)) && endsWith(globbing_mat, ".mat") && ~isemptytext(opts.globbing_var)
                ld = load(globbing_mat);
                subs = ld.(opts.globbing_var);
            else
                subs = globbing_mat;
            end
            subs = convertCharsToStrings(subs);
            subs = asrow(subs);
            fprintf("%s: ", stackstr()); disp(subs);

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
            fprintf("%s: loading ASHCPPar object from %s\n", stackstr(), mat_fqfn(1))
            as = mlraut.AnalyticSignalHCPPar.load(mat_fqfn(1), class="mlraut.AnalyticSignalHCPPar");
            as.out_dir = opts.out_dir;
            warning("on", "MATLAB:class:LoadDefinitionUpdated");

            % prepare theta_hat:  random phases ~ $[-\pi, \pi]$
            rand_lower = single(-pi);
            rand_range = single(2*pi);
            
            % sample subjects and tasks per subject (typically 2-4)
            for sub = asrow(subs)
                if isemptytext(sub)
                    continue
                end
                mat_fqfns = mglob( ...
                    fullfile(opts.out_dir, sub, sprintf("sub-%s_ses-*ASHCPPar*.mat", sub)));
                if any(contains(mat_fqfns, "-all"))
                    mat_fqfns = mat_fqfns(~contains(mat_fqfns, "-all"));  % don't double count
                end
                if isemptytext(mat_fqfns)
                    continue
                end

                img_ = nan([length(mat_fqfn), as.num_nodes, opts.num_scrambles]);  % N_task x N_go x N_scrambles
                nerror = 0;
                task_idx = 0;
                for mat_fqfn = mat_fqfns
                    tic
                    try
                        task_idx = task_idx + 1;
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

                        %% scrambling phase of phi with theta_hat; must be done at lowest level

                        scrambles_ = nan([size(psi, 2), opts.num_scrambles]);  % N_go x N_scrambles
                        for sidx = 1:opts.num_scrambles
                            theta_hat = rand_lower + rand_range * rand(size(phi), 'single');
                            theta_hat(1) = 0;  % don't phase scramble DC power
                            fft_phi_scrambled = abs(fft(phi)).*exp(1i*theta_hat);
                            phi_scrambled = ifft(fft_phi_scrambled);
                            measure_scrambled = as.(opts.measure)(psi, phi_scrambled);  % N_theta x N_go
                            scrambles_(:,sidx) = mean(measure_scrambled, 1);  % N_go x N_scrambles
                        end

                        img_(task_idx, :, :) = scrambles_;  % N_tasks x N_go x N_scrambles
                    catch ME
                        fprintf("while working with %s: ", mat_fqfn);
                        nerror = nerror + 1;
                        handwarning(ME)
                    end
                    fprintf("time working with %s: ", mat_fqfn);
                    toc
                end

                % average tasks
                samples = squeeze(mean(img_, 1, "omitnan"));  % N_go x N_scrambles
                samples = samples';  % N_scrambles x N_go

                % build parts of filenames
                mat_fqfp = myfileprefix(mat_fqfn(1));
                if ~strcmp(opts.anatomy, "ctx")
                    py_fqfp = strrep(mat_fqfp, as.source_physio, sprintf("%s-%s", opts.anatomy, opts.new_physio));
                else
                    py_fqfp = strrep(mat_fqfp, as.source_physio, opts.new_physio);
                end
                [pth_,fp_] = myfileparts(py_fqfp);
                pth_ = fullfile(myfileparts(pth_), opts.sub_folders, sub);
                ensuredir(pth_);
                py_fqfp = fullfile(pth_, opts.measure_name + "_as_" + regexprep(fp_, "_ses-\S+_", "_ses-all_"));
                py_fqfn = py_fqfp + ".npy";

                % save .npy, ensuring MATLAB is using local conda python
                nilpy = '/home/usr/jjlee/.conda/envs/matlab_py310/bin/python';
                chpcpy = fullfile(getenv('HOME'), 'miniconda3/envs/py3/bin/python');
                if isfile(nilpy)
                    pe = pyenv('Version', nilpy);
                elseif isfile(chpcpy)
                    pe = pyenv('Version', chpcpy);
                else
                    error("mlraut:RunTimeError", "%s: no python found", stackstr());
                end
                fprintf("%s: %s, %s\n", stackstr(), pe.Version, pe.Executable)
                if ~all(isnan(samples), "all")
                    py.numpy.save(py_fqfn, samples);  % N_scrambles x N_go
                end

                % apply stat function handle, then save out_fqfn as cifti, then return;
                % not for production runs
                if isempty(opts.stat)
                     continue
                end
                if isa(opts.stat, "function_handle") || contains(opts.stat, "cifti", IgnoreCase=true)
                    if isa(opts.stat, "function_handle")
                        try
                            stat_samples_ = opts.stat(samples, 1);
                        catch
                            try
                                stat_samples_ = opts.stat(samples, [], 1);
                            catch
                                stat_samples_ = opts.stat(samples);
                            end
                        end
                    else
                        stat_samples_ = samples;
                    end
                    if opts.real
                        stat_samples_ = real(stat_samples_);
                        fp = "re" + mybasename(py_fqfp);
                    else
                        stat_samples_ = imag(stat_samples_);
                        fp = "im" + mybasename(py_fqfp);
                    end
                    cii = cifti_read(fullfile(opts.out_dir, "template.dscalar.nii"));
                    if size(stat_samples_, 1) > 1
                        N_scrambles = size(stat_samples_, 1);
                        stat_samples_ = stat_samples_(1:min(100, N_scrambles), :);
                        cii.cdata = stat_samples_';  % N_go x 100 scrambles
                        cii.diminfo{2} = cifti_diminfo_make_series(min(size(stat_samples_)), 0, 1, 'SECOND');
                        ext = ".dtseries.nii";
                    else
                        cii.cdata = stat_samples_';  % N_go x 1
                        ext = ".dscalar.nii";
                    end
                    cii_fqfn =  fullfile(myfileparts(py_fqfp), fp + ext);
                    cifti_write(cii, char(cii_fqfn));
                end
            end
        end

        function mats_par_to_group_dscalar(mats_patt, opts)
            %% intended for, e.g.,
            %  fullfile("/Volumes/PrecunealSSD2/AnalyticSignalHCP", ...
            %      "BrainSmashAdapter_sample_subs_tasks_plvs_as_sub-56_ses-56_proc-RV-std-gsr0-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar-20251022225151_par16.mat")

            arguments
                mats_patt {mustBeTextScalar} = "BrainSmashAdapter_sample_subs_tasks_plvs_as_*_proc-RV-std*_par*.mat"
                opts.mat_varname {mustBeText} = ["measure_data", "samples"]
                opts.dscalar_fqfn {mustBeTextScalar} = ""
                opts.out_dir {mustBeTextScalar} = "brainsmash_adapter_groupwise"
            end
            assert(~contains(mats_patt, "tasks_physio"))
            if ~isfolder(opts.out_dir)
                mkdir(opts.out_dir)
            end
            
            % aggregate
            mats = mglob(mats_patt);
            samples_ = [];
            ntag_ = 0;
            nmats_ = 0;
            nsubs_ = 0;
            ntasks_ = 0;            
            for m = mats
                try
                    ld = load(m);
                    assert(isfield(ld, "measure_data"))                    
                    md = ld.measure_data;
                    assert(contains(m, md.measure_name))
                    samples_ = [samples_; md.samples]; %#ok<AGROW>
                    ntag_ = ntag_ + md.ntag;
                    nmats_ = nmats_ + 1;
                    nsubs_ = nsubs_ + length(md.subs);
                    ntasks_ = ntasks_ + size(md.samples, 1);
                catch ME
                    handwarning(ME);
                end
            end
            samples_ = mean(samples_, 1, "omitnan");
            samples_ = single(samples_);
            samples_ = real(samples_);
            samples_ = samples_';
            fprintf("%s:\n", stackstr());
            fprintf("\tntag: %g, nmats: %g, nsubs: %g, ntasks: %g\n", ...
                ntag_, nmats_, nsubs_, ntasks_);

            % write dscalar
            if isemptytext(opts.dscalar_fqfn)
                [~,fp] = myfileparts(mats(end));
                todrop = extractAfter(fp, "-subset-ASHCPPar");
                fp = extractBefore(fp, todrop);
                subPattern = "_sub-" + digitsPattern + "_";
                sesPattern = "_ses-" + digitsPattern + "_";
                fp = replace(fp, subPattern, "_sub-" + ntag_ + "_");
                fp = replace(fp, sesPattern, "_ses-" + ntag_ + "_");

                opts.dscalar_fqfn = fullfile(opts.out_dir, fp + ".dscalar.nii");
            end
            mlraut.Cifti.saveas_dscalar_nii(samples_, opts.dscalar_fqfn);            
        end

        function mats_par_to_subs_dscalar(mats_patt, opts)
            %% intended for, e.g.,
            %  fullfile("/Volumes/PrecunealSSD2/AnalyticSignalHCP", ...
            %      "BrainSmashAdapter_sample_subs_tasks_plvs_as_sub-56_ses-56_proc-RV-std-gsr0-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar-20251022225151_par16.mat")

            arguments
                mats_patt {mustBeTextScalar} = "BrainSmashAdapter_sample_subs_tasks_plvs_as_*_proc-RV-std*_par*.mat"
                opts.mat_varname {mustBeText} = ["measure_data", "samples"]
                opts.out_dir {mustBeTextScalar} = "brainsmash_adapter_subwise"
            end
            assert(~contains(mats_patt, "tasks_physio"))
            if ~isfolder(opts.out_dir)
                mkdir(opts.out_dir)
            end

            % glob
            mats = mglob(mats_patt);

            % sort files by extracting numerics after "_par" in mats
            parPattern = "_par" + digitsPattern;
            parMatches = extract(mats, parPattern);
            parNumbers = str2double(replace(parMatches, "_par", ""));
            [~, sortIdx] = sort(parNumbers);
            mats = mats(sortIdx);

            % segregate
            subid_ = 0;
            for m = mats
                try
                    ld = load(m);
                    assert(isfield(ld, "measure_data"))
                    md = ld.measure_data;
                    assert(contains(m, md.measure_name))

                    for idx_sample = 1:size(md.samples, 1)
                        try
                            dscalar_fqfn = "unknown";
                            samples_ = md.samples(idx_sample, :);
                            samples_ = single(samples_);
                            samples_ = real(samples_);
                            samples_ = samples_';

                            % write dscalar
                            subid_ = subid_ + 1;
                            ntag_ = sprintf("%06d", subid_);
                            [~,fp] = myfileparts(m);
                            todrop = extractAfter(fp, "-subset-ASHCPPar");
                            fp = extractBefore(fp, todrop);
                            subPattern = "_sub-" + digitsPattern + "_";
                            sesPattern = "_ses-" + digitsPattern + "_";
                            fp = replace(fp, subPattern, "_sub-" + ntag_ + "_");
                            fp = replace(fp, sesPattern, "_ses-leq4_");
                            dscalar_fqfn = fullfile(opts.out_dir, fp + ".dscalar.nii");

                            mlraut.Cifti.saveas_dscalar_nii(samples_, dscalar_fqfn);
                            fprintf("%s wrote %s\n", stackstr(), dscalar_fqfn)
                        catch ME
                            handwarning(ME)
                            fprintf("%s: failed writing %s\n", stackstr(), dscalar_fqfn)
                        end
                    end
                catch ME
                    handwarning(ME);
                    fprintf("%s: failed reading %s\n", stackstr(), m)
                end
            end
        end

        function globbing_mat = mats_to_jsons(globbing_mat, opts)
            %% collects opts.measure for all tasks for all sub
            arguments
                globbing_mat {mustBeText} = ...
                    fullfile( ...
                    getenv('SINGULARITY_HOME'), 'AnalyticSignalHCP', 'mlraut_AnalyticSignalHCPPar_globbing.mat')
                opts.globbing_var = "globbed"
                opts.out_dir {mustBeFolder} = "/home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalHCP"
            end

            import mlraut.BrainSmashAdapter.makeConcise

            % `globbing_mat` -> `subs`
            if isscalar(string(globbing_mat)) && endsWith(globbing_mat, ".mat") && ~isemptytext(opts.globbing_var)
                ld = load(globbing_mat);
                subs = ld.(opts.globbing_var);
            else
                subs = globbing_mat;
            end
            subs = convertCharsToStrings(subs);
            subs = asrow(subs);
            fprintf("%s: ", stackstr()); disp(subs);
            
            % sample subjects and tasks per subject (typically 2-4)
            for sub = asrow(subs)
                if isemptytext(sub)
                    continue
                end
                mat_fqfns = mglob( ...
                    fullfile(opts.out_dir, sub, sprintf("sub-%s_ses-*ASHCPPar*.mat", sub)));
                if any(contains(mat_fqfns, "-all"))
                    mat_fqfns = mat_fqfns(~contains(mat_fqfns, "-all"));  % don't double count
                end
                if isemptytext(mat_fqfns)
                    continue
                end

                for mat_fqfn = mat_fqfns
                    tic
                    try
                        ld = load(mat_fqfn);
                        this_subset = ld.this_subset;

                        for f_ = asrow(fieldnames(this_subset))
                            f = convertCharsToStrings(f_);
                            fcontent = this_subset.(f);
                            if isnumeric(fcontent) && numel(fcontent) > 100
                                this_subset.(f) = makeConcise(fcontent);
                            end
                            if strcmp(f, "HCP_signals")                                
                                this_subset.(f) = makeConcise(fcontent);
                            end
                            if strcmp(f, "digital_filter")
                                this_subset.(f) = class(fcontent);
                            end
                            if isa(fcontent, "containers.Map")
                                vals = fcontent.values;
                                this_subset.(f) = struct( ...
                                    'Count', fcontent.Count, ...
                                    'keys', fcontent.keys, ...
                                    'values', makeConcise(vals));
                            end
                        end
                        
                    catch ME
                        fprintf("while working with %s: ", mat_fqfn);
                        handwarning(ME)
                        this_subset = struct('filename', mat_fqfn, 'ME', ME);
                    end

                    % write json
                    json_fqfn = myfileprefix(mat_fqfn) + ".json";                    
                    writelines(jsonencode(this_subset, 'PrettyPrint', true), json_fqfn);

                    fprintf("time working with %s: ", mat_fqfn);
                    toc
                end
            end
        end

        function mat2dscalar(mat, opts)
            arguments
                mat {mustBeNonempty}
                opts.mat_varname {mustBeText} = ["measure_data", "samples"]
                opts.dscalar_fqfn {mustBeTextScalar} = ""
            end
            mat_input = "numeric";
            if istext(mat) && isfile(mat)
                mat_input = mat;
                ld = load(mat);
                assert(isfield(ld, opts.mat_varname(1)))
                mat = ld.(opts.mat_varname(1));
                for vidx = 2:length(opts.mat_varname)
                    mat = mat.(opts.mat_varname(vidx));
                end
            end

            % format mat for cifti
            mat = single(mat);
            mat = real(mat);
            dim_go = find(size(mat) == 91282);  % dim of grayordinate
            if 1 ~= dim_go
                mat = shiftdim(mat, dim_go - 1);
            end
            sz = size(mat);
            if prod(sz(2:end)) > 1
                for dim = 2:length(sz)
                    mat = mean(mat, dim, "omitnan");
                end
            end

            % write dscalar
            if isemptytext(opts.dscalar_fqfn)
                assert(isfile(mat_input))
                opts.dscalar_fqfn = strcat(myfileprefix(mat_input), '.dscalar.nii');
            end
            mlraut.Cifti.saveas_dscalar_nii(mat, opts.dscalar_fqfn);
        end

        function out_fqfn = write_cifti(measure_data, opts)
            arguments
                measure_data struct {mustBeNonempty}
                opts.prefix {mustBeTextScalar} = ""
                opts.out_dir {mustBeFolder} = "/home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalHCP"
                opts.as = []
            end

            measure_name = measure_data.measure_name;
            ntag = measure_data.ntag;  % abbrevs
            ptags = measure_data.ptags;
            out_fqfn =  fullfile(opts.out_dir, ...
                sprintf('%s%s_as_sub-%i_ses-%i_%s', ...
                opts.prefix, measure_name, ntag, ntag, ptags));
            if ~isempty(opts.as)
                opts.as.cifti.write_cifti(measure_data.samples, out_fqfn);
            else
                measure_data.as.cifti.write_cifti(measure_data.samples, out_fqfn);
            end
        end
    end

    %% Helpers

    methods (Static)
        function output = makeConcise(input)
            % This function recursively traverses the input (struct, cell, etc.)
            % and replaces large numeric arrays with their size.

            import mlraut.BrainSmashAdapter.makeConcise

            if isstruct(input)
                % --- It's a struct ---
                % Recurse for each field.
                output = input; % Start by copying
                fields = fieldnames(input);

                for i = 1:numel(input) % Handle struct arrays
                    for f = 1:length(fields)
                        fieldName = fields{f};
                        output(i).(fieldName) = makeConcise(input(i).(fieldName));
                    end
                end

            elseif iscell(input)
                % --- It's a cell array ---
                % Recurse for each element in the cell.
                output = cell(size(input));
                for i = 1:numel(input)
                    output{i} = makeConcise(input{i});
                end

            elseif isnumeric(input) && ~isscalar(input)
                % --- CORE LOGIC: It's a non-scalar numeric array ---
                % Replace it with a concise string representation.
                sz = size(input);

                % Create a string like "1000x1000x3"
                szStr = strjoin(arrayfun(@num2str, sz, 'UniformOutput', false), 'x');

                % Create the final descriptive string
                output = sprintf("[%s %s array]", szStr, class(input));

            else
                % --- It's anything else ---
                % (e.g., scalar, string, logical, char)
                % Keep it as-is.
                output = input;
            end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
