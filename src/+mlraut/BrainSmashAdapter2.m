classdef BrainSmashAdapter2 < handle

    methods
        function this = BrainSmashAdapter2(varargin)
        end
    end

    methods (Static)
        function [j,c] = cluster_mats_par_to_group_npy(mats_patt, opts)
            %% Args:
            %    mats_patt {mustBeText} = ...
            %        fullfile( ...
            %        getenv("SINGULARITY_HOME"), "AnalyticSignalHCP", ...
            %        "BrainSmashAdapter_sample_subs_tasks_MEASURENAME_as_*_proc-PHYSIO*_par*.mat")
            %    opts.new_physio {mustBeText} = ["RV-std", "iFV", "sFV", "3rdV", "latV", "centrumsemiovale", "dmn", "sms", "vis"]
            %    opts.measure_name {mustBeText} = ["plvs", "Zsup0tuned", "Zsup1tuned"]
            %    opts.out_dir {mustBeTextScalar} = "/scratch/jjlee/Singularity/AnalyticSignalHCP"
            %    opts.mempercpu char {mustBeTextScalar} = "8gb";
            %    opts.walltime char {mustBeTextScalar} = "24:00:00";
                
            arguments
                mats_patt {mustBeText} = ...
                    fullfile( ...
                    getenv("SINGULARITY_HOME"), "AnalyticSignalHCP", ...
                    "BrainSmashAdapter", ...
                    "sample_subs_tasks", ...
                    "BrainSmashAdapter_sample_subs_tasks_MEASURENAME_as_*_proc-PHYSIO*_par*.mat")
                opts.new_physio {mustBeText} = "RV-std"  % ["RV-std", "iFV", "sFV", "3rdV", "latV", "centrumsemiovale", "dmn", "sms", "vis"]
                opts.measure_name {mustBeText} = "plvs"  % ["plvs", "Zsup0tuned", "Zsup1tuned", "Zsup2tuned"]
                opts.out_dir {mustBeTextScalar} = ...
                    "/scratch/jjlee/Singularity/AnalyticSignalHCP/BrainSmashAdapter/sample_subs_tasks"
                opts.mempercpu char {mustBeTextScalar} = "16gb";
                opts.walltime char {mustBeTextScalar} = "3:00:00";
            end

            % additional internal params
            Ncol = 16;  % Ncol ~ num jobs; Nrow ~ ceil(1113/16) = 70 subs per job
            account_name = "joshua_shimony";

            for measure_idx = 1:length(opts.measure_name)
                for physio_idx = 1:length(opts.new_physio)

                    % mats_patt -> mats
                    mats_patt = strrep(mats_patt, "MEASURENAME", opts.measure_name{measure_idx});
                    mats_patt = strrep(mats_patt, "PHYSIO", opts.new_physio{physio_idx});

                    %% contact cluster slurm

                    warning('off', 'MATLAB:legacy:batchSyntax');
                    warning('off', 'parallel:convenience:BatchFunctionNestedCellArray');
                    warning('off', 'MATLAB:TooManyInputs');

                    account_name = char(account_name);
                    opts.mempercpu = convertStringsToChars(opts.mempercpu);
                    opts.walltime = convertStringsToChars(opts.walltime);
                    mats_patt = convertStringsToChars(mats_patt);
                    opts.out_dir = convertStringsToChars(opts.out_dir);

                    c = mlraut.CHPC3.propcluster_free(account, mempercpu=opts.mempercpu, walltime=opts.walltime);
                    try
                        j = c.batch( ...
                            @mlraut.BrainSmashAdapter2.mats_par_to_group_npy, ...
                            1, ...
                            {mats_patt, ...
                            'out_dir', opts.out_dir}, ...
                            'CurrentFolder', opts.out_dir, ...
                            'AutoAddClientPath', false);
                    catch ME
                        handwarning(ME)
                    end

                    warning('on', 'MATLAB:legacy:batchSyntax');
                    warning('on', 'parallel:convenience:BatchFunctionNestedCellArray');
                    warning('on', 'MATLAB:TooManyInputs');

                end
            end
        end

        function mats = mats_par_to_group_npy(mats_patt, opts)
            %% intended for, e.g.,
            %  fullfile("/Volumes/PrecunealSSD2/AnalyticSignalHCP", ...
            %      "BrainSmashAdapter_sample_subs_tasks_plvs_as_sub-56_ses-56_proc-RV-std-gsr0-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar-20251022225151_par16.mat")
            %
            %  Args:
            %    mats_patt may be pattern or string array of filenames
            %    opts.npy_fqfn {mustBeTextScalar} = ""
            %    opts.out_dir {mustBeTextScalar} = ""
            %
            %  Returns:
            %    mats

            arguments
                mats_patt {mustBeText} = "BrainSmashAdapter_sample_subs_tasks_plvs_as_*_proc-RV-std*_par*.mat"
                opts.npy_fqfn {mustBeTextScalar} = ""
                opts.out_dir {mustBeTextScalar} = ""
            end
            mats_patt = convertCharsToStrings(mats_patt);
            assert(~contains(mats_patt, "tasks_physio"))
            if ~isemptytext(opts.out_dir) && ~isfolder(opts.out_dir)
                mkdir(opts.out_dir)
            end
            
            % aggregate
            if isscalar(mats_patt) && contains(mats_patt, "*")
                mats = mglob(mats_patt);
            else
                mats = mats_patt;
            end            
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
            samples_ = single(samples_);
            samples_ = real(samples_);  % N_sub x N_go
            fprintf("%s:\n", stackstr());
            fprintf("\tntag: %g, nmats: %g, nsubs: %g, ntasks: %g\n", ...
                ntag_, nmats_, nsubs_, ntasks_);

            % write npy
            if isemptytext(opts.npy_fqfn)
                [~,fp] = myfileparts(mats(end));
                todrop = extractAfter(fp, "-subset-ASHCPPar");
                fp = extractBefore(fp, todrop);
                subPattern = "_sub-" + digitsPattern + "_";
                sesPattern = "_ses-" + digitsPattern + "_";
                fp = replace(fp, subPattern, "_sub-" + ntag_ + "_");
                fp = replace(fp, sesPattern, "_ses-" + ntag_ + "_");

                opts.npy_fqfn = fullfile(opts.out_dir, fp + ".npy");
            end            
            
            % Ensure MATLAB is using your conda Python
            nilpy = '/home/usr/jjlee/.conda/envs/matlab_py310/bin/python';
            chpcpy = fullfile(getenv('HOME'), 'miniconda3/envs/py3/bin/python');
            if isfile(nilpy)
                pe = pyenv('Version', nilpy);
            elseif isfile(chpcpy)
                pe = pyenv('Version', chpcpy);
            else
                error("mlraut:RunTimeError", "%s: no python found", stackstr());
            end
            disp(pe.Version)
            disp(pe.Executable)

            % Save npy
            py.numpy.save(opts.npy_fqfn, samples_);
        end        
    end
end
