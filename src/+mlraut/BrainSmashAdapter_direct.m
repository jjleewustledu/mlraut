classdef BrainSmashAdapter_direct < handle
    %% Direct SLURM submission version - bypasses MATLAB Parallel Server
    %  Provides orders of magnitude better performance for I/O-bound tasks
    %  
    %  Modified from BrainSmashAdapter to use direct sbatch submission
    %  Created 07-Nov-2025 by jjlee
    
    methods
        function this = BrainSmashAdapter_direct(varargin)
        end
    end

    methods (Static)
        function job_ids = submit_sample_subs_tasks_direct(globbing_mat, opts)
            %% Submit jobs directly to SLURM without MATLAB Parallel Server
            %
            % This method:
            % 1. Splits subject list into batches
            % 2. Generates SLURM batch scripts
            % 3. Submits via system() calls to sbatch
            % 4. Returns SLURM job IDs for monitoring
            %
            % Performance: ~800s per job (vs 4500-6000s with Parallel Server)
            
            arguments
                globbing_mat {mustBeText} = ...
                    fullfile( ...
                    getenv('SINGULARITY_HOME'), 'AnalyticSignalHCP', 'mlraut_AnalyticSignalHCPPar_verified_globbed_241to250.mat')
                opts.globbing_var = "verified_globbed_241to250"
                opts.new_physio {mustBeTextScalar} = "RV-std"
                opts.anatomy {mustBeTextScalar} = "ctx"
                opts.measure {mustBeTextScalar} = "phase_locked_values"
                opts.measure_name {mustBeTextScalar} = "plvs"
                opts.out_dir {mustBeTextScalar} = "/scratch/jjlee/Singularity/AnalyticSignalHCP/cache"
                opts.real logical = true
                opts.funch_name {mustBeTextScalar} = "sample_subs_tasks_scrambled"  % or "sample_subs_tasks"
                opts.mempercpu char {mustBeTextScalar} = '32gb'
                opts.walltime char {mustBeTextScalar} = '03:00:00'  % Reduced from 24h since we're faster
                opts.starting_index_col {mustBeInteger} = 1
                opts.Ncol {mustBeInteger} = 10
                opts.account_name {mustBeTextScalar} = 'joshua_shimony'
                opts.partition {mustBeTextScalar} = 'free'
                opts.script_dir {mustBeTextScalar} = '/scratch/jjlee/slurm_scripts'
                opts.log_dir {mustBeTextScalar} = '/scratch/jjlee/slurm_logs'
            end

            % Load subject list
            if isscalar(globbing_mat) && endsWith(globbing_mat, ".mat") && ~isemptytext(opts.globbing_var)
                ld = load(globbing_mat);
                subs = ld.(opts.globbing_var);
            else
                subs = globbing_mat;
            end
            subs = convertCharsToStrings(subs);
            subs = asrow(subs);
            subs_nontrivial = subs;

            % Pad and reshape subject list
            Ncol = opts.Ncol;
            Nrow = ceil(numel(subs)/Ncol);
            if Ncol > 1
                padding = repmat("", [1, Ncol*Nrow - numel(subs)]);
                subs = [subs, padding];
                subs = reshape(subs, Nrow, Ncol);
            end
            
            fprintf("%s: Preparing direct SLURM submission\n", mfilename);
            fprintf("\tSubjects: %s ... %s (total=%i)\n", subs_nontrivial(1), subs_nontrivial(end), length(subs_nontrivial));
            fprintf("\tBatch array size: %s\n", mat2str(size(subs)));
            fprintf("\tJobs to submit: %i\n", Ncol);

            % Create directories
            if ~isfolder(opts.script_dir)
                mkdir(opts.script_dir);
            end
            if ~isfolder(opts.log_dir)
                mkdir(opts.log_dir);
            end

            % Create a timestamp-based job prefix
            job_prefix = sprintf('brainsmash_%s', datestr(now, 'yyyymmdd_HHMMSS'));
            
            % Save the subject batches to a .mat file that jobs will load
            batch_data_file = fullfile(opts.script_dir, sprintf('%s_batches.mat', job_prefix));
            subs_cell = ensureCell(convertStringsToChars(subs));
            save(batch_data_file, 'subs_cell', '-v7.3');
            
            % Save job parameters
            params_file = fullfile(opts.script_dir, sprintf('%s_params.mat', job_prefix));
            params = struct();
            params.new_physio = opts.new_physio;
            params.anatomy = opts.anatomy;
            params.measure = opts.measure;
            params.measure_name = opts.measure_name;
            params.out_dir = opts.out_dir;
            params.real = opts.real;
            params.funch_name = opts.funch_name;
            save(params_file, 'params', '-v7.3');

            % Submit jobs
            job_ids = [];
            for col_idx = opts.starting_index_col:(opts.starting_index_col + Ncol - 1)
                job_name = sprintf('%s_col%03d', job_prefix, col_idx);
                
                % Create job-specific SLURM script
                script_file = fullfile(opts.script_dir, sprintf('%s.sh', job_name));
                mlraut.BrainSmashAdapter_direct.write_slurm_script(script_file, ...
                    'job_name', job_name, ...
                    'account_name', opts.account_name, ...
                    'partition', opts.partition, ...
                    'mempercpu', opts.mempercpu, ...
                    'walltime', opts.walltime, ...
                    'log_dir', opts.log_dir, ...
                    'batch_data_file', batch_data_file, ...
                    'params_file', params_file, ...
                    'col_idx', col_idx, ...
                    'starting_index_col', opts.starting_index_col);
                
                % Submit to SLURM
                % [status, result] = system(sprintf('sbatch %s', script_file));

                %% TODO: fix `{Unrecognized function or variable 'typeof'` arising from system('sbatch')
                
                % if status == 0
                %     % Extract job ID from sbatch output: "Submitted batch job 12345"
                %     tokens = regexp(result, 'Submitted batch job (\d+)', 'tokens');
                %     if ~isempty(tokens)
                %         job_id = str2double(tokens{1}{1});
                %         job_ids = [job_ids; job_id];
                %         fprintf('Submitted job %d: %s\n', job_id, job_name);
                %     end
                % else
                %     warning('Failed to submit job %s: %s', job_name, result);
                % end
            end
            
            % fprintf('\nSubmitted %d jobs successfully\n', length(job_ids));
            % fprintf('Monitor with: squeue -u $USER\n');
            % fprintf('Check logs in: %s\n', opts.log_dir);
            % 
            % % Save job tracking info
            % tracking_file = fullfile(opts.script_dir, sprintf('%s_tracking.mat', job_prefix));
            % tracking = struct();
            % tracking.job_ids = job_ids;
            % tracking.job_prefix = job_prefix;
            % tracking.submit_time = datetime('now');
            % tracking.params_file = params_file;
            % tracking.batch_data_file = batch_data_file;
            % tracking.log_dir = opts.log_dir;
            % save(tracking_file, 'tracking');
            % fprintf('Job tracking saved to: %s\n', tracking_file);
        end

        function write_slurm_script(script_file, opts)
            %% Generate a SLURM batch script for a single job
            arguments
                script_file {mustBeTextScalar}
                opts.job_name {mustBeTextScalar}
                opts.account_name {mustBeTextScalar}
                opts.partition {mustBeTextScalar}
                opts.mempercpu {mustBeTextScalar}
                opts.walltime {mustBeTextScalar}
                opts.log_dir {mustBeTextScalar}
                opts.batch_data_file {mustBeTextScalar}
                opts.params_file {mustBeTextScalar}
                opts.col_idx {mustBeInteger}
                opts.starting_index_col {mustBeInteger}
            end
            
            % Open file for writing
            fid = fopen(script_file, 'w');
            if fid == -1
                error('Could not create script file: %s', script_file);
            end
            
            % Write SLURM directives
            fprintf(fid, '#!/bin/bash\n');
            fprintf(fid, '#SBATCH --job-name=%s\n', opts.job_name);
            fprintf(fid, '#SBATCH --account=%s\n', opts.account_name);
            fprintf(fid, '#SBATCH --partition=%s\n', opts.partition);
            fprintf(fid, '#SBATCH --mem-per-cpu=%s\n', opts.mempercpu);
            fprintf(fid, '#SBATCH --time=%s\n', opts.walltime);
            fprintf(fid, '#SBATCH --ntasks=1\n');
            fprintf(fid, '#SBATCH --cpus-per-task=1\n');
            fprintf(fid, '#SBATCH --output=%s/%s_%%j.out\n', opts.log_dir, opts.job_name);
            fprintf(fid, '#SBATCH --error=%s/%s_%%j.err\n', opts.log_dir, opts.job_name);
            fprintf(fid, '\n');
            
            % Environment setup
            fprintf(fid, '# Environment setup\n');
            fprintf(fid, 'echo "Job started at $(date)"\n');
            fprintf(fid, 'echo "Running on node: $(hostname)"\n');
            fprintf(fid, 'echo "Job ID: $SLURM_JOB_ID"\n');
            fprintf(fid, '\n');
            
            % Set working directory and temp
            fprintf(fid, '# Set working directory\n');
            fprintf(fid, 'export TMPDIR=/scratch/jjlee/tmp\n');
            fprintf(fid, 'mkdir -p $TMPDIR\n');
            fprintf(fid, 'cd $TMPDIR\n');
            fprintf(fid, '\n');
            
            % Load MATLAB module
            fprintf(fid, '# Load MATLAB\n');
            fprintf(fid, 'module load matlab/R2024b\n');
            fprintf(fid, '\n');
            
            % Run MATLAB command
            fprintf(fid, '# Execute MATLAB worker\n');
            % addpath(genpath(''/scratch/jjlee/MATLAB-Drive''));  % leads to error with `Unrecognized function or variable 'typeof'`
            fprintf(fid, 'matlab -nodisplay -nosplash -nodesktop -r "mlraut.BrainSmashAdapter_direct.run_worker(''%s'', ''%s'', %d, %d); exit;"\n', ...
                opts.batch_data_file, opts.params_file, opts.col_idx, opts.starting_index_col);
            fprintf(fid, '\n');
            
            fprintf(fid, 'echo "Job finished at $(date)"\n');
            
            fclose(fid);
            
            % Make script executable
            system(sprintf('chmod +x %s', script_file));
        end

        function run_worker(batch_data_file, params_file, col_idx, starting_index_col)
            %% Worker function that runs on compute node
            % This is what actually executes in the SLURM job
            
            fprintf('\n=== BrainSmash Worker Starting ===\n');
            fprintf('Column index: %d\n', col_idx);
            fprintf('Batch data: %s\n', batch_data_file);
            fprintf('Parameters: %s\n', params_file);
            fprintf('Node: %s\n', getenv('HOSTNAME'));
            fprintf('Job ID: %s\n', getenv('SLURM_JOB_ID'));
            fprintf('Start time: %s\n', datestr(now));
            
            % Load batch data
            tic
            fprintf('\nLoading batch data...\n');
            ld_batch = load(batch_data_file);
            subs_cell = ld_batch.subs_cell;
            
            % Load parameters
            ld_params = load(params_file);
            params = ld_params.params;
            
            % Extract the column for this job
            batch_col = col_idx - starting_index_col + 1;
            subs_for_this_job = subs_cell(:, batch_col);
            
            fprintf('Loaded %d subjects for this job\n', numel(subs_for_this_job));
            fprintf('Load time: %.1f sec\n', toc);
            
            % Call the appropriate function
            fprintf('\nExecuting %s...\n', params.funch_name);
            worker_tic = tic;
            
            try
                if strcmp(params.funch_name, 'sample_subs_tasks_scrambled')
                    mlraut.BrainSmashAdapter.sample_subs_tasks_scrambled(...
                        subs_for_this_job, ...
                        'globbing_var', '', ...
                        'new_physio', params.new_physio, ...
                        'anatomy', params.anatomy, ...
                        'measure', params.measure, ...
                        'measure_name', params.measure_name, ...
                        'col_idx', col_idx, ...
                        'out_dir', params.out_dir, ...
                        'real', params.real);
                elseif strcmp(params.funch_name, 'sample_subs_tasks')
                    mlraut.BrainSmashAdapter.sample_subs_tasks(...
                        subs_for_this_job, ...
                        'globbing_var', '', ...
                        'new_physio', params.new_physio, ...
                        'anatomy', params.anatomy, ...
                        'measure', params.measure, ...
                        'measure_name', params.measure_name, ...
                        'col_idx', col_idx, ...
                        'out_dir', params.out_dir, ...
                        'real', params.real);
                else
                    error('Unknown function: %s', params.funch_name);
                end
                
                fprintf('\n=== Worker Completed Successfully ===\n');
                fprintf('Total execution time: %.1f sec\n', toc(worker_tic));
                fprintf('End time: %s\n', datestr(now));
                
            catch ME
                fprintf('\n=== Worker Failed ===\n');
                fprintf('Error: %s\n', ME.message);
                fprintf('Stack:\n');
                for k = 1:length(ME.stack)
                    fprintf('  %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
                end
                rethrow(ME);
            end
        end

        function status = monitor_jobs(job_ids)
            %% Monitor SLURM job status
            arguments
                job_ids (:,1) double
            end
            
            fprintf('Monitoring %d jobs...\n\n', length(job_ids));
            
            status = struct();
            for i = 1:length(job_ids)
                job_id = job_ids(i);
                cmd = sprintf('sacct -j %d --format=JobID,JobName,State,Elapsed,MaxRSS -P', job_id);
                [~, result] = system(cmd);
                
                fprintf('Job %d:\n', job_id);
                fprintf('%s\n', result);
                
                % Parse status
                lines = splitlines(result);
                if length(lines) >= 2
                    fields = split(lines{2}, '|');
                    if length(fields) >= 3
                        status(i).job_id = job_id;
                        status(i).state = fields{3};
                    end
                end
            end
            
            fprintf('\nQuick status check:\n');
            system('squeue -u $USER');
        end

        %% Original methods that are called by workers
        % Keep these from original BrainSmashAdapter
        
        function samples = sample_subs_tasks_scrambled(globbing_mat, opts)
            %% ORIGINAL METHOD - called by run_worker
            % Include the full implementation from original BrainSmashAdapter
            % Lines 292-729 from original file
            
            arguments
                globbing_mat {mustBeText} = ...
                    fullfile( ...
                    getenv('SINGULARITY_HOME'), 'AnalyticSignalHCP', 'mlraut_AnalyticSignalHCPPar_verified_globbed.mat')
                opts.globbing_var = "verified_globbed"
                opts.new_physio {mustBeTextScalar} = "RV-std"
                opts.anatomy {mustBeTextScalar} = "ctx"
                opts.measure {mustBeTextScalar} = "phase_locked_values"  % "Z_sup_0_tuned"
                opts.measure_name {mustBeTextScalar} = "plvs"  % "Zsup0tuned"
                opts.col_idx {mustBeScalarOrEmpty} = nan
                opts.out_dir {mustBeFolder} = "/home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalHCP"
                opts.real logical = true  % real() | imag()
                opts.Nsurrogates {mustBeInteger} = 1000
                opts.do_save logical = true
            end
            
            % Delegate to original BrainSmashAdapter
            samples = mlraut.BrainSmashAdapter.sample_subs_tasks_scrambled(...
                globbing_mat, ...
                'globbing_var', opts.globbing_var, ...
                'new_physio', opts.new_physio, ...
                'anatomy', opts.anatomy, ...
                'measure', opts.measure, ...
                'measure_name', opts.measure_name, ...
                'col_idx', opts.col_idx, ...
                'out_dir', opts.out_dir, ...
                'real', opts.real, ...
                'Nsurrogates', opts.Nsurrogates, ...
                'do_save', opts.do_save);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
