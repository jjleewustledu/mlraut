classdef AnalyticSignalGBMPar < handle & mlraut.AnalyticSignalGBM
    %% line1
    %  line2
    %  
    %  Created 29-Jun-2023 11:22:03 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 9.14.0.2286388 (R2023a) Update 3 for MACI64.  Copyright 2023 John J. Lee.
    
    methods (Static)
        function this = median_twistor(physio)
            arguments
                physio {mustBeTextScalar} = 'physio_iFV'
            end

            this = mlraut.AnalyticSignalGBMPar( ...
                subjects={'sub-I3CR0023'}, ...
                do_7T=false, ...
                do_resting=true, ...
                do_task=false, ...
                do_save=false, ...
                do_save_dynamic=false, ...
                do_save_ciftis=false, ...
                hp_thresh=0.01, ...
                lp_thresh=0.1, ...
                global_signal_regression=true, ...
                tags="AnalyticSignalGBMPar-median-twistor");

            mats = asrow(glob(fullfile(this.out_dir, physio, '*', 'sub-*_ses-*AnalyticSignalGBM*.mat')));
            n = length(mats);
            nx = this.num_nodes;

            X_ = zeros(1, nx);
            Y_ = zeros(1, nx);
            Z_ = zeros(1, nx);
            T_ = zeros(1, nx);

            for mat = mats
                % tic
                ld = load(mat{1});
                this_subset = ld.this_subset;
                try
                    ksi = this_subset.bold_signal;
                catch ME
                    if strcmp(ME.identifier, 'MATLAB:nonExistentField')
                        ksi = this_subset.analytic_signal.*this_subset.physio_signal;  % overly normalized in this_subset
                    else
                        rethrow(ME)
                    end
                end
                eta = this_subset.physio_signal;
                X = median(this.build_final_normalization((ksi.*conj(eta) + eta.*conj(ksi))/sqrt(2)), 1);
                Y = median(this.build_final_normalization((ksi.*conj(eta) - eta.*conj(ksi))/sqrt(2i)), 1);
                Z = median(this.build_final_normalization((ksi.*conj(ksi) - eta.*conj(eta))/sqrt(2)), 1);
                T = median(this.build_final_normalization((ksi.*conj(ksi) + eta.*conj(eta))/sqrt(2)), 1);

                X_ = X_ + X/n;
                Y_ = Y_ + Y/n;
                Z_ = Z_ + Z/n;
                T_ = T_ + T/n;
                % toc
            end

            this.write_ciftis( ...
                X_, sprintf('X_as_sub-all_ses-all_%s', this.tags), ...
                do_final_normalization=false, ...
                do_save_dynamic=false);
            this.write_ciftis( ...
                Y_, sprintf('Y_as_sub-all_ses-all_%s', this.tags), ...
                do_final_normalization=false, ...
                do_save_dynamic=false);
            this.write_ciftis( ...
                Z_, sprintf('Z_as_sub-all_ses-all_%s', this.tags), ...
                do_final_normalization=false, ...
                do_save_dynamic=false);
            this.write_ciftis( ...
                T_, sprintf('T_as_sub-all_ses-all_%s', this.tags), ...
                do_final_normalization=false, ...
                do_save_dynamic=false);
            this.write_ciftis( ...
                T_+Z_, sprintf('T+Z_as_sub-all_ses-all_%s', this.tags), ...
                do_final_normalization=false, ...
                do_save_dynamic=false);
        end
        
        function parcall(cores, opts)
            arguments
                cores {mustBeScalarOrEmpty} = 32
                opts.config_hemispheres {mustBeTextScalar} = "lesionR-CE" % "lesionR-iFV" "nolesion" "alllesion" ""
            end

            root_dir = '/home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalGBM/analytic_signal/dockerout/ciftify';
            out_dir = '/home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalGBM/analytic_signal/matlabout/lesionR-CE';
            ensuredir(out_dir);
            tasks = {'ses-1_task-rest_run-01_desc-preproc', 'ses-1_task-rest_run-02_desc-preproc'};

            g = glob(fullfile(root_dir, 'sub-*'));
            %g = flip(g);
            g = cellfun(@(x) basename(x), g, UniformOutput=false);
            %g = g(1:end);
            leng = length(g);
            parfor (idxg = 1:leng, cores)
                try
                    this = mlraut.AnalyticSignalGBMPar(subjects=g(idxg), ...
                        root_dir=root_dir, out_dir=out_dir, tasks=tasks, ...
                        source_physio='ROI', roi_fileprefix='CE_on_T1w');
                    this.config_hemispheres = opts.config_hemispheres; %#ok<PFBNS>
                    call(this);
                catch ME
                    handwarning(ME)
                end
            end
        end

        function parcall_physios(cores)   
            %% left insula, no midline shift

            arguments
                cores {mustBeScalarOrEmpty} = 3
            end

            root_dir = fullfile( ...
                getenv('SINGULARITY_HOME'), ...
                'AnalyticSignalGBM/analytic_signal/dockerout/ciftify');
            cd(root_dir);

            g = convertStringsToChars("sub-" + mlraut.AnalyticSignalGBMPar.SUBS);
            leng = length(g);
            for idxg = 1:1
            % parfor (idxg = 1:leng, cores)
                try
                    tic

                    out_dir = fullfile( ...
                        getenv('SINGULARITY_HOME'), ...
                        'AnalyticSignalGBM/analytic_signal/matlabout/physio_iFV');
                    ensuredir(out_dir);
                    this = mlraut.AnalyticSignalGBM( ...
                        subjects=g(idxg), ...
                        tasks={'ses-1_task-rest_run-all_desc-preproc'}, ...
                        do_save=true, ...
                        do_save_ciftis=true, ...
                        out_dir=out_dir, ...
                        source_physio='iFV');
                    call(this);

                    out_dir = fullfile( ...
                        getenv('SINGULARITY_HOME'), ...
                        'AnalyticSignalGBM/analytic_signal/matlabout/physio_CE');
                    ensuredir(out_dir);
                    ce = mlfourd.ImagingContext2( ...
                        fullfile(root_dir, g(idxg), 'MNINonLinear', 'CE_on_T1w.nii.gz'));
                    this = mlraut.AnalyticSignalGBM( ...
                        subjects=g(idxg), ...
                        tasks={'ses-1_task-rest_run-all_desc-preproc'}, ...
                        do_save=true, ...
                        do_save_ciftis=true, ...
                        out_dir=out_dir, ...
                        roi=ce);
                    call(this);

                    out_dir = fullfile( ...
                        getenv('SINGULARITY_HOME'), ...
                        'AnalyticSignalGBM/analytic_signal/matlabout/physio_WT');
                    ensuredir(out_dir);
                    wt = mlfourd.ImagingContext2( ...
                        fullfile(root_dir, g(idxg), 'MNINonLinear', 'WT_on_T1w.nii.gz'));
                    this = mlraut.AnalyticSignalGBM( ...
                        subjects=g(idxg), ...
                        tasks={'ses-1_task-rest_run-all_desc-preproc'}, ...
                        do_save=true, ...
                        do_save_ciftis=true, ...
                        out_dir=out_dir, ...
                        roi=wt);
                    call(this);

                    toc
                catch ME
                    handwarning(ME)
                end
            end

            % Elapsed time is ___ seconds.
        end
    end

    methods
        function this = AnalyticSignalGBMPar(varargin)
            this = this@mlraut.AnalyticSignalGBM(varargin{:});
        end

        function this = call(this)
            %% CALL all subjects

            % exclude subjects

            out_dir_ = this.out_dir;
            for s = 1:this.num_sub
                this.current_subject = this.subjects{s};
                if ~contains(out_dir_, this.current_subject)
                    proposed_dir = fullfile(out_dir_, this.current_subject);
                    ensuredir(proposed_dir);
                    this.out_dir = proposed_dir;
                end 
                this.call_subject(s);
            end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
