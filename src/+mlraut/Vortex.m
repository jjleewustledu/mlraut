classdef Vortex < handle & mlsystem.IHandle
    %  for verifying twistors in the context of spiral waves.
    %  
    %  Created 18-Sep-2025 21:01:39 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 25.1.0.2973910 (R2025a) Update 1 for MACA64.  Copyright 2025 John J. Lee.
    

    properties
        measure
        data_final

        flagRest
        % 1 ~ resting data, 0 ~ task data

        flagTask
        % 1 = language task, original 100 subjects
        % 2 = language task, additional 100 subjects
        % 3 = working memory task
        % for demonstration purpose, keep this as 1 (language task)

        flagSmooth
        % 0 = unsmoothed, raw data
        % 1 = temporally smoothed data (bandpass filtered)
        % 2 = spatiotemporally smoothed (bandpass filtered) data
        % specifies use of Spirals.spiral_detection_surfilt()

        subs
    end

    properties (Dependent)
        home
        Nsub  % number of subjects used for analysis, randomly selected from HCP database (S1200)
        Nt  % abbrev. for num_frames
        spirals
    end

    methods  % GET/SET
        function g = get.home(~)
            if isInParallelWorker
                g = "/scratch/jjlee/Singularity/BrainVortexToolbox";
            elseif contains(hostname, "twistor.attlocal")
                g = "/Users/jjlee/Singularity/BrainVortexToolbox";
            else
                g = "/home/usr/jjlee/Singularity/BrainVortexToolbox";
            end
            g = char(g);
        end

        function g = get.Nsub(this)
            g = length(this.subs);
        end

        function g = get.Nt(this)
            if strcmp(this.measure, "Atlas")
                g = 1199;
                return
            end
            if ~isempty(this.ashcp_)
                g = this.ashcp_.num_bins_angles;
                return
            end
            error("mlraut:RuntimeError", stackstr());
        end

        function g = get.spirals(this)
            g = this.spirals_;
        end
    end

    methods
        function this = Vortex(varargin)
            if ~isempty(varargin)
                this.ashcp_ = mlraut.AnalyticSignalHCPPar.load(varargin{:});
            end
            this.measure = "Atlas";  % "rezeta"
            this.spirals_ = mlraut.Spirals(this);
            this.data_ = mlraut.VortexData(this);
            this.surrogates_ = mlraut.VortexSurrogates(this);
            cd(this.home);

            this.flagRest = 1;
            this.flagTask = 1;
            this.flagSmooth = 2;

            this.subs = ["100206", "100610", ...
                "102614", "107220", "107725", "108020", "111009", "112516", "121315", "121416", "126628", "128026"];
        end

        function dts = dtseries(this, opts)
            arguments
                this mlraut.Vortex
                opts.subid {mustBeScalarOrEmpty} = ""
                opts.measure {mustBeText} = this.measure
                opts.task {mustBeText} = "rfMRI_REST1_LR"
                opts.physio {mustBeText} = "iFV"
                opts.as_fqfn logical = true
            end
            if isnumeric(opts.subid)
                opts.subid = string(opts.subid);
            end

            if strcmp(opts.measure, "Atlas")
                dts = opts.task + "_Atlas.dtseries.nii";
            else
                dts = sprintf( ...
                    "%s_as_sub-%s_ses-%s_proc-%s-gsr0-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar-Vortex.dtseries.nii", ...
                    opts.measure, opts.subid, opts.task, opts.physio);
            end
            if opts.as_fqfn
                dts = fullfile( ...
                    this.home, "Sample Data", "Resting", "Raw Data", opts.subid, dts);
            end
        end

        function s = surf(this, opts)
            arguments
                this mlraut.Vortex
                opts.subid {mustBeScalarOrEmpty} = ""
                opts.hemi {mustBeScalarOrEmpty} = 1
                opts.as_fqfn logical = true
            end
            if isnumeric(opts.subid)
                opts.subid = string(opts.subid);
            end
            if isnumeric(opts.hemi)
                if 1 == opts.hemi
                    hemi = "L";
                else
                    hemi = "R";
                end
            elseif istext(opts.hemi)
                hemi = upper(opts.hemi);                
            end

            s = sprintf("%s.%s.flat.32k_fs_LR.surf.gii", opts.subid, hemi);
            if opts.as_fqfn
                s = fullfile( ...
                    this.home, "Sample Data", "Resting", "Data Pos", opts.subid, s);
            end
        end

        function write_cifti(this, opts)
            arguments
                this mlraut.Vortex
                opts.kind {mustBeTextScalar} = "zeta"
                opts.path {mustBeFolder} = pwd
                opts.tags {mustBeTextScalar} = this.tags()
            end

            assert(~isemnpty(this.ashcp_))
            switch char(opts.kind)
                case "zeta"
                    z = this.aschp_.zeta(this.ashcp_.bold_signal, this.ashcp_.physio_signal);
                    fn = sprintf('rezeta_as_sub-%s_ses-%s_%s', ...
                        this.ashchp_.current_subject, this.ashcp_.current_task, opts.tags);
                    fqfn = fullfile(opts.path, fn);
                    this.ashcp_.write_cifti(real(z), fqfn);
                otherwise
                    error("mlraut:ValueError", "%s: %s unknown", stackstr(), opts.kind)
            end
        end

        function data_out = preprocessing_main(this, opts)
            % hemisphere: 1 for left hemisphere, 2 for right hemisphere, keep this as 1 for demonstration purpose

            arguments
                this mlraut.Vortex
                opts.sub_index {mustBeInteger} = 1
                opts.hemisphere {mustBeInteger} = 1
                opts.flagSmooth {mustBeInteger} = this.flagSmooth
                opts.flagSur logical = false
            end

            if ~opts.flagSur
                data_out = this.data_.preprocessing_main( ...
                    opts.sub_index, ...
                    opts.hemisphere, ...
                    0, ...
                    this.flagRest, ...
                    this.flagTask, ...
                    opts.flagSmooth);
            else
                data_out = this.surrogates_.preprocessing_main( ...
                    opts.sub_index, ...
                    opts.hemisphere, ...
                    1, ...
                    this.flagRest, ...
                    this.flagTask, ...
                    opts.flagSmooth);
            end
        end

        function data_final = results_main_preproc(this)
            %% Start here.

            % struct of cell arrays will gather final data
            data_final.DataOut = cell(1, this.Nsub);
            data_final.real_avg = cell(1, this.Nsub);
            data_final.sur_avg = cell(1, this.Nsub);

            %% Preprocessing of raw fMRI data from HCP
            
            for flagSur = 0:1  %  0 for real data, 1 to generate surrogate data, both are needed for further analysis
                for hemisphere = 2:2  % 1 for left hemisphere, 2 for right hemisphere

                    % ignoring flagTask

                    for flagSmooth_ = 0:2
                        for sidx = 1:this.Nsub
                            if 1 == flagSmooth_
                                data_final.DataOut{sidx} = ...
                                    this.preprocessing_main( ...
                                    sub_index=sidx, hemisphere=hemisphere, flagSur=flagSur, flagSmooth=flagSmooth_);
                            else
                                this.preprocessing_main( ...
                                    sub_index=sidx, hemisphere=hemisphere, flagSur=flagSur, flagSmooth=flagSmooth_);
                            end
                        end
                    end
                end
            end

            this.data_final = data_final;
            save( ...
                fullfile(this.home, "Sample Data", ...
                sprintf("data_final_%s.mat", string(datetime("now", Format="uuuuMMddHHmmss")))), ...
                "data_final", "-v7.3");
        end

        function data_final = results_main(this)

            %% Full-size Vortex detection followed by statistical testing against sprials detected in the null model (surrogate data)
            %  only sprials that contain similarity index higher than 95th percentile of
            %  the null model (surrogate data) are seleceted for further analysis

            data_final = this.data_final;
            real_avg = cell(1, this.Nsub);
            sur_avg = cell(1, this.Nsub);

            parfor sidx = 1:this.Nsub
                for hemisphere = 1:2
                    if hemisphere == 1
                        % load spatiotemporally bandpass filtered (smoothed) real data
                        folder_name = [this.home,'/Sample Data/Resting/Preprocessed Data'];
                        cd(folder_name)
                        filename = ['Preprocessed_spatiotemporalbandpass_data_resting_LEFT_sub',num2str(sidx),'.mat'];
                        ld = load(filename);
                        %                        DataIn_smooth = permute(data_final.DataOut(:,1:175,:,1:this.Nt),[2,3,4,1]);;;
                        DataIn_smooth = ld.DataOut(1:175,:,1:this.Nt);
                        % load temporally bandpass filtered (unsmoothed) real data
                        folder_name = [this.home,'/Sample Data/Resting/Preprocessed Data'];
                        cd(folder_name)
                        filename = ['Preprocessed_temporalbandpass_data_resting_LEFT_sub',num2str(sidx),'.mat'];
                        ld = load(filename);
                        DataIn_unsmooth = ld.DataOut(1:175,:,1:this.Nt);
                        % load spatiotemporally bandpass filtered (smoothed) surrogate data
                        folder_name = [this.home,'/Sample Data/Resting/Preprocessed Data'];
                        cd(folder_name)
                        filename = ['Preprocessed_spatiotemporalbandpass_data_sur_resting_LEFT_sub',num2str(sidx),'.mat'];
                        ld = load(filename);
                        DataIn_sur_smooth = ld.DataOut_smooth;
                        % load temporally bandpass filtered (unsmoothed) surrogate data
                        folder_name = [this.home,'/Sample Data/Resting/Preprocessed Data'];
                        cd(folder_name)
                        filename = ['Preprocessed_temporalbandpass_data_sur_resting_LEFT_sub',num2str(sidx),'.mat'];
                        ld = load(filename);
                        DataIn_sur_unsmooth = ld.DataOut_unsmooth;
                    elseif hemisphere == 2
                        % load spatiotemporally bandpass filtered (smoothed) real data
                        folder_name = [this.home,'/Sample Data/Resting/Preprocessed Data'];
                        cd(folder_name)
                        filename = ['Preprocessed_spatiotemporalbandpass_data_resting_RIGHT_sub',num2str(sidx),'.mat'];
                        ld = load(filename);
                        %                        DataIn_smooth = permute(ld.DataOut(:,1:175,:,1:this.Nt),[2,3,4,1]);;;
                        DataIn_smooth = ld.DataOut(1:175,:,1:this.Nt);
                        % load temporally bandpass filtered (unsmoothed) real data
                        folder_name = [this.home,'/Sample Data/Resting/Preprocessed Data'];
                        cd(folder_name)
                        filename = ['Preprocessed_temporalbandpass_data_resting_RIGHT_sub',num2str(sidx),'.mat'];
                        ld = load(filename);
                        DataIn_unsmooth = ld.DataOut(1:175,:,1:this.Nt);
                        % load spatiotemporally bandpass filtered (smoothed) surrogate data
                        folder_name = [this.home,'/Sample Data/Resting/Preprocessed Data'];
                        cd(folder_name)
                        filename = ['Preprocessed_spatiotemporalbandpass_data_sur_resting_RIGHT_sub',num2str(sidx),'.mat'];
                        ld = load(filename);
                        DataIn_sur_smooth = ld.DataOut_smooth;
                        % load temporally bandpass filtered (unsmoothed) surrogate data
                        folder_name = [this.home,'/Sample Data/Resting/Preprocessed Data'];
                        cd(folder_name)
                        filename = ['Preprocessed_temporalbandpass_data_sur_resting_RIGHT_sub',num2str(sidx),'.mat'];
                        ld = load(filename);
                        DataIn_sur_unsmooth = ld.DataOut_unsmooth;
                    end
                    [real_avg{sidx},sur_avg{sidx}] = ...
                        this.spirals.spiral_detection_surfilt( ...
                            sidx, ...
                            DataIn_smooth, ...
                            DataIn_unsmooth, ...
                            DataIn_sur_smooth, ...
                            DataIn_sur_unsmooth, ...
                            this.home, ...
                            this.flagRest, ...
                            this.flagSmooth, ...
                            this.flagTask, ...
                            hemisphere);
                end
            end

            data_final.real_avg = real_avg;
            data_final.sur_avg = sur_avg;
            this.data_final = data_final;
            save( ...
                fullfile(this.home, "Sample Data", ...
                sprintf("data_final_%s.mat", string(datetime("now", Format="uuuuMMddHHmmss")))), ...
                "data_final", "-v7.3");
        end

        function data_final = results_main_figs(this, opts)
            %% Selectively build figures.
            %  flagSur: 0 for real data, 1 to generate surrogate data

            arguments
                this mlraut.Vortex
                opts.hemisphere {mustBeInteger} = 1
                opts.flagSmooth {mustBeInteger} = this.flagSmooth
                opts.flagSur logical = false
            end

            data_final = this.data_final;

            %% Spiral distribution zscore map: Fig .2d & radii & lifetime & propagation speed

            [data_final.spiral_template_timeavg_accu_stdnorm, ...
                data_final.spiral_radius_accu_avg, ...
                data_final.spiral_duration_accu_avg, ...
                data_final.spiral_count_timeavg, ...
                data_final.spiral_transverse_speed_accu_avg] = ...
                this.spirals.spiral_distribution_zscore_map_speed_duration_radius_count( ...
                this.flagRest,opts.flagSur,opts.flagSmooth,this.flagTask,opts.hemisphere,this.home,this.Nsub);

            %% fMRI amplitude vs. distance from singularity: Fig 2e

            % data_final.distance_from_centre_amplitude_data_pos_accu = ...
            %     distance_vs_amplitude( ...
            %     this.Nsub,this.flagRest,this.flagTask,opts.flagSur,opts.flagSmooth,opts.hemisphere,this.home);

            %% sprial interaction statistics: Fig 4a

            % [data_final.count_repulsion_percent_avg, ...
            %     data_final.count_partialannihilation_percent_avg, ...
            %     data_final.count_fullannihilation_percent_avg] = ...
            %     spiral_interaction_statistics( ...
            %     this.flagRest,opts.flagSur,opts.hemisphere,this.home,this.Nsub,this.flagTask);

            %% extract task-label (language task and working memory task)

            % data_final.task_type_name = TaskLabel_Extract( ...
            %     this.flagRest,this.home,this.Nsub,this.flagTask);

            %% Task-specific trial-averaged spiral distribution map: Fig 5 (language task)

            % [data_final.spiral_distribution_math_listen_mid_avg, ...
            %     data_final.spiral_distribution_math_answer_end_avg, ...
            %     data_final.spiral_distribution_story_listen_mid_avg, ...
            %     data_final.spiral_distribution_story_answer_end_avg] = ...
            %     task_specific_spiral_distribution( ...
            %     this.flagRest,opts.flagSur,opts.hemisphere,this.home,this.Nsub,this.flagTask);

            %% Task-specific trial-averaged spiral distribution, contrast map: Fig 5 (language task)
            % and contrast significance distributions in 7 functional networks: Fig 6a
             
            % [data_final.p_value_negative_log_math_story_listen, ...
            %     data_final.p_value_negative_log_math_story_answer] = ...
            %     spiral_contrast_significance( ...
            %     this.flagRest,opts.flagSur,opts.hemisphere,this.home,this.Nsub,this.flagTask);

            %% sprial centre-based classifier: Fig 6b and c (language task)
             
            % [data_final.classification_accuracy_avg, ...
            %     data_final.classification_accuracy_stderr, ...
            %     data_final.classification_accuracy] = ...
            %     spiral_classifer_language_task( ...
            %     this.flagRest,opts.flagSur,opts.hemisphere,this.home,this.Nsub,this.flagTask);

            %% task-evoked trial-averaged unfiltered fMRI signal: Fig 7a-c (language task)

            % listen_or_answer = 1; % 1 for listening tasks; 2 for answering tasks
            % [data_final.temp1_phase, ...
            %     data_final.Vx_flowmap_norm_phase, ...
            %     data_final.Vy_flowmap_norm_phase] = ...
            %     task_evoked_unfiltered_fMRI_signal( ...
            %     opts.flagSur,opts.hemisphere,this.home,this.Nsub,this.flagTask,listen_or_answer);

            %% PCA analysis: Fig 7d (language task)

            % listen_or_answer = 1; % 1 for listening tasks; 2 for answering tasks
            % motor_or_PCC = 1; % 1 for motor (M1-PMd); 2 for PCC
            % [data_final.score,data_final.latent] = PCA( ...
            %     opts.flagSur,opts.hemisphere,this.home,this.Nsub,this.flagTask,listen_or_answer,motor_or_PCC);

            %% Region of coordination (ROC): Fig 8b & 8d (language task)

            % data_final.region_of_coordination = region_of_coordination_ROC( ...
            %     opts.flagSur,opts.hemisphere,this.home,this.Nsub,this.flagTask);

            %% local phase vector field based classifier: Fig 8c (language task)

            % [data_final.classification_accuracy_avg, ...
            %     data_final.classification_accuracy_stderr] = ...
            %     PhaseVectorAngle_local_classifier( ...
            %     opts.flagSur,opts.hemisphere,this.home,this.Nsub,this.flagTask);

            this.data_final = data_final;
            save( ...
                fullfile(this.home, "Sample Data", ...
                sprintf("data_final_%s.mat", string(datetime("now", Format="uuuuMMddHHmmss")))), ...
                "data_final", "-v7.3");
        end
    end

    methods (Static)
        function this = construct_sample_data()
            % target = "/home/usr/jjlee/mnt/CHPC_scratch/Singularity/BrainVortexToolbox/Sample Data/Resting/Raw Data";
            % assert(isfolder(target))
            % load('mats.mat')
            % for idx = 1:length(mats)
            %     v = mlraut.Vortex.load(mats(idx));
            %     sub = fileparts(mats(idx));
            %     ts = fullfile(target, sub);
            %     ensuredir(ts);
            %     v.write_cifti(path=ts);
            % end
        end
    end

    %% PRIVATE

    properties (Access = private)
        ashcp_
        data_
        Nsub_
        spirals_
        surrogates_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
