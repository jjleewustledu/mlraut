classdef SpiralsFig5 < handle & mlraut.Spirals
    %% line1
    %  line2
    %  
    %  Created 25-Sep-2025 22:57:36 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 24.2.0.2923080 (R2024b) Update 6 for MACA64.  Copyright 2025 John J. Lee.
    
    methods
        function this = SpiralsFig5(varargin)
            this = this@mlraut.Spirals(varargin{:});
        end

        % ,spiral_distribution_story_listen_mid_avg,spiral_distribution_story_answer_end_avg
        function [spiral_distribution_math_listen_mid_avg,spiral_distribution_math_answer_end_avg] = ...
                task_specific_spiral_distribution( ...
                this, flagRest,flagSur,hemisphere,main_folder,~,~)
            %% load task label (language task)

            main_folder = char(main_folder);

            % load task label of each subject
            foldername = [main_folder,'/Sample Data/Resting/Task Label'];
            cd(foldername)
            name = dir(pwd) ;
            file_name2 = ['LanguageOrigTaskLabelAllSubject.mat'];
            load (file_name2);
            for subject = 1:this.Nsub
                fullTime_allsubject{subject} = TaskLabel_AllSubject_language_orig{subject};
            end

            %% Task specific spiral centre distribution of each subject

            for flagSur = 0:1 % both real and surrogate data are needed for further analysis

                % define parameters
                session_duration = 20;
                disp(['extracting task specific spiral centre distributions...'])
                spiral_distribution_math_listen = [];
                spiral_distribution_math_answer = [];
                % spiral_distribution_story_listen = [];
                % spiral_distribution_story_answer = [];
                spiral_distribution_math_question = [];
                % spiral_distribution_story_question = [];
                for subject = 1:this.Nsub
                    foldername = [main_folder,'/Sample Data/Resting/Spiral Detected'];
                    cd(foldername)
                    filename = ['Spiral_detected_surfilt_resting_LEFT_sub',num2str(subject),'.mat'];
                    load(filename)
                    if flagSur == 0
                        spiral_filt_nega = spiral_filt_nega_real_95perc_extend;
                        spiral_filt_pos = spiral_filt_pos_real_95perc_extend;
                    elseif flagSur == 1
                        spiral_filt_nega = spiral_filt_nega_sur;
                        spiral_filt_pos = spiral_filt_pos_sur;
                    end

                    if flagRest == 0
                        temp1_accu_1subject_spiral = zeros(175,251,315);
                    elseif flagRest == 1
                        temp1_accu_1subject_spiral = zeros(175,251,this.Nt);
                    end
                    % clockwise spirals
                    for ipatt = 1:size(spiral_filt_nega,1)
                        for t = 1:size(spiral_filt_nega,2)
                            temp1 = full(spiral_filt_nega{ipatt,t});
                            if nansum(temp1(:))~=0
                                temp1_accu_1subject_spiral(:,:,t) = temp1_accu_1subject_spiral(:,:,t) + temp1;
                            end
                        end
                    end
                    % anticlockwise spirals
                    for ipatt = 1:size(spiral_filt_pos,1)
                        for t = 1:size(spiral_filt_pos,2)
                            temp1 = full(spiral_filt_pos{ipatt,t});
                            if nansum(temp1(:))~=0
                                temp1_accu_1subject_spiral(:,:,t) = temp1_accu_1subject_spiral(:,:,t) + temp1;
                            end
                        end
                    end

                    t_duration = size(temp1_accu_1subject_spiral,3)-session_duration+1; % last time point to be used for analysis
                    % extract task label for each subject
                    temp1 = fullTime_allsubject{subject};
                    temp1_time = temp1(:,3);
                    % find time points for math tasks (listening or answering)
                    count = find(temp1_time==4); % 4 for math listening tasks
                    start_end_time_math_listen = temp1(count,:);
                    count = find(temp1_time==5); % 4 for math listening tasks
                    start_end_time_math_question = temp1(count,:);
                    count = find(temp1_time==6); % 6 for math answering tasks
                    start_end_time_math_answer = temp1(count,:);

                    % time points for story tasks (listening or answering)
                    % count = find(temp1_time==1); % label-1 for story listening tasks
                    % start_end_time_story_listen = temp1(count,:);
                    % count = find(temp1_time==2); % label-1 for story listening tasks
                    % start_end_time_story_question = temp1(count,:);
                    % count = find(temp1_time==3); % label-6 for story answering tasks
                    % start_end_time_story_answer = temp1(count,:);

                    % extract task specifc spiral distribution data: math tasks
                    for trial = 1:size(start_end_time_math_listen,1);
                        temp2_start = start_end_time_math_listen(trial,1);
                        temp2_end = start_end_time_math_listen(trial,2);;
                        temp2_mid = round((temp2_end-temp2_start)./2)+temp2_start;
                        temp2_25perc = round((temp2_end-temp2_start).*0.25)+temp2_start;
                        temp2_75perc = round((temp2_end-temp2_start).*0.75)+temp2_start;
                        t_select = [temp2_start temp2_25perc temp2_mid temp2_75perc temp2_end];
                        temp2_window = [temp2_start:temp2_start+session_duration-1];
                        if temp2_start<1
                            continue
                        end
                        if temp2_end > t_duration
                            break
                        end
                        spiral_distribution_math_listen{trial,subject} = temp1_accu_1subject_spiral(:,:,t_select);
                        if temp2_window(end) < t_duration
                            spiral_distribution_math_listen_1to20{trial,subject} = temp1_accu_1subject_spiral(:,:,temp2_window);
                        end
                    end
                    for trial = 1:size(start_end_time_math_question,1);
                        temp2_start = start_end_time_math_question(trial,1);
                        temp2_end = start_end_time_math_question(trial,2);;
                        temp2_mid = round((temp2_end-temp2_start)./2)+temp2_start;
                        temp2_25perc = round((temp2_end-temp2_start).*0.25)+temp2_start;
                        temp2_75perc = round((temp2_end-temp2_start).*0.75)+temp2_start;
                        t_select = [temp2_start temp2_25perc temp2_mid temp2_75perc temp2_end];
                        temp2_window = [temp2_start:temp2_start+session_duration-1];
                        if temp2_start<1
                            continue
                        end
                        if temp2_end > t_duration
                            break
                        end
                        spiral_distribution_math_question{trial,subject} = temp1_accu_1subject_spiral(:,:,t_select);
                        if temp2_window(end) < t_duration
                            spiral_distribution_math_question_1to20{trial,subject} = temp1_accu_1subject_spiral(:,:,temp2_window);
                        end
                    end
                    for trial = 1:size(start_end_time_math_answer,1);
                        temp2_start = start_end_time_math_answer(trial,1);
                        temp2_end = start_end_time_math_answer(trial,2);;
                        temp2_mid = round((temp2_end-temp2_start)./2)+temp2_start;
                        temp2_25perc = round((temp2_end-temp2_start).*0.25)+temp2_start;
                        temp2_75perc = round((temp2_end-temp2_start).*0.75)+temp2_start;
                        t_select = [temp2_start temp2_25perc temp2_mid temp2_75perc temp2_end];
                        temp2_window = [temp2_start:temp2_start+session_duration-1];
                        if temp2_start<1
                            continue
                        end
                        if temp2_end > t_duration
                            break
                        end
                        spiral_distribution_math_answer{trial,subject} = temp1_accu_1subject_spiral(:,:,t_select);
                        if temp2_window(end) < t_duration
                            spiral_distribution_math_answer_1to20{trial,subject} = temp1_accu_1subject_spiral(:,:,temp2_window);
                        end
                    end
                    % extract task specifc spiral distribution data: story tasks
                    % for trial = 1:size(start_end_time_story_listen,1);
                    %     temp2_start = start_end_time_story_listen(trial,1);
                    %     temp2_end = start_end_time_story_listen(trial,2);;
                    %     temp2_mid = round((temp2_end-temp2_start)./2)+temp2_start;
                    %     temp2_25perc = round((temp2_end-temp2_start).*0.25)+temp2_start;
                    %     temp2_75perc = round((temp2_end-temp2_start).*0.75)+temp2_start;
                    %     t_select = [temp2_start temp2_25perc temp2_mid temp2_75perc temp2_end];
                    %     temp2_window = [temp2_start:temp2_start+session_duration-1];
                    %     if temp2_start<1
                    %         continue
                    %     end
                    %     if temp2_end > t_duration
                    %         break
                    %     end
                    %     spiral_distribution_story_listen{trial,subject} = temp1_accu_1subject_spiral(:,:,t_select);
                    %     if temp2_window(end) < t_duration
                    %         spiral_distribution_story_listen_1to20{trial,subject} = temp1_accu_1subject_spiral(:,:,temp2_window);
                    %     end
                    % end
                    % for trial = 1:size(start_end_time_story_question,1);
                    %     temp2_start = start_end_time_story_question(trial,1);
                    %     temp2_end = start_end_time_story_question(trial,2);;
                    %     temp2_mid = round((temp2_end-temp2_start)./2)+temp2_start;
                    %     temp2_25perc = round((temp2_end-temp2_start).*0.25)+temp2_start;
                    %     temp2_75perc = round((temp2_end-temp2_start).*0.75)+temp2_start;
                    %     t_select = [temp2_start temp2_25perc temp2_mid temp2_75perc temp2_end];
                    %     temp2_window = [temp2_start:temp2_start+session_duration-1];
                    %     if temp2_start<1
                    %         continue
                    %     end
                    %     if temp2_end > t_duration
                    %         break
                    %     end
                    %     spiral_distribution_story_question{trial,subject} = temp1_accu_1subject_spiral(:,:,t_select);
                    %     if temp2_window(end) < t_duration
                    %         spiral_distribution_story_question_1to20{trial,subject} = temp1_accu_1subject_spiral(:,:,temp2_window);
                    %     end
                    % end
                    % for trial = 1:size(start_end_time_story_answer,1);
                    %     temp2_start = start_end_time_story_answer(trial,1);
                    %     temp2_end = start_end_time_story_answer(trial,2);;
                    %     temp2_mid = round((temp2_end-temp2_start)./2)+temp2_start;
                    %     temp2_25perc = round((temp2_end-temp2_start).*0.25)+temp2_start;
                    %     temp2_75perc = round((temp2_end-temp2_start).*0.75)+temp2_start;
                    %     t_select = [temp2_start temp2_25perc temp2_mid temp2_75perc temp2_end];
                    %     temp2_window = [temp2_start:temp2_start+session_duration-1];
                    %     if temp2_start<1
                    %         continue
                    %     end
                    %     if temp2_end > t_duration
                    %         break
                    %     end
                    %     spiral_distribution_story_answer{trial,subject} = temp1_accu_1subject_spiral(:,:,t_select);
                    %     if temp2_window(end) < t_duration
                    %         spiral_distribution_story_answer_1to20{trial,subject} = temp1_accu_1subject_spiral(:,:,temp2_window);
                    %     end
                    % end
                end

                no_of_trial = 0;
                for subject = 1:size(spiral_distribution_math_listen,2)
                    for trial = 1:size(spiral_distribution_math_listen,1)
                        temp1 = full(spiral_distribution_math_listen{trial,subject});
                        if nansum(temp1(:))~=0
                            no_of_trial = no_of_trial+ 1;
                            spiral_distribution_math_listen_mid(:,:,no_of_trial) = temp1(:,:,3);
                        end
                    end
                    spiral_distribution_math_listen_mid_avg = nanmean(spiral_distribution_math_listen_mid,3);
                end

                no_of_trial = 0;
                for subject = 1:size(spiral_distribution_math_answer,2)
                    for trial = 1:size(spiral_distribution_math_answer,1)
                        temp1 = full(spiral_distribution_math_answer{trial,subject});
                        if nansum(temp1(:))~=0
                            no_of_trial = no_of_trial+ 1;
                            spiral_distribution_math_answer_end(:,:,no_of_trial) = temp1(:,:,5);
                        end
                    end
                    spiral_distribution_math_answer_end_avg = nanmean(spiral_distribution_math_answer_end,3);
                end

                % no_of_trial = 0;
                % for subject = 1:size(spiral_distribution_story_listen,2)
                %     for trial = 1:size(spiral_distribution_story_listen,1)
                %         temp1 = full(spiral_distribution_story_listen{trial,subject});
                %         if nansum(temp1(:))~=0
                %             no_of_trial = no_of_trial+ 1;
                %             spiral_distribution_story_listen_mid(:,:,no_of_trial) = temp1(:,:,3);
                %         end
                %     end
                %     spiral_distribution_story_listen_mid_avg = nanmean(spiral_distribution_story_listen_mid,3);
                % end
                % 
                % no_of_trial = 0;
                % for subject = 1:size(spiral_distribution_story_answer,2)
                %     for trial = 1:size(spiral_distribution_story_answer,1)
                %         temp1 = full(spiral_distribution_story_answer{trial,subject});
                %         if nansum(temp1(:))~=0
                %             no_of_trial = no_of_trial+ 1;
                %             spiral_distribution_story_answer_end(:,:,no_of_trial) = temp1(:,:,5);
                %         end
                %     end
                %     spiral_distribution_story_answer_end_avg = nanmean(spiral_distribution_story_answer_end,3);
                % end
                if flagSur == 0
                    foldername = [main_folder,'/Sample Data/Resting/Analysis/'];
                    filename = ['task_specific_spiral_distribution.mat'];
                    save([foldername,filename],'spiral_distribution_math_listen_mid','spiral_distribution_math_listen_mid_avg'...
                        ,'spiral_distribution_math_answer_end_avg','spiral_distribution_math_answer_end'...
                        ,'spiral_distribution_math_listen','spiral_distribution_math_answer'...
                        ,'spiral_distribution_math_listen_1to20','spiral_distribution_math_answer_1to20');
                        % ,'spiral_distribution_story_listen_mid','spiral_distribution_story_listen_mid_avg'...
                        % ,'spiral_distribution_story_answer_end','spiral_distribution_story_answer_end_avg'...
                        %,'spiral_distribution_story_listen','spiral_distribution_story_answer'
                        %,'spiral_distribution_story_listen_1to20','spiral_distribution_story_answer_1to20'
                elseif flagSur == 1
                    foldername = [main_folder,'/Sample Data/Resting/Analysis/'];
                    filename = ['task_specific_spiral_distribution_sur.mat'];
                    save([foldername,filename],'spiral_distribution_math_listen_mid','spiral_distribution_math_listen_mid_avg'...
                        ,'spiral_distribution_math_answer_end_avg','spiral_distribution_math_answer_end'...                        
                        ,'spiral_distribution_math_listen','spiral_distribution_math_answer'...
                        ,'spiral_distribution_math_listen_1to20','spiral_distribution_math_answer_1to20');
                        %,'spiral_distribution_story_listen_mid','spiral_distribution_story_listen_mid_avg'...
                        %,'spiral_distribution_story_answer_end','spiral_distribution_story_answer_end_avg'...
                        %,'spiral_distribution_story_listen','spiral_distribution_story_answer'
                        %,'spiral_distribution_story_listen_1to20','spiral_distribution_story_answer_1to20'
                end
                cd(main_folder)

                if flagSur == 0
                    figure()
                    subplot(2,2,1)
                    imagesc(spiral_distribution_math_listen_mid_avg)
                    set(gca,'ydir','normal')
                    colormap jet
                    colorbar
                    title(['math listen end'])
                    subplot(2,2,2)
                    imagesc(spiral_distribution_math_answer_end_avg)
                    set(gca,'ydir','normal')
                    colormap jet
                    colorbar
                    title(['math answer mid'])
                    % subplot(2,2,3)
                    % imagesc(spiral_distribution_story_listen_mid_avg)
                    % set(gca,'ydir','normal')
                    % colormap jet
                    % colorbar
                    % title(['story listen end'])
                    % subplot(2,2,4)
                    % imagesc(spiral_distribution_story_answer_end_avg)
                    % set(gca,'ydir','normal')
                    % colormap jet
                    % colorbar
                    % title(['story answer mid'])
                end
            end
        end

        % p_value_negative_log_math_story_listen,
        function [p_value_negative_log_math_story_answer] = ...
                spiral_contrast_significance( ...
                this, flagRest,flagSur,hemisphere,main_folder,~,~)
            %% spiral contrast significance map: real data (math vs story listen)
            
            main_folder = char(main_folder);

            % load real data spiral distribution
            foldername = [main_folder,'/Sample Data/Resting/Analysis'];
            cd(foldername)
            filename = ['task_specific_spiral_distribution.mat'];
            load(filename)

            cd(main_folder)
            if hemisphere == 1
                load('parcellation_template7.mat')
            elseif hemisphere == 2
                load('parcellation_template22_RightBrain_subject1-100.mat')
                parcellation_template = parcellation_template22_RightBrain_100sub(:,:,1);
            end

            spiral_distribution_math_listen_1sub_avg = nan(175,251,this.Nsub);
            % spiral_distribution_story_listen_1sub_avg = nan(175,251,this.Nsub);
            for subject = 1:this.Nsub

                % math listen mid (t=3/5)
                no_of_trial = 0;
                spiral_distribution_math_listen_1sub = [];
                for trial = 1:size(spiral_distribution_math_listen,1)
                    temp1 = spiral_distribution_math_listen{trial,subject};
                    if nansum(temp1(:))~=0
                        no_of_trial = no_of_trial+ 1;
                        spiral_distribution_math_listen_1sub(:,:,no_of_trial) = temp1(:,:,3);
                    end
                end
                if nansum(spiral_distribution_math_listen_1sub(:)) == 0
                    continue
                end
                spiral_distribution_math_listen_1sub_avg(:,:,subject) = nanmean(spiral_distribution_math_listen_1sub,3);

                % spiral_distribution_story_listen_1sub = [];
                % no_of_trial = 0;
                % for trial = 1:size(spiral_distribution_story_listen,1)
                %     temp1 = spiral_distribution_story_listen{trial,subject};
                %     if nansum(temp1(:))~=0
                %         no_of_trial = no_of_trial+ 1;
                %         spiral_distribution_story_listen_1sub(:,:,no_of_trial) = temp1(:,:,3);
                %     end
                % end
                % if nansum(spiral_distribution_story_listen_1sub(:)) == 0
                %     continue
                % end
                % spiral_distribution_story_listen_1sub_avg(:,:,subject) = nanmean(spiral_distribution_story_listen_1sub,3);

            end

            spiral_contrast_story_math_listen_allsub_real = spiral_distribution_math_listen_1sub_avg;

            spiral_contrast_story_math_listen_allsub_real_avg = nanmean(spiral_contrast_story_math_listen_allsub_real,3);
            spiral_contrast_story_math_listen_allsub_real_std = nanstd(spiral_contrast_story_math_listen_allsub_real,0,3);
            spiral_contrast_story_math_listen_allsub_real_avg_stdnorm = spiral_contrast_story_math_listen_allsub_real_avg./spiral_contrast_story_math_listen_allsub_real_std;

            %% spiral contrast significance map: real data (math vs story answer)

            spiral_distribution_math_answer_1sub_avg = nan(175,251,this.Nsub);
            % spiral_distribution_story_answer_1sub_avg = nan(175,251,this.Nsub);
            for subject = 1:this.Nsub
                % math answer mid (t=3/5)
                no_of_trial = 0;
                spiral_distribution_math_answer_1sub = [];
                for trial = 1:size(spiral_distribution_math_answer,1)
                    temp1 = spiral_distribution_math_answer{trial,subject};
                    if nansum(temp1(:))~=0
                        no_of_trial = no_of_trial+ 1;
                        spiral_distribution_math_answer_1sub(:,:,no_of_trial) = temp1(:,:,5);
                    end
                end
                if nansum(spiral_distribution_math_answer_1sub(:)) == 0
                    continue
                end
                spiral_distribution_math_answer_1sub_avg(:,:,subject) = nanmean(spiral_distribution_math_answer_1sub,3);

                % spiral_distribution_story_answer_1sub = [];
                % no_of_trial = 0;
                % for trial = 1:size(spiral_distribution_story_answer,1)
                %     temp1 = spiral_distribution_story_answer{trial,subject};
                %     if nansum(temp1(:))~=0
                %         no_of_trial = no_of_trial+ 1;
                %         spiral_distribution_story_answer_1sub(:,:,no_of_trial) = temp1(:,:,5);
                %     end
                % end
                % if nansum(spiral_distribution_story_answer_1sub(:)) == 0
                %     continue
                % end
                % spiral_distribution_story_answer_1sub_avg(:,:,subject) = nanmean(spiral_distribution_story_answer_1sub,3);
            end

            spiral_contrast_story_math_answer_allsub_real = spiral_distribution_math_answer_1sub_avg;

            spiral_contrast_story_math_answer_allsub_real_avg = nanmean(spiral_contrast_story_math_answer_allsub_real,3);
            spiral_contrast_story_math_answer_allsub_real_std = nanstd(spiral_contrast_story_math_answer_allsub_real,0,3);
            spiral_contrast_story_math_answer_allsub_real_avg_stdnorm = spiral_contrast_story_math_answer_allsub_real_avg./spiral_contrast_story_math_answer_allsub_real_std;

            %% spiral contrast significance map: surrogate data (math vs story listen)

            foldername = [main_folder,'/Sample Data/Resting/Analysis'];
            cd(foldername)
            filename = ['task_specific_spiral_distribution_sur.mat'];
            load(filename)

            spiral_distribution_math_listen_1sub_avg = nan(175,251,this.Nsub);
            % spiral_distribution_story_listen_1sub_avg = nan(175,251,this.Nsub);

            for subject = 1:this.Nsub
                % math listen mid (t=3/5)
                no_of_trial = 0;
                spiral_distribution_math_listen_1sub = [];
                for trial = 1:size(spiral_distribution_math_listen,1)
                    temp1 = spiral_distribution_math_listen{trial,subject};
                    if nansum(temp1(:))~=0
                        no_of_trial = no_of_trial+ 1;
                        spiral_distribution_math_listen_1sub(:,:,no_of_trial) = temp1(:,:,3);

                    end
                end
                if nansum(spiral_distribution_math_listen_1sub(:)) == 0
                    continue
                end
                spiral_distribution_math_listen_1sub_avg(:,:,subject) = nanmean(spiral_distribution_math_listen_1sub,3);

                % no_of_trial = 0;
                % spiral_distribution_story_listen_1sub = [];
                % for trial = 1:size(spiral_distribution_story_listen,1)
                %     temp1 = spiral_distribution_story_listen{trial,subject};
                %     if nansum(temp1(:))~=0
                %         no_of_trial = no_of_trial+ 1;
                %         spiral_distribution_story_listen_1sub(:,:,no_of_trial) = temp1(:,:,3);
                % 
                %     end
                % end
                % if nansum(spiral_distribution_story_listen_1sub(:)) == 0
                %     continue
                % end
                % spiral_distribution_story_listen_1sub_avg(:,:,subject) = nanmean(spiral_distribution_story_listen_1sub,3);

            end
            spiral_contrast_story_math_listen_allsub_sur = spiral_distribution_math_listen_1sub_avg;
            spiral_contrast_story_math_listen_allsub_sur_avg = nanmean(spiral_contrast_story_math_listen_allsub_sur,3);

            %% spiral contrast significance map: surrogate data (math vs story answer)

            spiral_distribution_math_answer_1sub_avg = nan(175,251,this.Nsub);
            % spiral_distribution_story_answer_1sub_avg = nan(175,251,this.Nsub);

            for subject = 1:this.Nsub
                % math answer mid (t=3/5)
                no_of_trial = 0;
                spiral_distribution_math_answer_1sub = [];
                for trial = 1:size(spiral_distribution_math_answer,1)
                    temp1 = spiral_distribution_math_answer{trial,subject};
                    if nansum(temp1(:))~=0
                        no_of_trial = no_of_trial+ 1;
                        spiral_distribution_math_answer_1sub(:,:,no_of_trial) = temp1(:,:,5);

                    end
                end
                if nansum(spiral_distribution_math_answer_1sub(:)) == 0
                    continue
                end
                spiral_distribution_math_answer_1sub_avg(:,:,subject) = nanmean(spiral_distribution_math_answer_1sub,3);

                % no_of_trial = 0;
                % spiral_distribution_story_answer_1sub = [];
                % for trial = 1:size(spiral_distribution_story_answer,1)
                %     temp1 = spiral_distribution_story_answer{trial,subject};
                %     if nansum(temp1(:))~=0
                %         no_of_trial = no_of_trial+ 1;
                %         spiral_distribution_story_answer_1sub(:,:,no_of_trial) = temp1(:,:,5);
                % 
                %     end
                % end
                % if nansum(spiral_distribution_story_answer_1sub(:)) == 0
                %     continue
                % end
                % spiral_distribution_story_answer_1sub_avg(:,:,subject) = nanmean(spiral_distribution_story_answer_1sub,3);

            end
            spiral_contrast_story_math_answer_allsub_sur = spiral_distribution_math_answer_1sub_avg;
            spiral_contrast_story_math_answer_allsub_sur_avg = nanmean(spiral_contrast_story_math_answer_allsub_sur,3);

            %% Calculate p-value between real and surrogate sprial contrast map

            p_value = nan(175,251);
            df = nan(175,251);
            for irow = 1:size(spiral_contrast_story_math_listen_allsub_sur,1)
                for icol = 1:size(spiral_contrast_story_math_listen_allsub_sur,2)
                    temp1_sur = spiral_contrast_story_math_listen_allsub_sur(irow,icol,:);
                    if nansum(temp1_sur(:))~=0
                        temp1_real = spiral_contrast_story_math_listen_allsub_real(irow,icol,:);
                        temp1_sur_filt = temp1_sur(temp1_sur~=0);
                        temp1_real_filt = temp1_real(temp1_real~=0);
                        [h,p,ci,stats] = ttest2(temp1_sur_filt(:),temp1_real_filt(:),'Tail','left');
                        p_value(irow,icol) = p;
                        df(irow,icol) = stats.df;
                    end
                end
            end

            p_value_negative_log_math_story_listen = -1.*log10(p_value);
            p_value_negative_log_math_story_listen(p_value_negative_log_math_story_listen==inf) = nan;

            p_value = nan(175,251);
            df = nan(175,251);
            for irow = 1:size(spiral_contrast_story_math_answer_allsub_sur,1)
                for icol = 1:size(spiral_contrast_story_math_answer_allsub_sur,2)
                    temp1_sur = spiral_contrast_story_math_answer_allsub_sur(irow,icol,:);
                    if nansum(temp1_sur(:))~=0
                        temp1_real = spiral_contrast_story_math_answer_allsub_real(irow,icol,:);
                        temp1_sur_filt = temp1_sur(temp1_sur~=0);
                        temp1_real_filt = temp1_real(temp1_real~=0);
                        [h,p,ci,stats] = ttest2(temp1_sur_filt(:),temp1_real_filt(:),'Tail','left');
                        p_value(irow,icol) = p;
                        df(irow,icol) = stats.df;
                    end
                end
            end

            p_value_negative_log_math_story_answer = -1.*log10(p_value);
            p_value_negative_log_math_story_answer(p_value_negative_log_math_story_answer==inf) = nan;

            foldername = [main_folder,'/Sample Data/Resting/Analysis/'];
            filename = ['spiral_contrast_significance.mat'];
            save([foldername,filename],'p_value_negative_log_math_story_listen','p_value_negative_log_math_story_answer')

            figure()
            subplot(1,2,1)
            pcolor(p_value_negative_log_math_story_listen)
            shading interp
            colorbar
            colormap jet
            caxis([1.3,9])
            hold on
            for parcellation_ID = 1:23
                parcellation_template_1par = parcellation_template7;
                parcellation_template_1par(isnan(parcellation_template_1par)) = 0;
                parcellation_template_1par(parcellation_template_1par~=parcellation_ID) = 0;
                B = bwboundaries(parcellation_template_1par,'noholes');
                for k = 1:length(B)
                    boundary = B{k};
                    plot(boundary(:,2),boundary(:,1),'-','linewidth',1,'color',[1,1,1])

                end
            end
            hold off
            subplot(1,2,2)
            pcolor(p_value_negative_log_math_story_answer)
            shading interp
            colorbar
            colormap jet
            caxis([1.3,9])
            hold on
            for parcellation_ID = 1:23
                parcellation_template_1par = parcellation_template7;
                parcellation_template_1par(isnan(parcellation_template_1par)) = 0;
                parcellation_template_1par(parcellation_template_1par~=parcellation_ID) = 0;
                B = bwboundaries(parcellation_template_1par,'noholes');
                for k = 1:length(B)
                    boundary = B{k};
                    plot(boundary(:,2),boundary(:,1),'-','linewidth',1,'color',[1,1,1])

                end
            end
            hold off

            %% parcellation hierachy summary: 7 parcellations (math vs story listen)

            parcellations = [1,2,3,4,5,6,7];
            temp1_ind_accu = [];
            for p = 1:size(parcellations,2)
                temp1_ind = p_value_negative_log_math_story_listen(parcellation_template7 ==parcellations(p));
                temp1_ind(temp1_ind==inf) = nan;
                temp1_ind_accu{p} = temp1_ind;
                Contrast_significance_listen_ParcelAvg(p) = nansum(temp1_ind(:))./size(temp1_ind,1);
                Contrast_significance_listen_ParcelAvg_95CI(p) = nanstd(temp1_ind(:))./sqrt(size(temp1_ind,1)).*1.96;
                Contrast_significance_listen_ParcelAvg_stderr(p) = nanstd(temp1_ind(:))./sqrt(size(temp1_ind,1));
            end

            errlow = - Contrast_significance_listen_ParcelAvg_stderr;
            errhigh =  Contrast_significance_listen_ParcelAvg_stderr;

            x = [1:7];
            figure()
            bar(x,Contrast_significance_listen_ParcelAvg)
            plotName = [{'VIS'},{'SMN'},{'AUD'},{'CON'},{'DAN'},{'FPN'},{'DMN'}] ;
            set(gca,'xticklabel',plotName)
            hold on
            er = errorbar(x,Contrast_significance_listen_ParcelAvg,errlow,errhigh);
            er.Color = [0 0 0];
            er.LineStyle = 'none';
            hold off
            % ylim([0,0.12])
            title(['spiral contrast significance in 7 networks, Story vs math listen'])

            parcellations = [1,2,3,4,5,6,7];
            temp1_ind_accu = [];
            for p = 1:size(parcellations,2)
                temp1_ind = p_value_negative_log_math_story_answer(parcellation_template7 ==parcellations(p));
                temp1_ind(temp1_ind==inf) = nan;
                temp1_ind_accu{p} = temp1_ind;
                Contrast_significance_answer_ParcelAvg(p) = nansum(temp1_ind(:))./size(temp1_ind,1);
                Contrast_significance_answer_ParcelAvg_95CI(p) = nanstd(temp1_ind(:))./sqrt(size(temp1_ind,1)).*1.96;
                Contrast_significance_answer_ParcelAvg_stderr(p) = nanstd(temp1_ind(:))./sqrt(size(temp1_ind,1));
            end

            errlow = - Contrast_significance_answer_ParcelAvg_stderr;
            errhigh =  Contrast_significance_answer_ParcelAvg_stderr;

            x = [1:7];
            figure()
            bar(x,Contrast_significance_answer_ParcelAvg)
            plotName = [{'VIS'},{'SMN'},{'AUD'},{'CON'},{'DAN'},{'FPN'},{'DMN'}] ;
            set(gca,'xticklabel',plotName)
            hold on
            er = errorbar(x,Contrast_significance_answer_ParcelAvg,errlow,errhigh);
            er.Color = [0 0 0];
            er.LineStyle = 'none';
            hold off
            % ylim([0,0.12])
            title(['spiral contrast significance in 7 networks, Story vs math answer'])

            foldername = [main_folder,'/Sample Data/Resting/Analysis/'];
            filename = ['spiral_contrast_significance_7networks.mat'];
            save([foldername,filename], ...
                'plotName','Contrast_significance_listen_ParcelAvg','Contrast_significance_listen_ParcelAvg_stderr','Contrast_significance_answer_ParcelAvg','Contrast_significance_answer_ParcelAvg_stderr')
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
