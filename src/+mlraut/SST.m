classdef SST
    %% line1
    %  line2
    %  
    %  Created 04-Apr-2022 15:30:06 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 9.12.0.1884302 (R2022a) for MACI64.  Copyright 2022 John J. Lee.
    
    methods (Static)
        function createFolders(varargin)
            %% Adheres to Ryan's conventions for 
            %  "Global waves synchronize the brain"s functional systems with fluctuating arousal"
            %  Args:
            %      subjectID (required text)
            %      toglob (required text)
            %      desc (text): to embed in folder names and ConditionLabel.  Default is 'unknown'.

            ip = inputParser;
            addRequired(ip, 'subjectID', @istext);
            addRequired(ip, 'toglob', @istext);
            addParameter(ip, 'desc', 'unknown', @istext);
            parse(ip, varargin{:});
            ipr = ip.Results;
            if ~startsWith(ipr.subjectID, 'sub-')
                ipr.subjectID = strcat('sub-', ipr.subjectID);
            end
            ipr.desc = strrep(ipr.desc, '_', '-');

            sub_path = fullfile(getenv('SINGULARITY_HOME'), 'arousal_waves', ipr.subjectID, '');
            ensuredir(sub_path);
            sesnum = 0;
            for fn = globT(ipr.toglob)
                % edfreadOri
                [hdr,record] = edfreadOri(fn{1});
                disp(hdr)
                %assert(all(hdr.samples == hdr.samples(1))) % 256
                %assert(all(hdr.frequency == hdr.frequency(1))) % 256

                % ensure dirs
                s = strcat(hdr.startdate, '-', hdr.starttime);
                dt = datetime(s, 'InputFormat', 'dd.MM.yy-HH.mm.ss');
                numch = size(record, 1);
                sesnum = sesnum + 1;
                ses_path = fullfile(sub_path, ...
                    sprintf('ses-%s_desc-%s_proc-mat_ecog', datestr(dt, 'yyyymmddHHMMSS'), ipr.desc));
                ensuredir(ses_path);
                
                %% save session

                pwd0 = pushd(ses_path);

                % save channels mat
                for ch = 1:numch
                    eval(sprintf('ECoGData_ch%i = record(%i,:);', ch, ch))
                    eval(sprintf("save('ECoG_ch%i.mat', 'ECoGData_ch%i');", ch, ch))
                end

                % save hdr mat
                hdr.filename = fullfile(pwd0, fn{1});
                save('hdr.mat', 'hdr')

                % save ECoGTime.mat
                dT = 1/hdr.frequency(1); % sec
                T = size(record,2)/hdr.frequency(1) - dT; % sec
                ECoGTime  = 0:dT:T;
                assert(length(ECoGTime) == size(record,2))
                save('ECoGTime.mat', 'ECoGTime');

                % save Condition.mat
                ConditionIndex = 1;
                ConditionLabel = {ipr.desc};
                ConditionTime = 0;
                save('Condition.mat', 'ConditionIndex', 'ConditionLabel', 'ConditionTime');

                popd(pwd0);
            end
        end
        function createRunsText()
        end
    end

    properties
        avg_alpha
        avg_delta
        avg_gamma
        avg_sst
        f
        Fs_smooth = 5
        null_sst
        sst_count
        subjectPath        
    end

    properties (Dependent)
    end

    methods
        function this = SST(varargin)
            %% SST 
            %  Args:
            %      subjectPath (folder)
            %      runsPatt (text): e.g., '*_desc-clinical_collection_*' in subjectPath.
            
            ip = inputParser;
            addParameter(ip, "subjectPath", pwd, @isfolder)
            addParameter(ip, "runsPatt", '*', @istext)
            parse(ip, varargin{:})
            ipr = ip.Results;
            this.subjectPath = ipr.subjectPath;
            if ~startsWith(ipr.runsPatt, this.subjectPath)
                ipr.runsPatt = fullfile(this.subjectPath, ipr.runsPatt);
            end
            
            %% Setup
            addpath(genpath(fullfile( ...
                getenv('HOME'), 'MATLAB-Drive', 'arousal-waves-main', 'Dependencies', 'chronux_2_12')));
            runs = glob(ipr.runsPatt);
            
            % Initialize SST variables
            avg_sst_ = zeros(201,101);
            avg_gamma_ = zeros(201,128);
            avg_alpha_ = avg_gamma_;
            avg_delta_ = avg_gamma_;
            avg_acfs_ = avg_gamma_;
            sst_count_ = 0;
            null_sst_ = zeros(201,101);
            null_count = 0;
            
            for r = 1:length(runs)
                disp(['Processing ' runs{r}])
                disp([num2str(r) ' out of ' num2str(numel(runs))])

                ld = load(fullfile(runs{r}, 'hdr.mat'));            
                num_nodes = ld.hdr.ns;
                Fs = ld.hdr.frequency(1);
                time_vec = importdata(fullfile(runs{r}, 'ECoGTime.mat'));
            
                start = 10*Fs; % int index, starting ten seconds in
                stop = length(time_vec); % int index, to end of data collection
                num_samples = length(start:stop);
            
                %% Load data
                disp('Loading data...')
                raw_data = single(zeros(num_samples,num_nodes));
                for i = 1:num_nodes
                    tic
                    temp = importdata(fullfile(runs{4}, strcat('ECoG_ch', num2str(i), '.mat')));
                    raw_data(:,i) = temp(start:stop);
                    toc
                end
            
                good = true(num_nodes,1);
                if isfield(ld.hdr, 'badchannels')
                    good(ld.hdr.badchannels) = false; % e.g., [53, 73]
                end
                raw_data(:,~good) = nan;
                raw_data = bsxfun(@minus, raw_data, nanmean(raw_data,2)); % all good channels
                
                %% Time-frequency analysis
                
                [b1,a1] = butter(4,[49,51]/(Fs/2),'stop'); % notch
                [b2,a2] = butter(4,[99,101]/(Fs/2),'stop'); % notch
                filt_data = filtfilt(b1,a1,double(raw_data));
                filt_data = filtfilt(b2,a2,double(filt_data));
                    
                params.Fs = Fs;
                params.fpass = [1 100];
                params.trialave = 0;
                params.tapers = [3 5];
                [S,~,this.f] = mtspecgramc(filt_data,[1,1/this.Fs_smooth],params);
            
                power = real(10*log10(S));
                power_zm = bsxfun(@minus,power,nanmean(power,1));
                global_power = nanmean(power_zm,3)';
                
            %     figure;
            %     surf(t,this.f,global_power,'EdgeColor','none');
            %     axis xy;axis tight;colormap(jet);view(0,90);caxis([-5,5])
            %     colorbar;c=colorbar;c.Label.String = 'Power (dB)';
            %     xlabel('Seconds')
            %     ylabel('Frequency (Hz)')
            
                %% SST detection
            
                % Get low-frequency power
                [~,ind_low] = min(abs(this.f-0));
                [~,ind_high] = min(abs(this.f-4));
                lf_power = squeeze(nanmean(power_zm(:,ind_low:ind_high,:),2));
            
                % Get mid-frequency power
                [~,ind_low] = min(abs(this.f-9));
                [~,ind_high] = min(abs(this.f-21));
                mf_power = squeeze(nanmean(power_zm(:,ind_low:ind_high,:),2));
            
                % Get gamma BLP
                [~,ind_low] = min(abs(this.f-42));
                [~,ind_high] = min(abs(this.f-87));
                hf_power = squeeze(nanmean(power_zm(:,ind_low:ind_high,:),2));
                
                % Find SSTs
                stdevs = nanstd(lf_power);
                means = nanmean(lf_power);
                sst_inds = sum(bsxfun(@gt,lf_power,(means+stdevs)),2)>(.4*sum(good));
               
                [b,a] = butter(2,[.01,.05]/(this.Fs_smooth/2)); % infra-slow envelope of gamma blp
                hf_power = filtfilt(b,a,hf_power);
                mf_power = filtfilt(b,a,mf_power);
                lf_power = filtfilt(b,a,lf_power);
            
            
                last_sst = -100;
                for i = 1:length(sst_inds)
                    if sst_inds(i) && i>101 && i<(length(sst_inds)-100) && (i-3*this.Fs_smooth)>last_sst
                        avg_sst_ = avg_sst_ + global_power(:,i-100:i+100)';
                        avg_gamma_ = avg_gamma_ + hf_power(i-100:i+100,:);
                        avg_alpha_ = avg_alpha_ + mf_power(i-100:i+100,:);
                        avg_delta_ = avg_delta_ + lf_power(i-100:i+100,:);
            
                        sst_count_ = sst_count_+1;
                        last_sst = i;
                    end
                end
            
                % Generate null distribution
                last_sst = -100;
                null_inds = sst_inds(randperm(length(sst_inds)));
                for i = 1:length(null_inds)
                    if null_inds(i) && i>101 && i<(length(null_inds)-100) && (i-3/.2)>last_sst && null_count <= sst_count_
                        null_sst_ = cat(3,null_sst_,global_power(:,i-100:i+100)');
                        null_count = null_count+1;
                        last_sst = i;
                    end
                end
            end
            
            avg_sst_ = avg_sst_./sst_count_;
            this.avg_gamma = avg_gamma_./sst_count_;
            this.avg_alpha = avg_alpha_./sst_count_;
            this.avg_delta = avg_delta_./sst_count_;
            
            % normalize by null
            avg_sst_ = avg_sst_-squeeze(mean(null_sst_,3));
            this.avg_sst = avg_sst_./std(null_sst_,0,3);            
            this.null_sst = null_sst_;
            this.sst_count = sst_count_;
        end
        function plot_avg_gamma(this)
            figure;
            plot(nanmean(this.avg_gamma,2))
        end
        function plot_spectogram(this, map_filename)
            % Plot SST spectrogram and gamma band time series

            assert(isfile(map_filename))
            load(map_filename) % e.g., [root_dir '/' subject 'Map.mat'], containing I, X, Y
            
            avg_sst2 = this.avg_sst.*std(this.null_sst,0,3) + squeeze(mean(this.null_sst,3));
            sst_count2 = this.sst_count;
            avg_sst2 = avg_sst2.*sst_count2;
            avg_sst_ = (avg_sst1+avg_sst2)./(sst_count1+sst_count2);
            
            figure;
            surf((1/this.Fs_smooth)*(-100:100),this.f,avg_sst_','EdgeColor','none');
            axis xy;axis tight;colormap(viridis);view(0,90);
            colorbar;caxis([-.25,.25])
            set(gca,'fontsize',15,'fontweight','bold')            
            
            % Plot lags
            %avg_gamma = SST1.avg_gamma;
            %load('XX/Chibi/ChibiMap.mat');
            
            lag_lim = 7;
            num_nodes = size(this.avg_gamma,2);
            data = this.avg_gamma;
            data1 = nanmean(data,2);
            data2 = data;
            data1 = data1 - nanmean(data1);
            data2 = bsxfun(@minus,data2,nanmean(data2));
            
            ccfs = zeros(1+2*lag_lim*this.Fs_smooth,num_nodes);
            for i = 1:num_nodes
                tic
                ccfs(:,i) = xcorr(data2(:,i),data1,lag_lim*this.Fs_smooth,'coeff');
                toc
            end
            
            % store
            [corrs,inds] = max(ccfs);
            inds(isnan(corrs)) = nan;
            inds = inds-(length(ccfs)-1)/2;
            inds = inds-nanmedian(inds);
            inds = inds/this.Fs_smooth;
            
            figure
            image(I);axis equal
            hold on
            scatter(X,Y,50,inds,'filled','linewidth',2);
            axis off
            
            
            % Plot time series
            figure;hold; box off
            set(gca,'fontsize',13,'fontweight','bold')
            jets = jet(128);
            [a,b] = sort(inds);
            avg_gamma_sort = this.avg_gamma(:,b);
            for n = 1:128
                p1 = plot(-20:1/this.Fs_smooth:20,avg_gamma_sort(:,n),'color',jets(n,:),'linewidth',1);
                p1.Color(4) = .75;
            end
            plot(-20:1/this.Fs_smooth:20,nanmean(this.avg_gamma,2),'k','linewidth',5)
            axis off
            
            obj = scalebar;
            obj.YLen = .1;
            obj.XLen = 5;
            obj.YUnit = 'dB';
            obj.XUnit = 's';
            obj.hTextX_Pos = [1.4 -.012];
            obj.Position = [-18 -.18];
        end
        function save(this)
            fn = fullfile(this.subjectPath, strcat(clientname(true, 3), '.mat'));
            save(fn, 'this');
        end
        function surf_avg_sst(this)
            figure;
            surf((1/this.Fs_smooth)*(-100:100),this.f,this.avg_sst','EdgeColor','none');
            axis xy;
            axis tight;
            colormap(jet);
            view(0,90);
            colorbar;
            c=colorbar; %c.Label.String = 'Z';
            xlabel('Seconds')
            ylabel('Frequency (Hz)')
            caxis([-1,1])
            set(gca,'fontsize',13,'fontweight','bold')
        end
    end

    %% PROTECTED

    properties (Access = protected)
    end        
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
