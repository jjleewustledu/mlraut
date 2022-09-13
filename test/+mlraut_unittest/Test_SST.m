classdef Test_SST < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 04-Apr-2022 16:22:07 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/test/+mlraut_unittest.
    %  Developed on Matlab 9.12.0.1884302 (R2022a) for MACI64.  Copyright 2022 John J. Lee.
    
    properties
        ryans_ecog_dir = '/Volumes/raichle/ryan/ecog'
        testObj
    end
    
    methods (Test)
        function test_afun(this)
            import mlraut.*
            this.assumeEqual(1,1);
            this.verifyEqual(1,1);
            this.assertEqual(1,1);
        end
    end
    
    methods (TestClassSetup)
        function setupSST(this)
            import mlraut.*
            this.testObj_ = SST();
        end
    end
    
    methods (TestMethodSetup)
        function setupSSTTest(this)
            this.testObj = this.testObj_;
            this.addTeardown(@this.cleanTestMethod)
        end
    end
    
    properties (Access = private)
        testObj_
    end
    
    methods (Access = private)
        function cleanTestMethod(this)
        end
    end

    methods
        function monkey_FC(this, monkey, cond)
            
            %% Inputs

            %monkey = 'Chibi'; % Chibi, George, Kin2, Su
            %cond = 2; % 1 = EO; 2 = EC; 3 = sleep; 4 = anesth
            
            %% Setup

            root_dir = fullfile(this.ryans_ecog_dir, monkey, '');
            cd(root_dir)
            load(fullfile(root_dir, [monkey 'Map.mat']));
            
            switch cond
                case 1
                    label1 = 'EO';
                    label2 = 'AwakeEyesOpened';
                case 2
                    label1 = 'EC';
                    label2 = 'AwakeEyesClosed';
                case 3
                    label1 = 'sleep';
                    label2 = 'Sleeping';
                case 4
                    label1 = 'anesth';
                    label2 = 'Anesthetized';
            end
            
            runs = textread([root_dir '/../' monkey '_' label1 '_runs.txt'],'%s');
            
            num_nodes = 128;
            Fs = 1000;
            Fs_dec = 200;
            
            corr_mats = nan(num_nodes,num_nodes,numel(runs));
            
            for r = 1:length(runs)
                disp(['Processing ' runs{r}])
                disp([num2str(r) ' out of ' num2str(numel(runs))])
                data_dir = [root_dir '/../' runs{r}];
                time_vec = importdata([data_dir '/ECoGTime.mat']);
                info = importdata([data_dir '/Condition.mat']);
                start = info.ConditionTime(strcmp(info.ConditionLabel,[label2 '-Start'])); % in seconds
                stop = info.ConditionTime(strcmp(info.ConditionLabel,[label2 '-End']));
            
                start = find(time_vec==start);
                if isempty(start)
                    start = 10*Fs; % start ten seconds in
                end
                stop = find(time_vec==stop);
                num_samples = ceil(length(start:stop)/5);
                
                %format = true(length(start:stop),1);
            
                %% Load data

                disp('Loading data...')
                raw_data = single(zeros(num_samples,num_nodes));
                for i = 1:num_nodes
                    tic
                    temp = importdata([data_dir '/ECoG_ch' num2str(i) '.mat']);
                    raw_data(:,i) = decimate(temp(start:stop),5);
                    toc
                end
                
                %% Bad channels

                good = true(num_nodes,1);                
                switch monkey  
                    case 'George'
                        good([53,73]) = false;
            
                    case 'Su'
                        good(50) = false;
                end
            
                %% Time-frequency QC    
                 
            %     [pxx,f] = pwelch(raw_data,window,[],[],Fs);
            %     figure;loglog(f,pxx(:,good))
            
                raw_data(:,~good) = nan;
            
                %% Filter
            
                disp('Filtering...')
                Fs = Fs_dec;
                [b1,a1] = butter(2,[49,51]/(Fs/2),'stop'); % notch
                [b2,a2] = butter(1,[40,99]/(Fs/2)); % gamma
                [b3,a3] = butter(1,[.01,.1]/(Fs/2)); % infra-slow envelope of gamma blp
                
                filt_data = filtfilt(b1,a1,double(raw_data));
                filt_data = filtfilt(b2,a2,filt_data);
                filt_data = filtfilt(b3,a3,abs(hilbert(filt_data))); % hilbert envelope
                
                format = true(size(filt_data,1),1);
                format(1:10000) = false;
                corr_mats(:,:,r) = corr(filt_data(format,:)); % 128 x 128 x 8; size(filt_data,2) == 128                
            end
            
            %save([outdir '/' monkey '_' label1 '_FC_mats.mat'],corr_mats,'-v7.3')            
            
            %% Diffusion embedding
            load(fullfile(root_dir, [monkey 'Map.mat']));
            %load('GeorgeMap.mat');
            FCmat = tanh(nanmean(atanh(corr_mats),3)); % Fisher z transform before averaging
            
            alpha = .5;
            good = ~isnan(FCmat(:,1));
            FCmat = FCmat(good,good);
            FCmat = bsxfun(@rdivide,FCmat,sqrt(sum(FCmat.^2,2)));
            FCmat = FCmat*FCmat';
            %FCmat(FCmat<0) = 0;
            %FCmat(FCmat>1) = 1;
            FCmat = 1-acos(FCmat)/pi;
            L = FCmat;
            D = sum(L,2).^-alpha;
            Lalpha = L.*(D*D'); % normalized graph Laplacian
            Dalpha = sum(Lalpha,2);
            M = bsxfun(@rdivide,Lalpha,Dalpha);
            [maps,S,~] = svd(M);
            
            map = nan(128,1);
            map(good) = maps(:,2);
            
            figure
            image(I);axis equal
            hold on
            scatter(X,Y,50,map,'filled','linewidth',2); % 50 window, 300 full screen
            axis off
            
            colormap(magma) 
        end
        function monkey_plotting(this)
        end
        function monkey_SSTs(this)
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
