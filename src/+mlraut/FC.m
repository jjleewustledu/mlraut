classdef FC
    %% line1
    %  line2
    %  
    %  Created 04-Apr-2022 15:33:38 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 9.12.0.1884302 (R2022a) for MACI64.  Copyright 2022 John J. Lee.

    methods (Static)
        function createMap(electrodes_mat)
            assert(isfile(electrodes_mat))

            [elepth,elefp] = myfileparts(electrodes_mat);
            elefp = strrep(elefp, ' ', '_');
            pwd0 = pushd(elepth);
            ld = load(electrodes_mat);
            N = length(ld.electrodes);
            msk = mlfourd.ImagingContext2(fullfile(getenv('REFDIR'), 'colin27_t1_tal_lin_mask.nii'));
            tal = mlfourd.ImagingContext2(fullfile(getenv('REFDIR'), 'colin27_t1_tal_lin.nii'));
            tal = tal .* msk;
            tal.filepath = elepth;
            tal1 = tal.thresh(2382000);
            tal1 = tal1.uthresh(4140000);
            origin = tal.imagingFormat.originator - 1;
            for i = 1:N
                x(i) = electrodes(i).x + origin(1); 
            end
            for j = 1:N
                y(j) = electrodes(j).y + origin(2); 
            end
            for k = 1:N
                z(k) = electrodes(k).z + origin(3); 
            end
            grid(:,1) = x;
            grid(:,2) = y;
            grid(:,3) = z;
            pcgrid = pointCloud(grid);

            ele = copy(tal);
            ele.setPointCloud(pcgrid); % dipmax ~ 1
            ele = flip(ele, 1);
            ele.fileprefix = elefp;
            ele.save();
            ele.show(tal);

            figure; 
            pcshow(tal1);
            hold on
            pcshow(pcgrid.Location, '*m', 'MarkerSize', 50);
            hold off

            popd(pwd0);
        end
    end

    properties
        avg_alpha
        avg_delta
        avg_gamma
        avg_sst
        decimate_r = 1;
        f
        Fs_desired
        Fs_smooth = 5
        null_sst
        sst_count
        subjectPath        
    end

    properties (Dependent)
    end

    methods
        function this = FC(varargin)
            %% FC 
            %  Args:
            %      subjectPath (folder)
            %      runsPatt (text): e.g., '*_desc-clinical_collection_*' in subjectPath.
            %      Fs_desired (scalar): frequency of sampling desired
            
            ip = inputParser;
            addParameter(ip, "subjectPath", pwd, @isfolder)
            addParameter(ip, "runsPatt", '*', @istext)
            addParameter(ip, "Fs_desired", 200, @isscalar)
            
            parse(ip, varargin{:})
            ipr = ip.Results;
            this.subjectPath = ipr.subjectPath;
            if ~startsWith(ipr.runsPatt, this.subjectPath)
                ipr.runsPatt = fullfile(this.subjectPath, ipr.runsPatt);
            end
            this.Fs_desired = ipr.Fs_desired;

            %% Setup
            addpath(genpath(fullfile( ...
                getenv('HOME'), 'MATLAB-Drive', 'arousal-waves-main', 'Dependencies', 'chronux_2_12')));
            runs = glob(ipr.runsPatt);

            pwd0 = pushd(this.subjectPath);

            corr_mats = cell(1, numel(runs));            
            for r = 1:length(runs)
                disp(['Processing ' runs{r}])
                disp([num2str(r) ' out of ' num2str(numel(runs))])

                ld = load(fullfile(runs{r}, 'hdr.mat'));            
                num_nodes = ld.hdr.ns;
                Fs = ld.hdr.frequency(1);
                if Fs > this.Fs_desired
                    this.decimate_r = ceil(Fs/this.Fs_desired);
                    Fs = this.Fs_desired;
                end
                time_vec = importdata(fullfile(runs{r}, 'ECoGTime.mat'));

                start = 10*Fs; % int index, starting ten seconds in
                stop = length(time_vec); % int index, to end of data collection
                num_samples = ceil(length(start:stop)/this.decimate_r);
            
                %% Load data

                disp('Loading data...')
                raw_data = single(zeros(num_samples,num_nodes));
                for i = 1:num_nodes
                    tic
                    temp = importdata(fullfile(runs{4}, strcat('ECoG_ch', num2str(i), '.mat')));
                    raw_data(:,i) = decimate(temp(start:stop), this.decimate_r);
                    toc
                end
                
                good = true(num_nodes,1);
                if isfield(ld.hdr, 'badchannels')
                    good(ld.hdr.badchannels) = false; % e.g., [53, 73]
                end
                raw_data(:,~good) = nan;
            
                %% Filter
            
                disp('Filtering...')
                [b1,a1] = butter(2,[49,51]/(Fs/2),'stop'); % notch
                [b2,a2] = butter(1,[40,99]/(Fs/2)); % gamma
                [b3,a3] = butter(1,[.01,.1]/(Fs/2)); % infra-slow envelope of gamma blp
                
                filt_data = filtfilt(b1,a1,double(raw_data));
                filt_data = filtfilt(b2,a2,filt_data);
                filt_data = filtfilt(b3,a3,abs(hilbert(filt_data))); % hilbert envelope
                
                format = true(size(filt_data,1),1);
                format(1:10*this.Fs_desired) = false;
                corr_mats{r} = corr(filt_data(format,:)); % 128 x 128; size(filt_data,2) == 128                
            end
            
            %save([outdir '/' monkey '_' label1 '_FC_mats.mat'],corr_mats,'-v7.3')            
            
            %% Diffusion embedding
            ld_map = load(fullfile(this.subjectPath, 'Map.mat'));
            corr_mats_R3 = this.cell2mat3(corr_mats);
            FCmat = tanh(nanmean(atanh(corr_mats_R3),3)); % Fisher z transform before averaging
            
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
            image(ld_map.I);axis equal
            hold on
            scatter(ld_map.X,ld_map.Y,50,map,'filled','linewidth',2); % 50 window, 300 full screen
            axis off
            
            colormap(magma) 

            popd(pwd0);
        end

        function mat3 = cell2mat3(cellobj)
            %  Args:
            %      cellobj (cell): e.g., { matrix_1(N,N), matrix_2(N,N), ..., matrix_Nr(N,N) }
            %  Returns:
            %      mat3: e.g., tensor(N,N,Nr)

            N = min(cell2mat(cellfun(@(x) length(x), cellobj)));
            assert(N > 1, 'mlraut:RuntimeError', 'cell2mat3 requires nontrivial cellobj')
            Nr = length(cellobj);
            mat3 = nan(N, N, Nr);
            for r = 1:length(cellobj)
                mat = cellobj{r};
                mat3(:,:,r) = mat(1:N, 1:N);
            end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
