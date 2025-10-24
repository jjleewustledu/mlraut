classdef VortexSurrogates < handle & mlraut.VortexData
    %% Surrogate data supports the null model.
    %  
    %  Created 20-Sep-2025 20:34:41 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 25.1.0.2973910 (R2025a) Update 1 for MACA64.  Copyright 2025 John J. Lee.
    
    methods
        function this = VortexSurrogates(varargin)
            this = this@mlraut.VortexData(varargin{:})
        end

        function DataOut = preprocessing_main(this, subject,hemisphere,flagSur,flagRest,flagTask,flagSmooth)
            %% Highest level of preprocessing.

            assert(flagSur)

            main_folder = this.home;

            % set basic parameters
            if flagSmooth == 2 || flagSmooth == 0  % no need to do smoothed surrogate data as it is done already when flagSmooth == 1
                DataOut = [];
                return
            end

            % load raw and position data files
            for iSub = subject:subject
                tic
                subid = this.subs(iSub);
                if flagRest == 1 % resting data
                    dataDir = this.dtseries(subid=subid);
                    posFile = this.surf(subid=subid, hemi=hemisphere);
                else
                    error("mlraut:ValueError", stackstr())
                end
                
                % process the data files
                cd(main_folder)
                surMethodNum = 7 ; flagVisBpSig = 0 ;
                DataOut = this.load_fMRI( ...
                    dataDir,posFile,surMethodNum,flagVisBpSig,flagRest,flagSmooth,hemisphere,flagTask);

                if flagRest == 1
                    if flagSmooth == 1  % temporal and bandpass filtered (both smooth and unsmoothed surrogate data in one go to ensure consistency within the same batch of randomization)
                        if hemisphere == 1 % left hemisphere
                            folder_name = [main_folder,'/Sample Data/Resting/Preprocessed Data'];
                            cd(folder_name)
                            DataOut_unsmooth = DataOut(:,:,:,2);  % unsmoothed
                            save(['Preprocessed_temporalbandpass_data_sur_resting_LEFT_sub',num2str(iSub),'.mat'],'DataOut_unsmooth')
                            DataOut_smooth = DataOut(:,:,:,1);  % smoothed
                            save(['Preprocessed_spatiotemporalbandpass_data_sur_resting_LEFT_sub',num2str(iSub),'.mat'],'DataOut_smooth')
                        elseif hemisphere == 2 % right hemisphere
                            folder_name = [main_folder,'/Sample Data/Resting/Preprocessed Data'];
                            cd(folder_name)
                            DataOut_unsmooth = DataOut(:,:,:,2);  % unsmoothed
                            save(['Preprocessed_temporalbandpass_data_sur_resting_RIGHT_sub',num2str(iSub),'.mat'],'DataOut_unsmooth')
                            DataOut_smooth = DataOut(:,:,:,1);  % smoothed
                            save(['Preprocessed_spatiotemporalbandpass_data_sur_resting_RIGHT_sub',num2str(iSub),'.mat'],'DataOut_smooth')
                        end
                    end
                else
                    error("mlraut:ValueError", stackstr())
                end

                disp(['finishing subject ', num2str(iSub)])
                toc
                cd(main_folder)
            end
        end

        function [dataOut, b, bandpass_reInterp,data_reInterp,posValid,nanChans,mask] = ...
                load_fMRI(...
                this, dataDir,posFile,surMethodNum,flagVisBpSig,flagRest,flagSmooth,hemisphere,flagTask)

            %% load_fMRI performs preprocessing immediately after loading.
            %  N.B.:  temporal smoothing is permanently disabled and must be performed prior to saving data.

            dataDir = char(dataDir);
            posFile = char(posFile);

            params_ = this.params(hemisphere);
            downSRate = params_.downSRate ;
            xCord = params_.xCord ;
            yCord = params_.yCord ;

            %% Import cortex data

            tic
            % addpath(genpath([pwd,'/ToolOthers/fMRI/cifti-matlab-master']))
            % load the surface data
            data = ft_read_cifti(dataDir) ;
            % load the position data
            posLC = gifti(posFile) ;

            fsTem = 1/0.72 ;
            % this function extract the cortex part of data from the HCP data
            [sigValid,posValid,nanChans] = preproc_fRMI(data,posLC,fsTem) ;
            toc
            % clearvars data posLC posRC

            % checkRegion(posValid)

            %% new interpolation based on the whole left cortex

            x = double(posValid.vertices(:,1));
            y = double(posValid.vertices(:,2));
            % k = convhull(x,y);
            k = alphaShape(x,y,4) ;
            [a, b] = k.boundaryFacets();

            bw = poly2mask(b(:,1)-min(xCord)+1,b(:,2)-min(yCord)+1,max(yCord)-min(yCord)+1,...
                max(xCord)-min(xCord)+1);

            flagPlotSpe = 0 ;
            flagCheckIntp = 0 ;

            mask = double(bw(1:downSRate:end,1:downSRate:end));
            mask(mask==0) = nan ;
            warning('off', 'MATLAB:griddata:DuplicateDataPoints');
            data_reInterp = spaceFreq_fMRI(posValid,sigValid,xCord,yCord,...
                flagPlotSpe, flagCheckIntp,mask);
            warning('on', 'MATLAB:griddata:DuplicateDataPoints');
            dataOut = data_reInterp;
            clearvars a bw flagPlotSpe flagCheckIntp

            %% temporal bandpass filtering
            if flagSmooth == 1 || flagSmooth == 2
                close all
                dataIn = dataOut;
                fLow = 0.01 ; % lower limit of the bandpass filter
                fHigh = 0.1 ; % upper limit of the bandpass filter
                dataOut_reshape = reshape(dataIn,size(dataIn,1)*size(dataIn,2),[]) ;
                [bandpasSig,~, ~, ~] = bandpa_fMRI(dataOut_reshape,fsTem,fLow,fHigh) ;
                %     bandpasSig = zscore(bandpasSig,[],2) ; % ***************
                bandpass_reInterp = reshape(bandpasSig,size(dataIn,1),size(dataIn,2),[]) ;
                dataOut = bandpass_reInterp;
            end

            %% generate surrogate data

            if flagRest == 0 && flagTask ~= 3 % task data
                dataIn = dataOut(1:175,1:251,1:315); % ensure odd number
            elseif flagRest == 1 % rest data
                dataIn = dataOut(1:175,1:251,1:this.Nt);
            end
            if flagRest == 0 && flagTask == 3
                dataIn = dataOut(1:175,1:251,1:405); % ensure odd number
            end

            surData = zeros(size(dataIn)) ;
            dataIn(isnan(dataIn)) = 0 ;
            midPoint = floor(size(dataIn)/2)+1 ;
            indMid = sub2ind(size(dataIn),midPoint(1),midPoint(2),midPoint(3)) ;
            randSeq = randperm(indMid - 1) ;
            freqData = fftshift(fftn(dataIn)) ; % 3D fourier transform
            phaseData = angle(freqData) ;
            absData = abs(freqData) ;
            phaseDataRand = phaseData;
            phaseDataRand(1:indMid-1) = rand(indMid-1,1).*2*pi-pi; % random phase values -pi~pi
            for iInd = 1:indMid - 2
                phaseDataRand(indMid+iInd) = -phaseDataRand(indMid-iInd) ;
            end
            freqDataNew = ifftshift(absData.*exp(1i.*phaseData + 1i*phaseDataRand)) ; % random phase shuflling
            dataOut = real(ifftn(freqDataNew)) ; % 3D inverse fourier tranform
            dataOut_unsmooth = dataOut;

            %% spatial bandpass filtering - Gaussian

            if flagSmooth == 1
                dataIn = dataOut(1:175,:,:);
                sigmScale = params_.sigmScale/downSRate ;
                sigLPass = [] ;
                sigBPass = [] ;
                % spatial bandpass filter
                numScale = length(sigmScale)-1 ;
                for iTime = 1:size(dataIn,3)
                    filtSigma = sigmScale(1);   % 0.6
                    filtWidth = ceil(3*filtSigma);
                    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
                    sigLPass(1,:,:,iTime) = nanconv(dataIn(:,:,iTime),imageFilter,'edge', 'nanout');
                    for iScale = 1:numScale
                        filtSigma = sigmScale(iScale+1);
                        filtWidth = ceil(3*filtSigma);
                        imageFilter=fspecial('gaussian',filtWidth,filtSigma);
                        sigLPass(iScale+1,:,:,iTime) = nanconv(dataIn(:,:,iTime),imageFilter,'edge', 'nanout');
                        sigBPass(iScale,:,:,iTime) = sigLPass(iScale+1,:,:,iTime) - sigLPass(iScale,:,:,iTime) ;
                    end
                end

                dataOut = permute(sigBPass,[2,3,4,1]);
            end

            %% ensure the edge of the output data is compatible with the flattenedd cortex

            if hemisphere == 1
                load('parcellation_template.mat')
            elseif hemisphere == 2
                load('parcellation_template22_RightBrain_subject1-100.mat')
                parcellation_template = parcellation_template22_RightBrain_100sub(:,:,1);
            end
            dataOut = dataOut(1:175,:,:).*parcellation_template(1:175,:)./parcellation_template(1:175,:);
            
            dataOut_unsmooth = dataOut_unsmooth.*parcellation_template(1:175,:)./parcellation_template(1:175,:);
            dataOut(:,:,:,1) = dataOut;
            dataOut(:,:,:,2) = dataOut_unsmooth; % record both smooth and unsmoothed surrogate data in the same output
        end
    end

    %% PRIVATE

    properties (Access = private)
        vortex_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
