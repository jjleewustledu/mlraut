classdef VortexData < handle & mlsystem.IHandle
    %% Data that is not surrogate data.
    %  
    %  Created 20-Sep-2025 20:37:46 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 25.1.0.2973910 (R2025a) Update 1 for MACA64.  Copyright 2025 John J. Lee.
    

    properties (Dependent)
        home
        Nt
        subs
    end

    methods  % GET/SET
        function g = get.home(this)
            g = this.vortex_.home;
        end

        function g = get.Nt(this)
            g = this.vortex_.Nt;
        end

        function g = get.subs(this)
            g = this.vortex_.subs;
        end
    end

    methods
        function this = VortexData(vortex)
            this.vortex_ = vortex;
        end

        function dts = dtseries(this, varargin)
            dts = this.vortex_.dtseries(varargin{:});
        end

        function s = surf(this, varargin)
            s = this.vortex_.surf(varargin{:});
        end

        function p = params(~, hemisphere)
            % set basic parameters
            p.sigmScale =  [29.35 14.93];           % bandpass 5 bandwidth ranges
            p.downSRate = 2 ;                       % downsample the re-interpolation
            if hemisphere == 1
                p.xCord = -250:p.downSRate:250 ;    % coordinate re-interpolation, left hemisphere
                p.yCord = -150:p.downSRate:200 ;
            elseif hemisphere == 2
                p.xCord = -270:p.downSRate:230 ;    % coordinate re-interpolation, right hemisphere
                p.yCord = -180:p.downSRate:170 ;
            end
            p.fsTem = 1/0.72 ;                      % temporal sampling rate
        end

        function DataOut = preprocessing_main(this, subject,hemisphere,flagSur,flagRest,flagTask,flagSmooth)
            %% Highest level of preprocessing.

            assert(~flagSur)

            main_folder = this.home;

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
                    if flagSmooth == 1 % temporal bandpass filtered
                        if hemisphere == 1 % left hemisphere
                            folder_name = [main_folder,'/Sample Data/Resting/Preprocessed Data'];
                            cd(folder_name)
                            save(['Preprocessed_temporalbandpass_data_resting_LEFT_sub',num2str(iSub),'.mat'],'DataOut')
                        elseif hemisphere == 2 % right hemisphere
                            folder_name = [main_folder,'/Sample Data/Resting/Preprocessed Data'];
                            cd(folder_name)
                            save(['Preprocessed_temporalbandpass_data_resting_RIGHT_sub',num2str(iSub),'.mat'],'DataOut')
                        end
                    elseif flagSmooth == 0 % raw data
                        if hemisphere == 1 % left hemisphere
                            folder_name = [main_folder,'/Sample Data/Resting/Preprocessed Data'];
                            cd(folder_name)
                            save(['Preprocessed_raw_data_resting_LEFT_sub',num2str(iSub),'.mat'],'DataOut')
                        elseif hemisphere == 2 % right hemisphere
                            folder_name = [main_folder,'/Sample Data/Resting/Preprocessed Data'];
                            cd(folder_name)
                            save(['Preprocessed_raw_data_resting_RIGHT_sub',num2str(iSub),'.mat'],'DataOut')
                        end
                    elseif flagSmooth == 2 % spatiotemporal bandpss filtered data
                        if hemisphere == 1 % left hemisphere
                            folder_name = [main_folder,'/Sample Data/Resting/Preprocessed Data'];
                            cd(folder_name)
                            save(['Preprocessed_spatiotemporalbandpass_data_resting_LEFT_sub',num2str(iSub),'.mat'],'DataOut')
                        elseif hemisphere == 2 % right hemisphere
                            folder_name = [main_folder,'/Sample Data/Resting/Preprocessed Data'];
                            cd(folder_name)
                            save(['Preprocessed_spatiotemporalbandpass_data_resting_RIGHT_sub',num2str(iSub),'.mat'],'DataOut')
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

            %% spatial bandpass filtering - Gaussian

            if flagSmooth == 2
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
        end
    end

    %% PRIVATE

    properties (Access = private)
        vortex_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
