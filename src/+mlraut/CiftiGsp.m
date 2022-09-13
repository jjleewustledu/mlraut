classdef CiftiGsp 
	%% CIFTIGSP  
    %  uses cifti implementations from Tim Coalson:
    %  https://github.com/Washington-University/cifti-matlab

	%  $Revision$
 	%  was created 26-Mar-2021 14:26:01 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
 	%% It was developed on Matlab 9.9.0.1592791 (R2020b) Update 5 for MACI64.  Copyright 2021 John Joowon Lee.
 	
	properties
        folder_4k = 'cifti_timeseries_normalwall_atlas_32k_2016SVD'
        Nsubjects = 1
        %num_nodes = 59412; % num cortical vertices
        num_nodes = 69206; % num grayordinates for GSP data
        outdir
        rootdir
        ryansdir
        subjects
        wb_dir
    end
    
    properties (Dependent)
        ciftiPrototype
        hcpPrototype
        GSP_32k
    end

	methods     
        function g = get.ciftiPrototype(this)
            g = cifti_read( ...
                    fullfile(this.GSP_32k, this.subjects{1}, this.folder_4k, ...
                             [this.subjects{1} '_LR_surf_32k_fs_LR_smooth1.6985_subcort_smooth1.6985.dtseries.nii']));
        end
        function g = get.GSP_32k(this)
            g = fullfile(this.ryansdir, '32k');
        end
        function g = get.hcpPrototype(~)
            g = cifti_read( ...
                    fullfile(getenv('REFDIR'), ...
                             'test_ryan_91k.dtseries.nii'));
        end
        function ciftis = patricksArgmax(this, carr)
            assert(iscell(carr))
            
            ciftis = cell(1,5); % 2:6 cluster analyses
            for clusti = 1:5
                cifti = this.ciftiPrototype;
                cifti.diminfo{2}.length = 1;
                cifti.cdata = ascol(carr{clusti,3});
                cifti_write(cifti, sprintf('patricksArgmax_%iclusters.dtseries.nii', clusti+1))
                ciftis{clusti} = cifti;
            end
        end
        function ciftis = template(this)
            
            % Subjects            
            for s = 1:this.Nsubjects
                tic
                
                subj = this.subjects{s};
                disp(['Processing subject ' num2str(s) ': ' subj]);
                
                % load BOLD time series
                ciftis{s} = cifti_read( ...
                    fullfile(this.GSP_32k, subj, this.folder_4k, ...
                             [subj '_LR_surf_32k_fs_LR_smooth1.6985_subcort_smooth1.6985.dtseries.nii'])); %#ok<AGROW>
                
                % do whatever you want below to include/exclude certain brain regions
%                mask = cifti.brainstructure(cifti.brainstructure>0); % exclude medial wall vertices, which have no data
                
%                mask1 = mask<3; %cortex
%                BOLD = cifti.cdata(mask1,:)';
                
%                motion_file = [root_dir '/' subj '/FCmaps/' subj '_faln_dbnd_xr3d_atl_g7_bpss_resid.vals'];
%                motion = dlmread(motion_file);
%                format = motion<=motion_thresh; 
            
%                cifti.cdata = zeros(length(mask),1);
%                cifti.hdr.dim(6) = size(BOLD,1); %must equal # columns (number of "timepoints" or maps in cifti file)
%                cifti.cdata(mask1,:) = BOLD'; % change to spacextime before putting into cifti file
%                cifti_write(fullfile(this.outdir, 'test.dtseries.nii'), cifti);
                % there are multiple cifti filetypes (e.g. dscalar, dtseries) but 
                % I would recommend just using dtseries for everything
                
                toc
            
            end
        end
        
        function ciftis = writeClusters(this, carr, fileprefix)
            assert(iscell(carr))
            assert(ischar(fileprefix))
            
            ciftis = cell(1,5); % 2:6 cluster analyses
            for clusti = 1:5
                cifti = this.hcpPrototype; %this.ciftiPrototype;
                cifti.diminfo{2}.length = 1;
                arr = carr{clusti, 3};
                cifti.cdata = ascol(arr(1,:));
                cifti_write(cifti, sprintf('%s_%iclusters.dtseries.nii', fileprefix, clusti+1))
                ciftis{clusti} = cifti;
            end
        end
		  
 		function this = CiftiGsp(varargin)
            ip = inputParser;
            addParameter(ip, 'outdir', fullfile(getenv('HOME'), 'Tmp', ''), @ischar)
            switch computer
                case 'GLNXA64'
                    addParameter(ip, 'rootdir', '/data/nil-bluearc/shimony/GSP', @isfolder)
                    addParameter(ip, 'ryansdir', '/data/nil-bluearc/raichle/ryan/GSP', @isfolder)
                    addParameter(ip, 'wb_dir', '/data/nil-bluearc/raichle/ryan/workbench/bin_linux64', @isfolder)
                case 'MACI64'
                    addParameter(ip, 'rootdir', fullfile(getenv('HOME'), 'Google Drive', 'GSP'), @isfolder)
                    addParameter(ip, 'ryansdir', fullfile(getenv('HOME'), 'Google Drive', 'GSP'), @isfolder)
                    addParameter(ip, 'wb_dir', fullfile('Applications', 'workbench', 'bin_macosx64', ''), @isfolder)
                otherwise
                    error('mlraut:ValueError', 'CiftiGSP.ctor')
            end
            parse(ip, varargin{:})
            ipr = ip.Results;
            this.outdir = ipr.outdir;            
            if ~isfolder(this.outdir)
                mkdir(this.outdir)
            end
            this.rootdir = ipr.rootdir;
            this.ryansdir = ipr.ryansdir;
            this.wb_dir = ipr.wb_dir;
            this.subjects = textread(fullfile(this.GSP_32k, 'patids_2BOLDs.lst'), '%s'); %#ok<DTXTRD>
            
 			addpath(genpath(fullfile(getenv('HOME'), 'MATLAB-Drive', 'cifti-matlab-master')))
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end

