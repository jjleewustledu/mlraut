classdef Gifti < handle & mlsystem.IHandle
    %% line1
    %  line2
    %  
    %  Created 10-Feb-2024 13:47:30 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    
    methods
        function g = aparc_a2009s_label_gii(this, sub, hemi)
            %  e.g.:
            %  cd('/Volumes/PrecunealSSD2/HCP/AWS/hcp-openaccess/HCP_1200/100307/MNINonLinear/fsaverage_LR32k')
            %  g = gifti('100307.L.aparc.a2009s.32k_fs_LR.label.gii')
            
            arguments
                this mlraut.Gifti
                sub {mustBeTextScalar} = this.ihcp_.current_subject
                hemi {mustBeTextScalar} = 'L'
            end
            if contains(hemi, 'L', IgnoreCase=true)
                hemi = 'L';
            end
            if contains(hemi, 'R', IgnoreCase=true)
                hemi = 'R';
            end

            pth = fullfile(this.ihcp_.root_dir, sub, 'MNINonLinear', 'fsaverage_LR32k');
            fn = sprintf('%s.%s.aparc.a2009s.32k_fs_LR.label.gii', sub, hemi);
            fn = fullfile(pth, fn);
            g = gifti(convertStringsToChars(fn));
        end

        function this = Gifti(ihcp)
            arguments
                ihcp mlraut.HCP
            end

            this.ihcp_ = ihcp;
        end
    end

    %% PROTECTED

    properties (Access = protected)
        ihcp_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
