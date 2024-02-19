classdef CorticalData < handle & mlraut.NetworkData
    %% line1
    %  line2
    %  
    %  Created 10-Feb-2024 15:11:28 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    

    properties (Dependent)
        anat_mask
    end

    methods %% GET
        function g = get.anat_mask(this)
            if ~isempty(this.anat_mask_)
                g = this.anat_mask_;
                return
            end

            ld = load(fullfile(this.waves_dir, 'supporting_files', 'mask_ctx_HCP.mat'));
            this.anat_mask_ = ld.mask_ctx;
            g = this.anat_mask_;
        end
    end

    methods
        function this = CorticalData(varargin)
            this = this@mlraut.NetworkData(varargin{:});
        end
    end

    %% PRIVATE

    properties (Access = private)
        anat_mask_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
