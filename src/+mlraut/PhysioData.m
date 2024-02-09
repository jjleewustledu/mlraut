classdef PhysioData < handle & mlsystem.IHandle
    %% line1
    %  line2
    %  
    %  Created 08-Feb-2024 19:45:15 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    
    methods
        function this = PhysioData(ihcp)
            arguments
                ihcp mlraut.HCP {mustBeNonempty}
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
