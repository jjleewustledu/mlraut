classdef Propagator < handle
    %% line1
    %  line2
    %  
    %  Created 19-Jul-2025 18:16:41 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 25.1.0.2943329 (R2025a) for MACA64.  Copyright 2025 John J. Lee.
    
    methods
        function this = Propagator(varargin)
            %% PROPAGATOR 
            %  Args:
            %      arg1 (its_class): Description of arg1.
            
            ip = inputParser;
            addParameter(ip, "arg1", [], @(x) true)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
