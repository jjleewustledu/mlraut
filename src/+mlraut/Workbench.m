classdef Workbench < handle & mlsystem.IHandle
    %% line1
    %  line2
    %  
    %  Created 10-Feb-2024 13:47:38 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    
    methods
        function this = Workbench(varargin)
            %% WORKBENCH 
            %  Args:
            %      arg1 (its_class): Description of arg1.
            
            this = this@mlsystem.IHandle(varargin{:});
            
            ip = inputParser;
            addParameter(ip, "arg1", [], @(x) true)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end