classdef AnalyticSignalHCPAgingPar < handle & mlraut.AnalyticSignalHCPAging
    %% line1
    %  line2
    %  
    %  Created 07-Feb-2024 23:37:28 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    
    methods (Static)
        function parcall(cores, opts)
            arguments
                cores {mustBeScalarOrEmpty} = 96
                opts.N_sub {mustBeScalarOrEmpty} = 725
            end

            %root_dir = '/home/usr/jjlee/mnt/CHPC_hcpdb/packages/unzip/HCP_1200';
            root_dir = '/home/usr/jjlee/mnt/CHPC_scratch/Singularity/HcpAging/HCPAgingRec/fmriresults01';
            out_dir = '/vgpool02/data2/jjlee/AnalyticSignalHcpAging';
            ensuredir(out_dir);
            %tasks = {'rfMRI_REST1_RL','rfMRI_REST2_RL'};

            g = glob(fullfile(root_dir, '*'));
            g = flip(g); % examine more recent ones
            g = cellfun(@(x) basename(x), g, UniformOutput=false);
            g = g(~contains(g, 'manifests'));
            g = g(1:opts.N_sub);
            leng = length(g);
            for idxg = 1:leng
            %parfor (idxg = 1:leng, cores)
                try
                    this = mlraut.AnalyticSignalPar(subjects=g(idxg), ...
                        root_dir=root_dir, out_dir=out_dir);
                    call(this);
                catch ME
                    handwarning(ME)
                end
            end
        end
    end

    methods
        function this = AnalyticSignalHCPAgingPar(varargin)
            this = this@mlraut.AnalyticSignalHCPAging(varargin{:});
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
