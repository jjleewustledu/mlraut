classdef AnalyticSignalPar < handle & mlraut.AnalyticSignal
    %% line1
    %  line2
    %  
    %  Created 13-Apr-2023 02:11:52 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 9.14.0.2206163 (R2023a) for MACI64.  Copyright 2023 John J. Lee.
    
    methods (Static)
        function parcall(cores, opts)
            arguments
                cores {mustBeScalarOrEmpty} = 4
                opts.N_sub {mustBeScalarOrEmpty} = 1113
                opts.flip_globbed logical = true
            end

            root_dir = '/home/usr/jjlee/mnt/CHPC_hcpdb/packages/unzip/HCP_1200';
            %root_dir = '/home/usr/jjlee/mnt/CHPC_scratch/Singularity/HcpAging/HCPAgingRec/fmriresults01';
            out_dir = '/home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignal';
            %out_dir = '/vgpool02/data2/jjlee/AnalyticSignalHcpAging';

            g = glob(fullfile(root_dir, '*'));
            if opts.flip_globbed
                g = flip(g); % examine more recent ones
            end
            g = cellfun(@(x) basename(x), g, UniformOutput=false);
            g = g(~contains(g, 'manifests'));
            g = g(1:opts.N_sub);
            leng = length(g);
            %for idxg = 1:1
            %parfor (idxg = 1:2, 2)
            parfor (idxg = 1:leng, cores)
                try
                    if isfolder(fullfile(out_dir, g(idxg)))
                        continue
                    end
                    this = mlraut.AnalyticSignalPar( ...
                        subjects=g(idxg), ...                     
                        do_save=false, ...
                        do_save_ciftis=true, ...
                        do_save_ciftis_of_diffs=false, ...
                        tags="AnalyticSignalPar-parcall");
                    call(this);
                catch ME
                    handwarning(ME)
                end
            end
        end
    end

    methods
        function this = AnalyticSignalPar(varargin)
            this = this@mlraut.AnalyticSignal(varargin{:});
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
