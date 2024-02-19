classdef AnalyticSignalPar < handle & mlraut.AnalyticSignal
    %% line1
    %  line2
    %  
    %  Created 13-Apr-2023 02:11:52 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 9.14.0.2206163 (R2023a) for MACI64.  Copyright 2023 John J. Lee.
    
    methods (Static)
        function parcall(cores, opts)
            arguments
                cores {mustBeScalarOrEmpty} = 45
                opts.N_sub {mustBeScalarOrEmpty} = 1113
            end

            root_dir = '/home/usr/jjlee/mnt/CHPC_hcpdb/packages/unzip/HCP_1200';
            %root_dir = '/home/usr/jjlee/mnt/CHPC_scratch/Singularity/HcpAging/HCPAgingRec/fmriresults01';
            %out_dir = '/home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignal';
            %out_dir = '/vgpool02/data2/jjlee/AnalyticSignalHcpAging';

            g = glob(fullfile(root_dir, '*'));
            g = flip(g); % examine more recent ones
            g = cellfun(@(x) basename(x), g, UniformOutput=false);
            g = g(~contains(g, 'manifests'));
            g = g(1:opts.N_sub);
            leng = length(g);
            %for idxg = 1:1
            parfor (idxg = 1:leng, cores)
                try
                    this = mlraut.AnalyticSignalPar( ...
                        subjects=g(idxg), ...                     
                        do_save=true, ...
                        do_save_ciftis=false, ...
                        tags=stackstr(use_dashes=true));
                    call(this);
                catch ME
                    handwarning(ME)
                end
            end
        end
    end

    methods
        function save(this, s, t)

            the_out_dir_ = this.out_dir;
            the_tags_ = this.tags;

            % reduce file size
            this.roi = [];
            this.bold_data_ = [];
            this.cohort_data_ = [];
            this.cifti_last_ = [];
            
            save(fullfile(the_out_dir_, ...
                sprintf("sub-%s_ses-%s_%s.mat", this.subjects{s}, strrep(this.tasks{t}, "_", "-"), the_tags_)), ...
                'this');
        end

        function this = AnalyticSignalPar(varargin)
            this = this@mlraut.AnalyticSignal(varargin{:});
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
