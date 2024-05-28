classdef AnalyticSignalHCPAgingPar < handle & mlraut.AnalyticSignalHCPAging
    %% line1
    %  line2
    %  
    %  Created 07-Feb-2024 23:37:28 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    
    methods (Static)
        function parcall(cores, opts)
            arguments
                cores {mustBeScalarOrEmpty} = 8
                opts.N_sub {mustBeScalarOrEmpty} = 725
            end

            % root_dir = '/home/usr/jjlee/mnt/CHPC_hcpdb/packages/unzip/HCP_1200';
            root_dir = fullfile(getenv('SINGULARITY_HOME'), 'HCPAging', 'HCPAgingRec', 'fmriresults01');

            g = glob(fullfile(root_dir, 'HCA*'));
            g = strip(g, filesep);
            g = flip(g); % examine more recent ones first
            g = cellfun(@(x) basename(x), g, UniformOutput=false);
            g = g(1:opts.N_sub);
            leng = length(g);
            %for idxg = 1:1
            parfor (idxg = 1:leng, cores)
                try
                    tic
                    this = mlraut.AnalyticSignalHCPAging( ...
                        subjects=g(idxg), ...
                        tasks={'fMRI_CONCAT_ALL'}, ...
                        do_save=true, ...
                        do_save_ciftis=true, ...
                        tags=stackstr(use_dashes=true));
                    call(this)
                    toc
                catch ME
                    handwarning(ME)
                end
            end
        end
        function [j,c] = parcluster(globbing_mat)
            arguments
                globbing_mat {mustBeFile} = fullfile( ...
                    getenv('SINGULARITY_HOME'), 'HCPAging', 'HCPAgingRec', 'fmriresults01', ...
                    'mlraut_AnalyticSignalHCPAgingPar_globbing.mat')
            end
            ld = load(globbing_mat);
            globbed = asrow(ld.globbed);
            % globbed = globbed(1:3);

            mlraut.CHPC3.propcluster();
            c = parcluster;
            disp(c.AdditionalProperties)
            for g = globbed
                try
                    j = c.batch( ...
                        @mlraut.AnalyticSignalHCPAgingPar.construct_and_call, ...
                        1, ...
                        {'subjects', g(1)}, ...
                        'CurrentFolder', '.', ...
                        'AutoAddClientPath', false);
                catch ME
                    handwarning(ME)
                end
            end
        end
        function duration = construct_and_call(opts)
            arguments
                opts.subjects cell
                opts.tasks cell = {'fMRI_CONCAT_ALL'}
                opts.tags {mustBeTextScalar} = "AnalyticSignalHCPAgingPar"
            end

            % setenv('TMPDIR', '/scratch/jjlee/tmp') % worker nodesk
            setenv('SINGULARITY_HOME', '/scratch/jjlee/Singularity')
            tic
            this = mlraut.AnalyticSignalHCPAgingPar( ...
                subjects=opts.subjects, ...
                tasks=opts.tasks, ...
                do_save=false, ...
                do_save_ciftis=true, ...
                tags=opts.tags);
            call(this);
            duration = toc;
        end
    end

    methods
        function this = AnalyticSignalHCPAgingPar(varargin)
            this = this@mlraut.AnalyticSignalHCPAging(varargin{:});
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
