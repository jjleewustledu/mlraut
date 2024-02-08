classdef AnalyticSignalPar < handle & mlraut.AnalyticSignal
    %% line1
    %  line2
    %  
    %  Created 13-Apr-2023 02:11:52 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 9.14.0.2206163 (R2023a) for MACI64.  Copyright 2023 John J. Lee.
    
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
        function this = AnalyticSignalPar(varargin)
            this = this@mlraut.AnalyticSignal(varargin{:});
        end
        
        function this = call(this)
            %% CALL all subjects

            % exclude subjects
            this.subjects = this.subjects(~contains(this.subjects, '_7T'));
            %this.subjects = this.subjects(~contains(this.subjects, 'sub-'));

            out_dir_ = this.out_dir;
            for s = 1:this.num_sub
                this.current_subject = this.subjects{s};
                if ~contains(out_dir_, this.current_subject)
                    proposed_dir = fullfile(out_dir_, this.current_subject);
                    ensuredir(proposed_dir);
                    this.out_dir = proposed_dir;
                end 
                this.call_subject(s);
            end
        end
        function this = call_subject(this, s)
            arguments
                this mlraut.AnalyticSignal
                s double
            end

            for t = 1:this.num_tasks
                this.current_task = this.tasks{t};   

                % BOLD
                try
                    bold = this.task_dtseries(); 
                    assert(~isempty(bold))
                    bold = this.omit_late_frames(bold);
                catch ME
                    disp([this.current_subject ' ' this.current_task ' BOLD missing or defective:']);
                    handwarning(ME)
                    continue
                end

                % Global signal
                gs = this.global_signal(bold);
                 
                % Physio
                try
                    physio = this.task_physio(bold);
                    assert(~isempty(physio))
                catch ME
                    disp([this.current_subject ' ' this.current_task ' physio missing or defective:']);
                    handwarning(ME)
                    continue
                end

                % Analytic signal
                bold_ = this.center_and_rescale(this.band_pass(bold - gs));
                physio_ = this.center_and_rescale(this.band_pass(physio)); % removes gs as needed
                bold_ = hilbert(bold_);
                physio_ = hilbert(physio_);
                as = conj(physio_).*bold_; % <psi_p|BOLD_operator|psi_p> ~ <psi_p|psi_b>, not unitary
                as = this.normalize_all(as);
        
                % Store reduced analytic signal, real(), imag(), abs(), angle()
                save(fullfile(this.out_dir, sprintf('%s_bold-gs%s_%i_%i', stackstr(2), this.tags, s, t)), 'bold_');
                save(fullfile(this.out_dir, sprintf('%s_as%s_%i_%i', stackstr(2), this.tags, s, t)), 'as');
            end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
