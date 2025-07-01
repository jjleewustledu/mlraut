classdef NetworkData < handle & mlsystem.IHandle
    %% line1
    %  line2
    %  
    %  Created 10-Feb-2024 15:04:15 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    

    % properties (Abstract)
    %     anat_mask  % e.g. cbm, ctx, str, thal
    % end

    properties (Constant)
        NETWORKS_YEO_NAMES = ...
            {'visual', 'somatomotor', 'dorsal attention', 'ventral attention', 'limbic', 'frontoparietal', ...
             'default mode', ...
             'task+', 'task-'}
    end

    properties (Dependent)
        networks_Yeo
        num_frames
        waves_dir
    end

    methods %% GET
        function g = get.networks_Yeo(this)
            if ~isempty(this.networks_Yeo_)
                g = this.networks_Yeo_;
                return
            end

            ld = load(fullfile(this.waves_dir, 'supporting_files', 'networks_HCP.mat'));
            this.networks_Yeo_ = ld.assns2;
            g = this.networks_Yeo_;
        end
        function g = get.num_frames(this)
            g = this.ihcp_.num_frames;
        end
        function g = get.waves_dir(this)
            g = this.ihcp_.waves_dir;
        end
    end

    methods
        function signal = build_anat_signal(this, psi)
            %% For each anatomical partition, average over nodes for a given subject and given task.
            %  Args:
            %       psi {mustBeNumeric} : complex analytic signal
            %  Returns:
            %       signals complex ~ this.num_frames

            arguments
                this mlraut.NetworkData
                psi {mustBeNumeric} = this.psi_
            end

            signal = complex(nan(this.num_frames, 1));

            try
                msk = this.anat_mask;  % subclassed by anatomical subclasses
                signal = mean(psi(:,msk), 2, 'omitnan');
            catch ME
                handerror(ME)
            end
        end

        function signals = build_Yeo_signals(this, psi)
            %% For each Yeo subnet, average over nodes for a given subject and given task.
            %  Args:
            %       psi {mustBeNumeric} : complex analytic signal
            %  Returns:
            %       signals complex ~ this.num_frames x Nrsn

            arguments
                this mlraut.NetworkData
                psi {mustBeNumeric} = this.psi_
            end

            Nyeo = sum(~contains(this.NETWORKS_YEO_NAMES, 'task'));
            signals = complex(nan(this.num_frames, Nyeo));

            for n = 1:Nyeo % Yeo's 7x RSNs
                try
                    msk = this.networks_Yeo == n & this.anat_mask;
                    signals(:,n) = mean(psi(:,msk), 2, 'omitnan');
                catch ME
                    handerror(ME)
                end
            end

            task_pos = 1 <= this.networks_Yeo & this.networks_Yeo <= 6;
            task_neg = this.networks_Yeo == 7;

            % task+
            idx = contains(this.NETWORKS_YEO_NAMES, 'task+');
            msk = task_pos & this.anat_mask;
            signals(:,idx) = mean(psi(:,msk), 2, 'omitnan');

            % task-
            idx = contains(this.NETWORKS_YEO_NAMES, 'task-');
            msk = task_neg & this.anat_mask;
            signals(:,idx) = mean(psi(:,msk), 2, 'omitnan');
        end

        function [psi1,shapes] = reshape_by_anatomy(this, psi)
            ld = load(fullfile(this.waves_dir, 'supporting_files', 'mask_ctx_HCP.mat'));
            mask_ctx = ld.mask_ctx;
            ld = load(fullfile(this.waves_dir, 'supporting_files', 'mask_cbm_HCP.mat'));
            mask_cbm = ld.mask_cbm;
            ld = load(fullfile(this.waves_dir, 'supporting_files', 'mask_str_HCP.mat'));
            mask_str = ld.mask_str;
            ld = load(fullfile(this.waves_dir, 'supporting_files', 'mask_thal_HCP.mat'));
            mask_thal = ld.mask_thal;
            indices_rsn = this.networks_Yeo;

            % aufbau
            psi1 = [];
            shapes = struct();
            accum = 0;
            for rsn = 1:7
                msk = logical(mask_ctx .* (indices_rsn == rsn));
                psi1 = [psi1, psi(:, msk)]; %#ok<*AGROW>
                shapes.ctx.("rsn"+rsn) = sum(msk);
                accum = accum + sum(msk);
            end
            shapes.ctx.accum = accum;
            accum = 0;
            for rsn = 1:7
                msk = logical(mask_cbm .* (indices_rsn == rsn));
                psi1 = [psi1, psi(:, msk)]; %#ok<*AGROW>
                shapes.cbm.("rsn"+rsn) = sum(msk);
                accum = accum + sum(msk);
            end
            shapes.cbm.accum = accum;
            accum = 0;
            for rsn = 1:7
                msk = logical(mask_str .* (indices_rsn == rsn));
                psi1 = [psi1, psi(:, msk)]; %#ok<*AGROW>
                shapes.str.("rsn"+rsn) = sum(msk);
                accum = accum + sum(msk);
            end
            shapes.str.accum = accum;
            accum = 0;
            for rsn = 1:7
                msk = logical(mask_thal .* (indices_rsn == rsn));
                psi1 = [psi1, psi(:, msk)]; %#ok<*AGROW>
                shapes.thal.("rsn"+rsn) = sum(msk);
                accum = accum + sum(msk);
            end
            shapes.thal.accum = accum;
        end

        function this = NetworkData(ihcp, psi)
            arguments
                ihcp mlraut.AnalyticSignal {mustBeNonempty}
                psi {mustBeNumeric}
            end

            this.psi_ = psi;
            this.ihcp_ = ihcp;
        end
    end

    %% PROTECTED

    properties (Access = protected)
        psi_
        ihcp_
        networks_Yeo_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
