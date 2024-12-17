classdef Lee2024 < handle
    %% line1
    %  line2
    %  
    %  Created 10-Feb-2024 13:43:26 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    
    methods
        function this = Lee2024(varargin)
        end

        function write_metric_stats(this, metric_lbl)
            arguments
                this mlraut.Lee2024 %#ok<INUSA>
                metric_lbl {mustBeTextScalar} = "abs_as_"
            end

            out_dir = '/home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalHCPAging';
            %out_dir = '/home/usr/jjlee/mnt/CHPC_scratch/Singularity/AnalyticSignalHCP';
            cd(out_dir);

            g = glob(fullfile(out_dir, 'HCA*', sprintf('%s*_avgt.dscalar.nii', metric_lbl)));
            %g = glob(fullfile(out_dir, '*', 'AnalyticSignalPar_call_subject_as_proc-norm_xyzt-ROI_*.mat'));
            leng = length(g);
            metric_mu = zeros(91282, 1);
            metric_sig2 = zeros(91282, 1);

            for gidx = 1:leng
                try
                    % tic
                    the_cifti = cifti_read(g{gidx});
                    metric_mu = metric_mu + the_cifti.cdata/leng; % mean
                    % toc % 0.8 - 1.6 sec
                catch ME
                    handwarning(ME)
                end
            end
            save(metric_lbl+"_mu.mat", "metric_mu");
            the_cifti.cdata = metric_mu;
            cifti_write(the_cifti, convertStringsToChars(metric_lbl+"_mu.dscalar.nii"));
        
            for gidx = 1:leng
                try
                    the_cifti = cifti_read(g{gidx});
                    metric_sig2 = metric_sig2 + (the_cifti.cdata - metric_mu).^2/leng; % var
                catch ME
                    handwarning(ME)
                end
            end
            metric_sig = metric_sig2.^(0.5);
            save(metric_lbl+"_sig.mat", "metric_sig");
            the_cifti.cdata = metric_sig;
            cifti_write(the_cifti, convertStringsToChars(metric_lbl+"_sig.dscalar.nii"));

            % coeff. of var.
            cov = the_cifti;
            cov.cdata  = metric_sig ./ metric_mu;
            cifti_write(cov, convertStringsToChars(metric_lbl+"_cov.dscalar.nii"));

            % snr
            snr = the_cifti;
            snr.cdata  = metric_mu ./ metric_sig;
            cifti_write(snr, convertStringsToChars(metric_lbl+"_snr.dscalar.nii"));
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
