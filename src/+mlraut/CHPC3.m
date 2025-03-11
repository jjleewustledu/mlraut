classdef CHPC3
    %% Provides support functions for using the Matlab Parallel Server at CHPC3.
    %  
    %  Created 07-Apr-2022 16:13:53 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 9.10.0.1851785 (R2021a) Update 6 for MACI64.  Copyright 2022 John J. Lee.
    
    methods (Static)
        function clean_tempdir()
            try
                deleteExisting(fullfile(tempdir, '*.nii*'));
                deleteExisting(fullfile(tempdir, '*.save'));
            catch ME
                disp(ME)
            end
        end
        function [t, A] = parallel_example(iter)
            if nargin==0
                iter = 8;
            end
            disp('Start sim')

            t0 = tic;
            parfor idx = 1:iter
                A(idx) = idx;
                pause(2)
                idx
            end
            t = toc(t0);

            disp('Sim completed')
            save RESULTS A
        end

        function c = propcluster(account_name, opts)
            arguments
                account_name = 'aristeidis_sotiras'  % 'joshua_shimony' 'manu_goyal' 'jjlee'
                opts.partition = 'tier2_cpu'
            end
            if ~strcmp(account_name, 'aristeidis_sotiras')
                opts.partition = 'tier1_cpu';
            end
            account_name = convertStringsToChars(account_name);
            opts.partition = convertStringsToChars(opts.partition);

            c = parcluster;
            c.AdditionalProperties.AccountName = account_name;
            c.AdditionalProperties.AdditionalSubmitArgs = sprintf('--account=%s', account_name);
            c.AdditionalProperties.ClusterHost = 'login3.chpc.wustl.edu';
            c.AdditionalProperties.EmailAddress = '';
            c.AdditionalProperties.EnableDebug = false;
            c.AdditionalProperties.GpusPerNode = 0;
            c.AdditionalProperties.MemPerCPU = '100gb';
            % c.AdditionalProperties.Node = '';
            c.AdditionalProperties.Partition = opts.partition;
            c.AdditionalProperties.RemoteJobStorageLocation = '/home/jjlee/.matlab/3p_cluster_jobs/chpc/twistor.attlocal.net.dhcp.wustl.edu/R2024b/nonshared';
            c.AdditionalProperties.UseIdentityFile = false;
            c.AdditionalProperties.UseSmpd = false;
            c.AdditionalProperties.Username = 'jjlee';
            c.AdditionalProperties.WallTime = '04:00:00';
            disp(c.AdditionalProperties)
        end

        function c = propcluster_tiny(account_name, opts)
            arguments
                account_name = 'aristeidis_sotiras'  % 'joshua_shimony' 'manu_goyal' 'jjlee'
                opts.partition = 'tier2_cpu'
            end
            if ~strcmp(account_name, 'aristeidis_sotiras')
                opts.partition = 'tier1_cpu';
            end
            account_name = convertStringsToChars(account_name);
            opts.partition = convertStringsToChars(opts.partition);

            c = parcluster;
            c.AdditionalProperties.AccountName = account_name;
            c.AdditionalProperties.AdditionalSubmitArgs = sprintf('--account=%s', account_name);
            c.AdditionalProperties.ClusterHost = 'login3.chpc.wustl.edu';
            c.AdditionalProperties.EmailAddress = '';
            c.AdditionalProperties.EnableDebug = true;
            c.AdditionalProperties.GpusPerNode = 0;
            c.AdditionalProperties.MemPerCPU = '16gb';
            % c.AdditionalProperties.Node = '';
            c.AdditionalProperties.Partition = opts.partition;
            c.AdditionalProperties.RemoteJobStorageLocation = '/home/jjlee/.matlab/3p_cluster_jobs/chpc/twistor.attlocal.net.dhcp.wustl.edu/R2024b/nonshared';
            c.AdditionalProperties.UseIdentityFile = false;
            c.AdditionalProperties.UseSmpd = false;
            c.AdditionalProperties.Username = 'jjlee';
            c.AdditionalProperties.WallTime = '100:00:00';
            disp(c.AdditionalProperties)
        end

        function setenvs()
            [~,r] = system('hostname');
            % if ~contains(r, 'cluster')
            %     return
            % end

            setenv('TMPDIR', '/scratch/jjlee/tmp') % worker nodesk

            setenv('SINGULARITY_HOME', '/scratch/jjlee/Singularity')
            setenv('AFNIPATH', '/export/afni/afni-20.3.03/linux_openmp_64')
            setenv('ANTSPATH', '/export/ants/ants-2.3.5/bin')
            setenv('DEBUG', '');
            setenv('FREESURFER_HOME', '/export/freesurfer/freesurfer-7.4.1')
            setenv('FSLDIR', '/export/fsl/fsl-6.0.7.8')

            setenv('FSLOUTPUTTYPE', 'NIFTI_GZ')
            setenv('FSLMULTIFILEQUIT', 'TRUE')
            setenv('FSLTCLSH', fullfile(getenv('FSLDIR'),'bin','fsltclsh'))
            setenv('FSLWISH', fullfile(getenv('FSLDIR'),'bin','fslwish'))
            setenv('FSLLOCKDIR', '')
            setenv('FSLMACHINELIST', '')
            setenv('FSLREMOTECALL', '')
            setenv('FSLREMOTECALL', 'cuda.q')
            setenv('PYOPENGL_PLATFORM', 'osmesa')

            setenv('REFDIR', '/home/jjlee/.local/atlas')
            setenv('RELEASE', '/home/jjlee/.local/lin64-tools')            
            setenv('PATH', ...
                strcat(getenv('RELEASE'), ':', ...
                       getenv('AFNIPATH'), ':', ...
                       fullfile(getenv('FREESURFER_HOME'), 'bin'), ':', ...
                       fullfile(getenv('FSLDIR'), 'bin'), ':', ...
                       '/export/singularity/singularity-3.9.0/bin', ':', ...
                       '/usr/bin', ':', ...
                       getenv('PATH')))
            setenv('LD_LIBRARY_PATH', ...
                strcat('/usr/lib64', ':', getenv('LD_LIBRARY_PATH'))) % need libOSMesa.so.8 for fsleyes render
                   
            %disp("mlraut.CHPC3.setenvs():getenv('PATH'):")
            %disp(getenv('PATH'))
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
