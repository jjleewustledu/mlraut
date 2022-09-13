classdef ECoG
    %% line1
    %  line2
    %  
    %  Created 22-Apr-2022 12:53:31 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 9.12.0.1884302 (R2022a) for MACI64.  Copyright 2022 John J. Lee.
    
    methods (Static)
        function carls_coords2atl(varargin)
            %  Args:
            %      coords_file (text): contains coords (mm, from 4dfp center) on the pial surface of a 711-2B atlas
            %                          obtained by projection of electrodes visible on CT.  
            %  See also:
            %      Hacker/ECoG/atlas/mk_brainsurf_mesh_pts_test.csh
            %      imgsurf_4dfp, burn_sphere_4dfp            %      
            %      Hacker, et al. 2017.  https://doi.org/10.1016/j.neuroimage.2017.01.054

            ip = inputParser;
            addRequired(ip, 'coords_file', @isfile)
            addParameter(ip, 'atlas', '111', @istext)
            addParameter(ip, 'radius', 1.65, @istext) % 2.3 mm exposure of platinum electrode
            addParameter(ip, 'out_file', 'carls_coords2atl', @istext)
            parse(ip, varargin{:});
            ipr = ip.Results;
            if contains(ipr.atlas, '111')
                atl = fullfile(getenv('REFDIR'), '711-2B_111.4dfp.hdr');
            end
            if contains(ipr.atlas, '222')
                atl = fullfile(getenv('REFDIR'), '711-2B_222.4dfp.hdr');
            end
            if contains(ipr.atlas, '333')
                atl = fullfile(getenv('REFDIR'), '711-2B_333.4dfp.hdr');
            end
            
            cmd = sprintf('imgsurf_4dfp %s %s', atl, ipr.coords_file);
            mlbash(cmd);
            assert(isfile(strcat(ipr.coords_file, '.surf')))

            tbl = readtable(strcat(ipr.coords_file, '.surf'), 'FileType', 'text');
            tbl.Properties.VariableNames = {'x', 'y', 'z'};
            temp = tempname(pwd);
            tidx = 1;
            x = tbl.x(tidx);
            y = tbl.y(tidx);
            z = tbl.z(tidx);
            cmd = sprintf('burn_sphere_4dfp %g %g %g %s %s -o%g -v%i', ...
                x, y, z, atl, temp, ipr.radius, tidx);
            mlbash(cmd);
            for tidx = 2:size(tbl, 1)
                x = tbl.x(tidx);
                y = tbl.y(tidx);
                z = tbl.z(tidx);
                cmd = sprintf('burn_sphere_4dfp %g %g %g %s %s -a -o%g -v%i', ...
                    x, y, z, temp, ipr.out_file, ipr.radius, tidx);
                mlbash(cmd);
                mlfourdfp.FourdfpVisitor.copyfile_4dfp(myfileprefix(ipr.out_file), temp);
            end
            deleteExisting(temp);
        end
    end

    methods
        function this = ECoG(varargin)
            %% ECOG 
            %  Args:
            %      arg1 (its_class): Description of arg1.
            
            ip = inputParser;
            addParameter(ip, "arg1", [], @(x) false)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
