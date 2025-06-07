classdef Twistors < handle & mlsystem.IHandle
    %% line1
    %  line2
    %  
    %  Created 11-Dec-2024 22:23:24 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 24.2.0.2773142 (R2024b) Update 2 for MACA64.  Copyright 2024 John J. Lee.
    
    properties        
    end

    properties (Dependent)
        bold_signal  % hilbert(bold) ~ Nt x Nx
        physio_signal  % hilbert(physio) ~ Nt x 1
        tr
        use_neg_ddt
        v_physio_is_inf % v_physio reaches head diameter in time << tr
    end

    methods  %% GET
        function g = get.bold_signal(this)
            g = this.ihcp_.bold_signal;
        end
        function g = get.physio_signal(this)
            g = this.ihcp_.physio_signal;
        end
        function g = get.tr(this)
            g = this.ihcp_.tr;
        end
        function g = get.use_neg_ddt(this)
            g = this.ihcp_.use_neg_ddt;
        end
        function g = get.v_physio_is_inf(this)
            g = this.ihcp_.v_physio_is_inf;
        end
    end

    methods
        function binned = bin_by_physio_angle(this, varargin)
            binned = this.ihcp_.cifti.bin_by_physio_angle(varargin{:});
        end

        function [pos,v] = center_of_mass_position(this, roi)
            %% pos ~ 3 x 1, in mm
            %  v ~ 3 x 1, indices of imaging img

            arguments
                this mlraut.Twistors
                roi mlfourd.ImagingContext2 = this.ihcp_.roi
            end

            if isempty(roi) || strcmp(roi.stateTypeclass, "mlfourd.TrivialTool")
                pos = [0; 0; 0];
                v = [1; 1; 1];
                return
            end

            % Gemini Advanced 1.5 Pro, 20241211 230000
            % I have a 3D array in matlab representing a grid of positions in space. The array contains numbers
            % describing a cluster of massive objects. How can I quickly calculate position of the center of mass of the
            % cluster?

            mass_array = double(roi);

            % Define grid dimensions
            x_dim = size(mass_array, 1);
            y_dim = size(mass_array, 2);
            z_dim = size(mass_array, 3);

            % Create coordinate grids
            [y_coords, x_coords, z_coords] = meshgrid(1:y_dim, 1:x_dim, 1:z_dim);

            % Total mass of the cluster
            total_mass = sum(mass_array(:)); 

            % Calculate the weighted coordinates, center of mass
            x_cm = sum(x_coords(:) .* mass_array(:)) / total_mass;
            y_cm = sum(y_coords(:) .* mass_array(:)) / total_mass;
            z_cm = sum(z_coords(:) .* mass_array(:)) / total_mass;

            v = [x_cm; y_cm; z_cm];
            pos = this.voxel_indices_to_position(v, roi);
        end

        function psi = connectivity(this, varargin)
            psi = this.ihcp_.connectivity(varargin{:});
        end

        function pos = grayordinate_positions(this)
            %% pos ~ 3 x N_grayords, in mm

            % cifti has models
            task = cifti_read(convertStringsToChars(this.ihcp_.task_dtseries_fqfn));
            models = task.diminfo{1}.models;

            % left cortex
            gl = gifti(convertStringsToChars(this.ihcp_.cohort_data.surf_gii_fqfn("L")));
            gl_pos = gl.vertices';
            gl_pos = gl_pos(:, models{1}.vertlist+1);  % 0-start to 1-start; remove bad vertices
            % e.g., task.diminfo{1}.models{1} ~ 'CORTEX_LEFT', count = 29696, numvert = 32492

            % right cortex
            gr = gifti(convertStringsToChars(this.ihcp_.cohort_data.surf_gii_fqfn("R")));
            gr_pos = gr.vertices';
            gr_pos = gr_pos(:, models{2}.vertlist+1);  % 0-start to 1-start; remove bad vertices
            % e.g., task.diminfo{1}.models{2} ~ 'CORTEX_RIGHT', count = 29716, numvert = 32492

            % subcortical
            outinfo = cifti_diminfo_dense_get_volume_all_info(task.diminfo{1});
            ic = mlfourd.ImagingContext2(this.ihcp_.wmparc_fqfn);
            sc_pos = this.voxel_indices_to_position(outinfo.voxlist1, ic);

            % qc:
            % extracted = zeros(outinfo.voldims, 'single');
            % extracted(cifti_vox2ind(outinfo.voldims, outinfo.voxlist1)) = sc.cdata(outinfo.ciftilist, 1);
            % ifc = ic.imagingFormat; ifc.img = extracted;
            % ifc.view

            pos = [gl_pos, gr_pos, sc_pos];  % [3x29696, 3x29716, 3x31870] ~ [3x91282]
        end

        function propagated_signal = propagate_physio(this, physio_signal, opts)
            %% Propagate physio_signal from physio_pos radially according to propagation velocity v.
            %  Prior to arrival of physio signal, the array for propagated_signal will be Inf or one
            %  according to unvisited_is_inf.
            %  bold_signal is only needed for array sizes.
            %
            %  Args:
            %      this mlraut.Twistors
            %      physio_signal double  % N_t x 1
            %      opts.size_bold_signal double  % [N_t, N_x]
            %      opts.physio_pos double = this.center_of_mass_position()  % 3 x 1
            %      opts.unvisited_is_inf logical = false            
            %  Returns:
            %      propagated_signal (N_t, N_x) double 

            arguments
                this mlraut.Twistors
                physio_signal double  % N_t x 1
                opts.size_bold_signal double  % [N_t, N_x]
                opts.physio_pos double = this.center_of_mass_position()  % 3 x 1
                opts.unvisited_is_inf logical = false
            end
            
            if this.v_physio_is_inf
                propagated_signal = ones(opts.size_bold_signal).*physio_signal;
                return
            end

            % privately, use mm/s
            v = 1e3*this.ihcp_.v_physio;

            if opts.unvisited_is_inf
                big = 1e3*max(physio_signal, [], "all");
                propagated_signal = big*ones(opts.size_bold_signal);
            else
                propagated_signal = physio_signal(1)*ones(opts.size_bold_signal);
            end
            bold_pos = this.grayordinate_positions();  % 3 x N_x
            x = abs(vecnorm(bold_pos - opts.physio_pos));  % Nx x 1, in mm

            Nt = size(propagated_signal, 1);
            tr = this.ihcp_.tr;
            times = ascol(0:tr:(Nt - 1)*tr);
            times = times + tr/2;  % mid-frame times
            for tidx_phys = 1:Nt
                t = (tidx_phys - 0.5)*tr;  % mid-frame time
                tplus = t + tr/2;
                tminus = t - tr/2;
                selection = tminus*v <= x & x <= tplus*v;  % frame boundaries overlap
                propagated_signal(:, selection) = ...
                    this.forward(physio_signal, times, dt=x(selection)/v, unvisited_is_inf=opts.unvisited_is_inf);
            end
        end

        %% Twistor elements & related

        function psi = X(this, psi, phi)
            %% of twistor

            if this.use_neg_ddt
                [psi,phi] = this.neg_ddt(psi, phi);
                psi = this.ihcp_.build_centered_and_rescaled(psi);
            end
            psi = (psi.*conj(phi) + phi.*conj(psi))/sqrt(2);
        end

        function psi = Y(this, psi, phi)
            %% of twistor

            if this.use_neg_ddt
                [psi,phi] = this.neg_ddt(psi, phi);
                psi = this.ihcp_.build_centered_and_rescaled(psi);
            end
            psi = (psi.*conj(phi) - phi.*conj(psi))/sqrt(2i);
        end

        function psi = Z(this, psi, phi)
            %% of twistor

            if this.use_neg_ddt
                [psi,phi] = this.neg_ddt(psi, phi);
                psi = this.ihcp_.build_centered_and_rescaled(psi);
            end
            psi = (psi.*conj(psi) - phi.*conj(phi))/sqrt(2);
        end

        function psi = T(this, psi, phi)
            %% of twistor

            if this.use_neg_ddt
                [psi,phi] = this.neg_ddt(psi, phi);
                psi = this.ihcp_.build_centered_and_rescaled(psi);
            end
            psi = (psi.*conj(psi) + phi.*conj(phi))/sqrt(2);
        end

        function theta = angle(this, psi, phi)
            if this.use_neg_ddt
                [psi,phi] = this.neg_ddt(psi, phi);
                psi = this.ihcp_.build_centered_and_rescaled(psi);
            end
            theta = angle(psi.*phi);            
        end

        function psi = neg_dbold_dt(this, psi, phi)
            %% of twistor

            psi = this.neg_ddt(psi, phi);
            psi = this.ihcp_.build_centered_and_rescaled(psi);
        end

        function plvs = phase_locked_values(~, psi, phi)
            %% equivalent to Ryan's expression
            %  plvs(:,s,t) = nanmean(exp(1i*(bsxfun(@minus,unwrap(angle(h1)),unwrap(angle(h2))))));
            %  in physio_phase_mapping.m
            
            as = psi .* conj(phi);
            plvs = conj(as) ./ abs(as);
        end

        function theta = unwrap(this, psi, phi)
            %% floating-point equivalent to:
            %  theta = unwrap(angle(psi)) - unwrap(angle(phi))

            theta = unwrap(this.angle(psi, phi));
        end

        %% ctor

        function this = Twistors(ihcp)
            arguments
                ihcp mlraut.HCP
            end

            this.ihcp_ = ihcp;
        end
    end

    methods (Static)
        function signal1 = forward(signal, times, opts)
            arguments
                signal double  % Nt x 1
                times double  % Nt x 1
                opts.dt double  % 1 x Nx
                opts.unvisited_is_inf logical = false
            end
            Nt = length(times);
            Nx = length(opts.dt);
            signal1 = ones(Nt, Nx);
            if opts.unvisited_is_inf 
                boundary = 1e3*max(signal, [], "all");
            else
                boundary = signal(1);
            end
            for idx = 1:Nx
                signal1(:,idx) = interp1(times, signal, times-opts.dt(idx), "linear", boundary);
            end
        end

        function [psi,phi] = neg_ddt(psi, phi, opts)
            %% scale should be the approximate period of an oscillatory signal; 
            %  see also Test_Twistors.test_neg_dtt().

            arguments
                psi {mustBeNumeric}
                phi {mustBeNumeric}
                opts.scale {mustBeScalarOrEmpty} = 1
            end

            Nt = size(phi, 1);
            rep = repmat({':'}, 1, ndims(phi) - 1);
            phi = phi(1:Nt-1, rep{:});

            psi = -opts.scale*diff(psi, 1);  % 1st deriv. of time
            psi(psi < 0) = 0;  % see Fultz et al. 2019
        end

        function psi = selective_X(mat, opts)
            %% of twistor

            arguments
                mat {mustBeFile}
                opts.selection {mustBeTextScalar} = "nonpositivity"
                opts.save_cifti logical = true
                opts.out_dir {mustBeTextScalar} = ""
            end

            ld = load(mat);
            psi = ld.this.bold_signal;
            phi = ld.this.physio_signal;

            psi = (psi.*conj(phi) + phi.*conj(psi))/sqrt(2);

            switch opts.selection
                case "negativity"
                    select = real(mean(psi, 2)) < 0;  % mean over all vertices
                case "positivity"
                    select = real(mean(psi, 2)) > 0;  % mean over all vertices
                case "nonnegativity"
                    select = real(mean(psi, 2)) >= 0;  % mean over all vertices
                case "nonpositivity"
                    select = real(mean(psi, 2)) <= 0;  % mean over all vertices
                otherwise
                    return
            end

            if opts.save_cifti
                if ~isemptytext(opts.out_dir)
                    ld.this.out_dir = opts.out_dir;
                end
                tags = ld.this.tags(opts.selection);
                this.ihcp_.cifti.write_ciftis( ...
                    psi, ...
                    sprintf('X_as_sub-%s_ses-%s_%s', ld.this.current_subject, ld.this.current_task, tags), ...
                    partitions=select, ...
                    do_save_dynamic=false);
            end
        end

        function pos = voxel_indices_to_position(v, ic)
            %% v ~ 3 x N voxel indices (starting at 1), N >= 1
            %  ic ~ mlfourd.ImagingContext2
            %  pos ~ 3 x N positions in scanner coords, centered at anterior commissure

            arguments
                v {mustBeNumeric}
                ic mlfourd.ImagingContext2
            end

            sz = size(v);
            if sz(1) == 3
                v1 = [v - 1; ones(1, sz(2))];
            end
            assert(size(v1, 1) == 4)

            roi = ic.imagingFormat;
            T = [roi.hdr.hist.srow_x; roi.hdr.hist.srow_y; roi.hdr.hist.srow_z; 0 0 0 1];
            pos = T*v1;
            pos = pos(1:3, :);
        end
    end

    %% PRIVATE

    properties (Access=private)
        ihcp_  % mlraut.HCP
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
