classdef Twistors < handle & mlsystem.IHandle
    %% line1
    %  line2
    %  
    %  Created 11-Dec-2024 22:23:24 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 24.2.0.2773142 (R2024b) Update 2 for MACA64.  Copyright 2024 John J. Lee.
    
    properties (Constant)
        ANT_COMMISSURE_COORD = [90, 128, 68]  % Conte69_AverageT1w.nii.gz
        PRECUNEUS_COORD = [90, 73, 116]
        LOCUS_CERULEUS_COORD = [90, 94, 60]
        V1_L_COORD = [81, 30, 69]
        V1_R_COORD = [99, 30, 69]
        DACC_COORD = [90, 158, 95]
        SCAN1_L_COORD = [60, 115, 137]  % Penfield shoulder
        SCAN2_L_COORD = [43, 117, 120]  % Penfield neck
        SCAN3_L_COORD = [33, 122, 85]  % Penfield tongue
        SCAN1_R_COORD = [120, 115, 137]
        SCAN2_R_COORD = [137, 117, 120]
        SCAN3_R_COORD = [147, 122, 85]

        theta_berry = 4*pi  % 2*pi
    end

    properties
        mask_ctx
        mask_cbm
        mask_HC
        mask_other
        mask_str
        mask_thal
    end

    properties (Dependent)
        bold_signal  % hilbert(bold) ~ Nt x Nx
        hp_thresh  % low freq. bound, Ryan ~ 0.01 Hz
        lp_thresh  % high freq. bound, Ryan ~ 0.05-0.1 Hz, for CSF studies ~ 0.1 Hz
        num_bins_angles  % this.ihcp_.num_bins_angles, needed by bin_by_physio_angle
        physio_signal  % hilbert(physio) ~ Nt x 1
        tr
        use_neg_ddt
        v_physio_is_inf % v_physio reaches head diameter in time << tr
        waves_dir
    end

    methods  %% GET
        function g = get.bold_signal(this)
            g = this.ihcp_.bold_signal;
        end
        function g = get.hp_thresh(this)
            g = this.ihcp_.hp_thresh;
        end
        function g = get.lp_thresh(this)
            g = this.ihcp_.lp_thresh;
        end
        function g = get.num_bins_angles(this)
            g = this.ihcp_.num_bins_angles;
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
        function g = get.waves_dir(this)
            g = this.ihcp_.waves_dir;
        end
    end

    methods
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

        function pos = greyordinate_concatenate(this, poset)
            %% [left, right, subcortex] ~ [Nx29696, Nx29716, Nx31870] ~ [Nx91282]; split subcortext x5

            assert(size(poset{1}, 2) ==  29696)
            assert(size(poset{2}, 2) == 29716)
            assert(size(poset{3}, 2) == 17853)
            assert(size(poset{4}, 2) == 1559)
            assert(size(poset{5}, 2) == 3553)
            assert(size(poset{6}, 2) == 2536)
            assert(size(poset{7}, 2) == 6369)
            pos(:,1:29696) = poset{1};
            pos(:,29697:59412) = poset{2};
            pos(:,this.mask_cbm) = poset{3};
            pos(:,this.mask_HC) = poset{4};
            pos(:,this.mask_str) = poset{5};
            pos(:,this.mask_thal) = poset{6};
            pos(:,this.mask_other) = poset{7};
        end

        function poset = greyordinate_split(this, pos)
            %% [left, right, subcortex] ~ [Nx29696, Nx29716, Nx31870] ~ [Nx91282]; split subcortext x5

            assert(size(pos, 2) == 91282)
            poset{1} = pos(:,1:29696);
            poset{2} = pos(:,29697:59412);
            poset{3} = pos(:,this.mask_cbm);
            poset{4} = pos(:,this.mask_HC);
            poset{5} = pos(:,this.mask_str);
            poset{6} = pos(:,this.mask_thal);
            poset{7} = pos(:,this.mask_other);
        end

        function pos = grayordinate_positions(this, precision, opts)
            %% pos ~ 3 x N_grayords, in mm

            arguments
                this mlraut.Twistors
                precision {mustBeTextScalar} = ""
                opts.center_coord {mustBeNumeric} = this.ANT_COMMISSURE_COORD  % Conte69_AverageT1w.nii.gz
                opts.type {mustBeTextScalar} = "midthickness"
                opts.qc logical = false
            end

            % cifti has models
            task = cifti_read(convertStringsToChars(this.ihcp_.task_dtseries_fqfn));
            models = task.diminfo{1}.models;

            % left cortex
            gl = gifti(convertStringsToChars(this.ihcp_.cohort_data.surf_gii_fqfn("L", type=opts.type)));
            gl_pos = gl.vertices';
            gl_pos = gl_pos(:, models{1}.vertlist+1);  % 0-start to 1-start; remove bad vertices
            % e.g., task.diminfo{1}.models{1} ~ 'CORTEX_LEFT', count = 29696, numvert = 32492

            % right cortex
            gr = gifti(convertStringsToChars(this.ihcp_.cohort_data.surf_gii_fqfn("R", type=opts.type)));
            gr_pos = gr.vertices';
            gr_pos = gr_pos(:, models{2}.vertlist+1);  % 0-start to 1-start; remove bad vertices
            % e.g., task.diminfo{1}.models{2} ~ 'CORTEX_RIGHT', count = 29716, numvert = 32492

            % subcortical
            outinfo = cifti_diminfo_dense_get_volume_all_info(task.diminfo{1});
            ic = mlfourd.ImagingContext2(this.ihcp_.wmparc_fqfn);
            sc_pos = this.voxel_indices_to_position(outinfo.voxlist1, ic);

            pos = [gl_pos, gr_pos, sc_pos];  % [3x29696, 3x29716, 3x31870] ~ [3x91282]
            pos = pos + ascol(this.ANT_COMMISSURE_COORD) - ascol(opts.center_coord);

            if strcmpi(precision, "single") && ~isa(pos, "single")
                pos = single(pos);
            end
            if strcmpi(precision, "double") && ~isa(pos, "double")
                pos = double(pos);
            end

            % qc
            if opts.qc
                this.visualize_density(pos);

                leftinfo = cifti_diminfo_dense_get_surface_info(task.diminfo{1}, 'CORTEX_LEFT');
                rightinfo = cifti_diminfo_dense_get_surface_info(task.diminfo{1}, 'CORTEX_RIGHT');
                info.cortex_left_idx = leftinfo.ciftilist;
                info.cortex_right_idx = rightinfo.ciftilist;
                info.subcortical_idx = outinfo.ciftilist;
                this.visualize_greyordinate_subsets(pos, info);
            end
        end

        function this = load_anatomy_masks(this)
            ld = load(fullfile(this.waves_dir, 'supporting_files', 'mask_ctx_HCP.mat'));
            this.mask_ctx = ld.mask_ctx;
            ld = load(fullfile(this.waves_dir, 'supporting_files', 'mask_cbm_HCP.mat'));
            this.mask_cbm = ld.mask_cbm;
            ld = load(fullfile(this.waves_dir, 'supporting_files', 'mask_HC_HCP.mat'));
            this.mask_HC = ld.mask_HC;
            ld = load(fullfile(this.waves_dir, 'supporting_files', 'mask_str_HCP.mat'));
            this.mask_str = ld.mask_str;
            ld = load(fullfile(this.waves_dir, 'supporting_files', 'mask_thal_HCP.mat'));
            this.mask_thal = ld.mask_thal;

            this.mask_other = ~this.mask_ctx & ~this.mask_cbm & ~this.mask_HC & ~this.mask_str & ~this.mask_thal;
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
            bold_pos = this.grayordinate_positions(type="sphere");  % 3 x N_x
            x = abs(vecnorm(bold_pos - opts.physio_pos));  % Nx x 1, in mm

            Nt = size(propagated_signal, 1);
            tr_ = this.ihcp_.tr;
            times = ascol(0:tr_:(Nt - 1)*tr_);
            times = times + tr_/2;  % mid-frame times
            for tidx_phys = 1:Nt
                t = (tidx_phys - 0.5)*tr_;  % mid-frame time
                tplus = t + tr_/2;
                tminus = t - tr_/2;
                selection = tminus*v <= x & x <= tplus*v;  % frame boundaries overlap
                propagated_signal(:, selection) = ...
                    this.forward(physio_signal, times, dt=x(selection)/v, unvisited_is_inf=opts.unvisited_is_inf);
            end
        end

        %% Twistor elements & related

        function binned = bin_by_physio_angle(this, psi, phi, opts)
            %% Average psi(t) into bins of angles theta <- angle(phi), whereby {theta_i, i \in \mathbb{N}} 
            %  span [-2 pi, 2 pi], allowing assessment of Berry's phase.  This allows generation of figures 
            %  similar to Ryans' Science Adv. (2021) Fig. 4.
            %
            %  Args:
            %      this mlraut.Cifti
            %      psi {mustBeNumeric}  % already gsr, centered, filtered, rescaled, analytic; Nt x Ngo
            %      phi {mustBeNumeric}  % already centered, filtered, rescaled, analytic; Nt x 1
            %      opts.smoothing {mustBeInteger} = 3  % set to [] to not smooth
            %      opts.num_bins_angles {mustBeScalarOrEmpty} = 40
            %  Returns:
            %      binned numeric ~ num_bins_angles x Ngo; analytic

            arguments
                this mlraut.Twistors
                psi {mustBeNumeric}  % already gsr, centered, filtered, rescaled, analytic; Nt x Ngo
                phi {mustBeNumeric}  % already centered, filtered, rescaled, analytic; Nt x 1
                opts.smoothing {mustBeInteger} = []  % unused, but option is in the API
                opts.num_bins_angles {mustBeScalarOrEmpty} = this.num_bins_angles
            end
            phi = phi(1:size(psi, 1), :);
            num_bins_angles_ = opts.num_bins_angles;
            binlim = asrow(linspace(-this.theta_berry/2, this.theta_berry/2, num_bins_angles_ + 1));

            % init
            binned = zeros(num_bins_angles_, size(psi, 2));

            % wrapped physio is not unwrapped
            if size(phi, 2) > 1
                phi = mean(phi, 2);
            end
            theta = unwrap(angle(phi));
            Nberry = this.theta_berry/(2*pi);
            wrapped_theta = Nberry*angle(exp(1i*theta/Nberry));  % in [-2 pi, 2 pi] for this.theta_berry = 4 pi

            % average bold by phase bins
            for b = 2:num_bins_angles_+1
                selected = binlim(b-1) < wrapped_theta & wrapped_theta < binlim(b);
                binned(b-1,:) = mean(psi(selected, :), 1, "omitnan");
            end
        end

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
            psi = -1i*(psi.*conj(phi) - phi.*conj(psi))/sqrt(2);
        end

        function psi = Z(this, psi, phi)
            %% of twistor

            if this.use_neg_ddt
                [psi,phi] = this.neg_ddt(psi, phi);
                psi = this.ihcp_.build_centered_and_rescaled(psi);
            end
            psi = (psi.*conj(psi) - phi.*conj(phi))/sqrt(2);
        end

        function [Z_sup_alpha_,x_sup_AAp] = Z_sup_alpha(this, pi_sub_0p, pi_sub_1p, opts)
            %% S & ST, Penrose & Rindler, sec. 6.2.
            %  Args:
            %      pi_sub_0p {mustBeNumeric} \pi_{0'} ~ psi ~ -dbold/dt
            %      pi_sub_1p {mustBeNumeric} \pi_{1'} ~ phi ~ physio(t)
            %      opts.type {mustBeTextScalar} = "midthickness"  % "sphere", "midthickness"
            %      opts.center_coord {mustBeNumeric} = [90, 73, 116]  % mm, for precuneus
            %      opts.go_qc logical = false  % greyordinate qc
            %      opts.c {mustBeScalarOrEmpty} = this.lp_thresh  % speed limit for Poincare invariance ~ 1/timescale
            %      opts.use_E4 logical = true  % Ward & Wells sec. 8.1
            %      opts.s {mustBeScalarOrEmpty} = 0  % helicity in Penrose & Rindler eq. 6.2.7
            %  Returns:
            %      Z_sup_alpha_ (4-cell of single|double \mathsf{Z}^{\alpha}) ~ 
            %                   N_t x N_go for {\mathsf{Z}^0, \mathsf{Z}^1, \mathsf{Z}^2, \mathsf{Z}^3},
            %                   but since \mathsf{Z}^3 is the spatially uniform arousal, exclude to save storage.
            %      x_sup_AAp (2x2-cell of single|double x^{AA'}

            arguments
                this mlraut.Twistors
                pi_sub_0p {mustBeNumeric}
                pi_sub_1p {mustBeNumeric}
                opts.type {mustBeTextScalar} = "sphere"  % "sphere", "midthickness"
                opts.center_coord {mustBeNumeric} = this.LOCUS_CERULEUS_COORD  % mm, for precuneus
                opts.go_qc logical = false  % greyordinate qc
                opts.c {mustBeScalarOrEmpty} = this.lp_thresh  % speed limit for Poincare invariance ~ 0.1 ~ 1/timescale
                opts.use_E4 logical = false  % Ward & Wells sec. 8.1
                opts.s {mustBeScalarOrEmpty} = -1/2^(5/2)  % helicity in Penrose & Rindler eq. 6.2.7; only influences ~opts.use_E4
            end

            % consider using -dbold/dt
            if this.use_neg_ddt
                [pi_sub_0p,pi_sub_1p] = this.neg_ddt(pi_sub_0p, pi_sub_1p);
                pi_sub_0p = this.ihcp_.build_centered_and_rescaled(pi_sub_0p);
            end

            % init
            N_t = size(pi_sub_0p, 1);
            N_go = size(pi_sub_0p, 2);
            if ~all(size(pi_sub_0p) == size(pi_sub_1p))  % ensure matched sizes for psi & phi
                assert(all(size(pi_sub_1p) == [N_t, 1]))
                pi_sub_1p = repmat(pi_sub_1p, [1, N_go]);
            end

            % Euclidean coords
            go = this.grayordinate_positions(type=opts.type, qc=opts.go_qc, center_coord=opts.center_coord);  % 3 x N_go (mm)
            assert(all(size(go) == [3, N_go]))            
            func = str2func(this.ihcp_.rescaling);  % e.g., iqr
            try
                L = func(go, "all");
            catch ME
                handwarning(ME)
                L = func(go, [], "all");  % characteristic length scale of hemisphere (mm)
            end
            c = opts.c;
            x = go(1, :) / L;
            y = go(2, :) / L;
            z = go(3, :) / L;
            t = c * ascol(linspace(0, (N_t - 1)*this.tr, N_t));
            bulk = zeros(N_t, 1);  % forces shape with N_t rows

            % Use \mathbb{E}^4 according to Ward & Wells, pg. 388
            if opts.use_E4
                x_sup_00p = (-1i * t + z) / sqrt(2);
                x_sup_01p = (x - 1i * y + bulk) / sqrt(2);
                x_sup_10p = (x + 1i * y + bulk) / sqrt(2);
                x_sup_11p = (-1i * t - z) / sqrt(2);
                omega_sup_0 = 1i * ( ...
                    x_sup_00p .* pi_sub_0p + x_sup_01p .* pi_sub_1p);
                omega_sup_1 = 1i * ( ...
                    x_sup_10p .* pi_sub_0p + x_sup_11p .* pi_sub_1p);
            else
                % Penrose & Rindler equations 6.1.10, 6.2.7
                omega_o_sup_0 = 0;  % Twistor angular momentum ~ [\mathring{omega}^0; \mathring{omega}^1]
                omega_o_sup_1 = opts.s;  % \mathring{\omega}^1 ~ helicity; trying anti-self dual setting
                x_sup_00p = (t + z) / sqrt(2);
                x_sup_01p = (bulk + x + 1i * y) / sqrt(2);
                x_sup_10p = (bulk + x - 1i * y) / sqrt(2);
                x_sup_11p = (t - z) / sqrt(2);
                omega_sup_0 = omega_o_sup_0 - 1i * ( ...
                    x_sup_00p .* pi_sub_0p + x_sup_01p .* pi_sub_1p);
                omega_sup_1 = omega_o_sup_1 - 1i * ( ...
                    x_sup_10p .* pi_sub_0p + x_sup_11p .* pi_sub_1p);
            end

            x_sup_AAp = {x_sup_00p, x_sup_01p; x_sup_10p, x_sup_11p};
            Z_sup_alpha_ = {omega_sup_0, omega_sup_1, pi_sub_0p, pi_sub_1p};
        end

        function psi = T(this, psi, phi)
            %% of twistor

            if this.use_neg_ddt
                [psi,phi] = this.neg_ddt(psi, phi);
                psi = this.ihcp_.build_centered_and_rescaled(psi);
            end
            psi = (psi.*conj(psi) + phi.*conj(phi))/sqrt(2);
        end

        function x_ = x(this, psi, phi)
            x_ = this.X(psi, phi) ./ this.T(psi, phi);
        end

        function y_ = y(this, psi, phi)
            y_ = this.Y(psi, phi) ./ this.T(psi, phi);
        end

        function z_ = z(this, psi, phi)
            z_ = this.Z(psi, phi) ./ this.T(psi, phi);
        end

        function z = zeta(this, xi, eta)
            %% of twistor

            arguments
                this mlraut.Twistors
                xi {mustBeNumeric}
                eta {mustBeNumeric}
            end
            assert(~isempty(eta))
            assert(size(xi, 1) == size(eta, 1))

            if this.use_neg_ddt
                [xi,eta] = this.neg_ddt(xi, eta);
                xi = this.ihcp_.build_centered_and_rescaled(xi);
            end

            z = xi ./ eta;
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

            arguments
                this mlraut.Twistors
                psi {mustBeNumeric}
                phi {mustBeNumeric} = []
            end

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
            % this.load_anatomy_masks()
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
            %
            %  Returns:
            %      psi <- -dpsi/dt
            %      phi <- phi(1:end-1) to match new size of psi

            arguments
                psi {mustBeNumeric}
                phi {mustBeNumeric} = []
                opts.scale {mustBeScalarOrEmpty} = 1
            end

            psi = -opts.scale*diff(psi, 1);  % 1st deriv. of time
            % psi(psi < 0) = 0;  % see Fultz et al. 2019

            if ~isempty(phi)
                Nt = size(phi, 1);
                rep = repmat({':'}, 1, ndims(phi) - 1);
                phi = phi(1:Nt-1, rep{:});
            end
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
        
        function visualize_density(GC)
            positions = GC';

            % Calculate local density using k-nearest neighbors
            k = 50;  % Number of neighbors to consider
            [idx, distances] = knnsearch(positions, positions, 'K', k+1);

            % Use mean distance to k nearest neighbors as inverse density measure
            mean_distances = mean(distances(:, 2:end), 2);
            density_values = 1 ./ mean_distances;
            density_normalized = (density_values - min(density_values)) / ...
                (max(density_values) - min(density_values));

            % Create color map based on density
            colors = hot(256);  % Use hot colormap for density
            color_indices = round(density_normalized * 255) + 1;
            point_colors = colors(color_indices, :);

            % Create and display density-colored point cloud
            ptCloud = pointCloud(positions, 'Color', uint8(point_colors * 255));

            figure('Position', [100, 100, 1000, 800]);
            pcshow(ptCloud, 'MarkerSize', 5);
            title('Greyordinate Density Visualization');
            colormap(hot);

            % Correct way to set colorbar label
            cb = colorbar;
            cb.Label.String = 'Relative Density';
            cb.Label.FontSize = 12;

            axis equal;
            view(45, 20);
        end
        
        function visualize_greyordinate_subsets(GC, greyord_info)
            % Create subplots for different anatomical regions
            figure('Position', [50, 50, 1600, 500]);

            % Left cortex only
            subplot(1, 3, 1);
            left_coords = GC(:, greyord_info.cortex_left_idx)';
            ptCloud_left = pointCloud(left_coords);
            pcshow(ptCloud_left, 'MarkerSize', 4);
            title('Left Cortex');
            view(180, 0);  % Lateral view
            axis equal; grid on;

            % Right cortex only
            subplot(1, 3, 2);
            right_coords = GC(:, greyord_info.cortex_right_idx)';
            ptCloud_right = pointCloud(right_coords);
            pcshow(ptCloud_right, 'MarkerSize', 4);
            title('Right Cortex');
            view(0, 0);  % Lateral view
            axis equal; grid on;

            % Subcortical structures only
            subplot(1, 3, 3);
            subcort_coords = GC(:, greyord_info.subcortical_idx)';
            ptCloud_subcort = pointCloud(subcort_coords);
            pcshow(ptCloud_subcort, 'MarkerSize', 8);
            title('Subcortical Structures');
            view(45, 20);
            axis equal; grid on;

            sgtitle('Greyordinate Components', 'FontSize', 16);
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
