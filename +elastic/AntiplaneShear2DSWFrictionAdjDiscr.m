classdef AntiplaneShear2DSWFrictionAdjDiscr < noname.Discretization
% Discretizes the 1D anti-plane shear problem
% 
% u_tt = 1/rho ( (mu*u_x)_x + S)
% Psi_t = g(V, Psi)
% 
% for velocity potential u, and state variable Psi, where
% S is a dirac delta source term, and g a state evolution law.
properties
    name   = "AntiplaneShear2DSWFrictionAdjDiscr"
    description  = 'Discretizes the 1D wave equation for antiplane shear with a slip weakening friction interface. Solves for displacement u, velocity v and slip s'
    order                % Order of accuracy
    grid                 % Grid
    material             % Material parameters
    
    % Discretization matrix. Formulate 1st-order in time system
    % v = u_t
    % w = [u; v; Psi]
    % such that
    % D = [0      I     0;
    %      1/rho D2(mu) 0;
    %      0      0     0];  
    % F = [0; 1/rho(S + SAT); g]^T
    % w_t = D*w + F(w)
    D

    % Quadratures & norms
    H  % Spatial quadature of the SBP operator

    % DiffOp
    mbDiffOp

    % Component operators for first-order in time system
    % E.u*w = u, Eu'*u = [u; 0; 0]; 
    % E.v*w = v, Ev'*v = [0; v; 0];
    % E.Psi*w = Psi, EPsi'*Psi = [0; 0; Psi];
    E

    % Initial conditions
    w0 % w0 = [u0; v0; Psi0]
        
    % Friction and operators
    friction             % Struct holding parameters and friction functions/data
    e_fault_m            % Fault restriciton operator on minus side
    e_fault_p            % Fault restriciton operator on plus side
    penalty_fault        % Fault data penalty operator
    fault_traction       % Fault traction forcing
    %%% modified by R.F.
    % state_evo            % State evolution forcing
    % ↓
    strength_evo            % Strength evolution forcing
    %%%
    V_star
    tau_l_fun
    eta

    % Sources
    sources         % Struct holding parameters source functions/data
    S               % Continuous source term
    dirac_deltas    % dirac delta functions
    %%% added by R.F.
    bidm              
    H_flt
    flt_pos

    fault_sources
    adj_traction
    fault_dirac_deltas
    %%%
end

methods
    
    %   m          --      Number of grid points for each domain
    %   lims       --      Limits [x_l, x_i, x_r], specifying left/right boundaries and interface positions.
    %   order      --      order of accuracy
    %   material   --      struct for material properties
    %       material.rho   cell array with function rho(x) for density in each domain                 
    %       material.mu    cell array with function mu(x) for shear modulus in each domain                 
    %   bc         --      struct for boundary conditions
    %       bc.type             --  cell array with boundary conditions for left/right boundary
    %                               'dirichlet': Dirichlet conditions
    %                               'neumann': neumann conditions
    %                               'outflow bc': u_t +/- c*u_x = 0
    %       bc.data             --  cell array with boundary data for left/right boundary
    %                               empty for homogeneous data. Otherwise function data(t,u,ut).
    %   friction   --       struct for friction
    %       friction.fault_fun      -- function @(t,V,Psi)
    %       friction.state_evo_fun  -- function @(Psi,ut)
    %   sources    --      struct for point sources
    %       sources.x           -- array with source positions
    %       sources.source_fun  -- cell array with source functions @(t)
    %
    %   initialcondition    --      struct for initial conditions.
    %       initialcondition.u    -- function @(x,y) for potential condition
    %       initialcondition.ut   -- function @(x,y) for velocity condition
    %
    %   bc_method           --      string for numerical method used to impose boundary conditions
    %                                'standard' or 'erickson2022'
    %
    %   ic_method           --      string for numerical method used to impose nonlinear interface condition
    %                                'standard' or 'erickson2022'

%%% modified by R.F.
    function obj = AntiplaneShear2DSWFrictionAdjDiscr(opset, domain, m, order, material, bc, friction, sources)
    % ↓
    % function obj = AntiplaneShear2DSWFrictionAdjDiscr(opset, domain, m, order, material, bc, friction, sources, fault_sources)
%%%    
        % ---- Default values --------------
        default_arg('opset',@sbp.D2Variable);
        default_arg('m',101);
        default_arg('order',4);
        
        default_arg('material', struct);
        default_field(material, 'rho', {@(x,y) 0*x + 1, @(x,y) 0*x + 1});
        default_field(material, 'mu',  {@(x,y) 0*x + 1, @(x,y) 0*x + 1});
        
        default_arg('sources',[]);
        default_arg('initialcondition',[]);
        %%% added by R.F.
        default_arg('fault_sources',[]);
        %%%
    
        % Get grid and multiblock diff op
        if isequal(opset,@sbp.D2Nonequidistant) % Boundary-optimized operators 
            mbGrid = domain.getGrid(m, 'boundaryopt', order, 'acc');
        else
            mbGrid = domain.getGrid(m, 'equidist');
        end

        nBlocks = mbGrid.nBlocks();
        doParam = cell(nBlocks, 1);
        for i = 1:nBlocks
            if iscell(material.rho)
                doParam{i} = {@(x,y) 1./material.rho{i}(x,y), material.mu{i}, opset};
            elseif isa(rho, 'function_handle')
                doParam{i} = {@(x,y) 1./material.rho(x,y), material.mu, opset};
            else
                error('Inconsistent format for mu and rho.');
            end
        end

        mbDiffOp = multiblock.DiffOp(@scheme.LaplaceCurvilinearNew, mbGrid, order, doParam);
        
        % Spatial operators
        %%% modified by R.F.
        % [D, E, ~, ic_struct] = elastic.discrs.antiplaneshear.discrOps2D(domain, mbGrid, mbDiffOp, material, bc);
        % ↓
        [D, E, ~, ic_struct] = elastic.discrs.antiplaneshear.discrOps2Dsw(domain, mbGrid, mbDiffOp, material, bc);
        %%%

        % Sources
        if strcmp(sources.x,'surface')
            surf = domain.boundaryGroups.surface;
            nSurfBlocks= length(surf);
            blockIds = [];
            src_pos = {};
            for i = 1:nSurfBlocks
                blockId = surf{i}{1};
                pos = mbGrid.getBoundary(surf{i});
                [nSurf,~] = size(pos);
                for j = 1:nSurf
                    src_pos{end+1,1} = pos(j,:);
                    blockIds(end+1,1) = blockId;
                end
            end
            sources.x = src_pos;
            sources.blockIds = blockIds;
        end
        dirac_deltas = elastic.discrs.antiplaneshear.diracDelta2D(sources, mbGrid, order, {opset, opset});

        % ---- Fault interface -----
        Zm = ic_struct.Zm;
        Zp = ic_struct.Zp;
        Wm = ic_struct.Wm;
        Wp = ic_struct.Wp;
        eta = Zp.*Zm./(Zp + Zm);
        % Function that computes tau_l given the entire solution vector
        tau_l_fun = @(U) (Zp.*(Wm*U) - Zm.*(Wp*U)) ./ (Zp + Zm);


         % ---- Initial conditions -----
         [n, ns, nsip, nsim] = elastic.discrs.antiplaneshear.nunknowns2D(domain, mbGrid, mbDiffOp, bc);
        u0 = zeros(n,1);
        v0 = zeros(n,1);
        %%% modified by R.F.
        % Psi0 = zeros(nsim, 1);
        % ↓
        s0 = zeros(nsim, 1);
        %%%
        us0 = zeros(ns,1);
        usim0 = zeros(nsim,1);
        usip0 = zeros(nsip,1);
        %Vs0 = zeros(nsim,1);

        %%% modified by R.F.
        % w0 = [u0; v0; Psi0; us0; usim0; usip0];%; Vs0];
        % ↓
        w0 = [u0; v0; s0; us0; usim0; usip0];%; Vs0];
        %%%

        % ---- Set properties -----
        obj.grid = mbGrid;
        obj.order = order;
        obj.material = material;
        
        % Operators for full domain
        obj.D = D;
        obj.H = mbDiffOp.H;
        obj.E = E;
        obj.mbDiffOp = mbDiffOp;

        % Initial condition
        obj.w0 = w0; 
        
        % Friction properties
        obj.friction = friction;
        obj.e_fault_m = ic_struct.em;
        obj.e_fault_p = ic_struct.ep;
        obj.penalty_fault = ic_struct.penalty;
        obj.eta = eta;
        obj.tau_l_fun = tau_l_fun;
        
        % Point sources
        obj.sources = sources;
        obj.dirac_deltas = dirac_deltas;

        %%% added by R.F.
        obj.bidm = domain.boundaryGroups.fault_minus;
        obj.H_flt = obj.mbDiffOp.getBoundaryQuadrature(obj.bidm);
        obj.flt_pos = obj.grid.getBoundary(obj.bidm);
        %%%
    end
    
    % Prints some info about the discretisation
    function printInfo(obj)
        fprintf('Name: %s\n',obj.name);
        fprintf('Size: %d\n',obj.size());
    end
    
    % Return the number of DOF
    function n = size(obj)
        n = obj.grid.nPoints;
    end
    
    function V = fault_jump(obj, v)
        V = (obj.e_fault_p-obj.e_fault_m)'*v;
    end

    %%% modified by R.F.

    % function obj = setFaultTraction(obj)  
    %     % Friction functions
    %     F_V = obj.friction.data.F_V;
    %     G_V = obj.friction.data.G_V;
        
    %     % Operators 
    %     penalty = obj.penalty_fault;
    %     E = obj.E;
    %     eta = obj.eta;

    %     % Function that computes V* given tau_l and Psi
    %     V_star_from_tau_l = @(i_t, tau_l, Psi) -1./((eta + F_V(:,i_t))).*(tau_l + G_V(:,i_t).*Psi);
    %     obj.V_star = @(i_t, U) V_star_from_tau_l(i_t, obj.tau_l_fun(U), E.Psi*U);
    %     % Fault traction function
    %     obj.fault_traction = @(i_t, U, Vs) penalty*(F_V(:,i_t).*Vs + G_V(:,i_t).*(E.Psi*U));
    % end
        
    % function obj = setStateEvolution(obj)
    %     E = obj.E;
    %     F_Psi = obj.friction.data.F_Psi;
    %     G_Psi = obj.friction.data.G_Psi;
    %     obj.state_evo = @(i_t, U, Vs) E.Psi'*(G_Psi(:,i_t).*(E.Psi*U) + F_Psi(:,i_t).*Vs);
    % end
    
    % ↓

    function obj = setFaultTraction(obj)  
        % Friction functions
        swF_s = obj.friction.data.swF_s;
        defrupt = obj.friction.data.defrupt;
        swFaultVadj = obj.friction.funs.swFaultVadj;
        swFaulttauadj = obj.friction.funs.swFaulttauadj;

        % Operators 
        penalty = obj.penalty_fault;
        E = obj.E;
        eta = obj.eta;

        %%% modified for variation:
        % % % Function that computes V* given tau_l and taustr_adj
        % V_star_from_tau_l = @(i_t, tau_l, str) swFaultVadj(defrupt(:,i_t), str, eta, tau_l); % taustr_adj = s 
        % tau_star_from_tau_l = @(i_t, tau_l, str) swFaulttauadj(defrupt(:,i_t), str, eta, tau_l);
        %↓
        adj_traction = obj.adj_traction;
        has_traction_source = ~isempty(adj_traction);

        % Note that tau_l is negative sign
        if has_traction_source
            V_star_from_tau_l = @(i_t, tau_l, str) swFaultVadj(defrupt(:,i_t), str+ adj_traction(i_t), eta, tau_l); % taustr_adj = s 
            tau_star_from_tau_l = @(i_t, tau_l, str) swFaulttauadj(defrupt(:,i_t), str + adj_traction(i_t), eta, tau_l);
        else
            V_star_from_tau_l = @(i_t, tau_l, str) swFaultVadj(defrupt(:,i_t), str, eta, tau_l); % taustr_adj = s 
            tau_star_from_tau_l = @(i_t, tau_l, str) swFaulttauadj(defrupt(:,i_t), str, eta, tau_l);
        end
        %%%

        obj.V_star = @(i_t, U) V_star_from_tau_l(i_t, obj.tau_l_fun(U), E.s*U);
        obj.fault_traction = @(i_t, U) penalty*(tau_star_from_tau_l(i_t, obj.tau_l_fun(U), E.s*U));
        % obj.V_star = @(i_t, U) swFaultVadj(defrupt(:,i_t), E.s*U, eta, obj.tau_l_fun(U));
        % obj.fault_traction = @(i_t, U) penalty*(swFaulttauadj(defrupt(:,i_t), E.s*U, eta, obj.tau_l_fun(U)));
    end
        
    function obj = setStrengthEvolution(obj)
        E = obj.E;
        swF_s = obj.friction.data.swF_s;
        defrupt = obj.friction.data.defrupt;
        obj.strength_evo = @(i_t, U, Vs) E.s'*(swF_s(:,i_t).*Vs);
        % obj.strength_evo = @(i_t, U, Vs) E.s'*(swF_s(:,i_t).*Vs.*defrupt(:,i_t));
    end

    %%%

    function obj = setFaultSources(obj, fault_sources)
        defaultField(fault_sources, 'x', []);
        % defaultField(fault_sources, 'blockIds', []);
        defaultField(fault_sources, 'data', []);

        % 1. define delta functions on the fault
        % TO DO: deltafunction for curved fault
        flt_pos = obj.flt_pos;
        H_flt = obj.H_flt;        
        [Nflt, ~] = size(flt_pos);
        flt_grd = grid.Cartesian(flt_pos(:,1));

        fault_dirac_deltas = [];
        smoothness_conditions = 0;
        if ~isempty(fault_sources)
            ns = length(fault_sources.x); % number of fault sources
            fault_dirac_deltas = sparse(Nflt, ns);
            for i = 1:ns
                x_i = fault_sources.x{i};
                delta_fun = diracDiscr(flt_grd, x_i, obj.order, smoothness_conditions, H_flt);
                fault_dirac_deltas(:,i) = delta_fun;
            end
        end
        obj.fault_dirac_deltas = fault_dirac_deltas;
        obj.fault_sources = fault_sources;

        % 2. the delta functions f(x) x source function g(t)
        % -> adjoint traction source h(x,t)
        if isempty(fault_sources)
            adj_traction = [];
        else
            fault_source_data = obj.fault_sources.data;
            fault_deltas = obj.fault_dirac_deltas;
            [~,ns] = size(fault_deltas);
            Nt = length(fault_source_data{1});
            data = zeros(ns,Nt);
            % rho = multiblock.evalOn(obj.grid,obj.material.rho);
            % scaledReceivers = obj.E.v'*spdiag(1./rho)*deltas;
            for i = 1:ns
                data_i = fault_source_data{i};
                if iscolumn(data_i)
                    data_i = transpose(data_i);
                end
                data(i,:) = data_i;
            end
            adj_traction = @(i_t) fault_deltas*data(:,i_t);
        end
        obj.adj_traction = adj_traction;
    end

    % set fault sources as a function for x (not delta function)
    function obj = setFaultSourcesFunc(obj, fault_sources)
        defaultField(fault_sources, 'xfunc', []);
        % defaultField(fault_sources, 'blockIds', []);

        if isempty(fault_sources)
            adj_traction = [];
        else
            % 1. Get the functions on the fault
            obj.fault_sources = fault_sources;

            flt_pos = obj.flt_pos;
            % H_flt = obj.H_flt;        
            % [Nflt, ~] = size(flt_pos);
            % flt_grd = grid.Cartesian(flt_pos(:,1));
            if ~isempty(fault_sources.xfunc)
                source_vec = fault_sources.xfunc(flt_pos(:,1));
            end

            % 2. the spatial functions f(x) x source function g(t)
            % -> adjoint traction source h(x,t)
            fault_source_data = fault_sources.data;
            adj_traction = @(i_t) source_vec * fault_source_data(:,i_t);
            % Nt = length(fault_source_data);
            % data = zeros(1,Nt);
            % data_i = fault_source_data{i};
            % if iscolumn(data_i)
            %     data_i = transpose(data_i);
            % end
            % data(i,:) = data_i;
            % adj_traction = @(i_t) source_vec*data(:,i_t);
        end
        obj.adj_traction = adj_traction;
    end

    function obj = setPointSources(obj)
        if isempty(obj.sources)
            F = [];
        else
            source_data = obj.sources.data;
            deltas = obj.dirac_deltas;
            [~,ns] = size(deltas);
            Nt = length(source_data{1});
            data = zeros(ns,Nt);
            rho = multiblock.evalOn(obj.grid,obj.material.rho);
            scaledReceivers = obj.E.v'*spdiag(1./rho)*deltas;
            for i = 1:ns
                data_i = source_data{i};
                if iscolumn(data_i)
                    data_i = transpose(data_i);
                end
                data(i,:) = data_i;
            end
            F = @(i_t) scaledReceivers*data(:,i_t);
        end
        obj.S = F;
    end

            
    % Assembles the forcing terms to a single function F(t, U)
    %%% modified by R.F.
    % function F_tot = assembleForcing(obj)
    %     has_source = ~isempty(obj.S); % TBD: Keep this? Should always be true.
        
    %     % Assemble rhs F(i_t, U, Vs)
    %     if has_source
    %         F = @(i_t, U, Vs) obj.fault_traction(i_t, U, Vs) + obj.state_evo(i_t, U, Vs) + obj.S(i_t);
    %     else
    %         F = @(i_t, U, Vs) obj.fault_traction(i_t, U, Vs) + obj.state_evo(i_t, U, Vs);
    %     end
        
    %     % % Compute V_star once, pass along to other forcing terms
    %     % % and store in returned vector r
    %     % function r = F_precomp_Vs(i_t,U)
    %     %     Vs = obj.V_star(i_t, U);
    %     %     r = obj.E.Vs'*Vs + F(i_t, U, Vs);
    %     % end
    %     % F_tot = @F_precomp_Vs;
    %     % Compute V_star once, pass along to other forcing terms
    %     F_tot = @(t, U) F(t, U, obj.V_star(t, U));
    % end

    % ↓

    function F_tot = assembleForcing(obj)
        has_source = ~isempty(obj.S);
        
        % Assemble rhs F(i_t, U, Vs)
        if has_source
            F = @(i_t, U, Vs) obj.fault_traction(i_t, U) + obj.strength_evo(i_t, U, Vs) + obj.S(i_t);
        else
            F = @(i_t, U, Vs) obj.fault_traction(i_t, U) + obj.strength_evo(i_t, U, Vs);
        end
        
        F_tot = @(t, U) F(t, U, obj.V_star(t, U));
    end 
    %%%
    
    % Returns a timestepper for integrating the discretisation in time
    %     method is a string that states which timestepping method should be used.
    %          The implementation should switch on the string and deliver
    %          the appropriate timestepper. It should also provide a default value.
    %     time_align is a time that the timesteps should align with so that for some
    %                integer number of timesteps we end up exactly on time_align
    function [ts, N] = getTimestepper(obj, method, time_align, k)
        default_arg('method','struct');
        default_field(method,'order',4);
        default_arg('time_align',[]);
        default_arg('k',[]);
        
        t = 0;
        switch method.order
        case 3
            cfl = 0.1;
        case 4
            cfl = 0.1;
        case 5
            cfl = 0.1;
        case 6
            cfl = 0.1;
        otherwise
            error('RK of order %d not implemented',method.order);
        end
        
        if isempty(k)
            k = obj.getTimestep(method,cfl);
        end
        
        if ~isempty(time_align)
            [k, N] = alignedTimestep(k, time_align);
        end
        F = obj.assembleForcing();
        
        ts = time.ExplicitRungeKuttaDiscreteData(obj.D, [], F, k, t, obj.w0, method.order);
    end
    
    
    function k = getTimestep(obj, method, cfl)
        default_arg('cfl',[]);
        if isempty(cfl)
            error('Specify CFL')
        end
        
        nBlocks = obj.grid.nBlocks;
        k = inf;
        for i = 1:nBlocks
            a = diag(obj.mbDiffOp.diffOps{i}.a);
            b = diag(obj.mbDiffOp.diffOps{i}.b);

            % Physical wave speed
            v = sqrt(a.*b);

            % Reference grid spacing
            h = obj.mbDiffOp.diffOps{i}.h;
            h_xi = h(1);
            h_eta = h(2);

            x_xi = obj.mbDiffOp.diffOps{i}.x_u;
            x_eta = obj.mbDiffOp.diffOps{i}.x_v;
            y_xi = obj.mbDiffOp.diffOps{i}.y_u;
            y_eta = obj.mbDiffOp.diffOps{i}.y_v;

            % Approximate physical grid spacings
            h1 = sqrt( (x_xi*h_xi).^2 + (y_xi*h_xi).^2 );
            h2 = sqrt( (x_eta*h_eta).^2 + (y_eta*h_eta).^2 );

            k1 = cfl*min(h1./v);
            k2 = cfl*min(h2./v);

            % Time-step suggested by current block
            k_temp = min(k1,k2);

            % Use smallest time steps of all blocks
            k = min(k, k_temp);
        end
    end
    
    function r = getTimeSnapshot(obj, ts)
        if ts == 0
            r.t = 0;
            r.i_t = 1;
            r.u = obj.E.u*obj.w0;
            r.v = obj.E.v*obj.w0;
            %%% modified by R.F.
            % r.Psi = obj.E.Psi*obj.w0;
            % ↓
            r.s = obj.E.s*obj.w0;
            %%%
            r.us = obj.E.us*obj.w0;
            r.usim = obj.E.usim*obj.w0;
            r.usip = obj.E.usip*obj.w0;
            %r.Vs = obj.E.Vs*obj.w0;
            return
        end
        r.t = ts.t;
        r.i_t = ts.n*ts.s; % TODO: Make sure that r.i_t is the correct time index
        w = ts.getV();
        r.u = obj.E.u*w;
        r.v = obj.E.v*w;
        %%% modified by R.F.
        % r.Psi = obj.E.Psi*w;
        % ↓
        r.s = obj.E.s*w;
        %%%
        r.usim = obj.E.usim*w;
        r.usip = obj.E.usip*w;
        r.us = obj.E.us*w;
        %r.Vs = obj.E.Vs*w;
    end
    
    % Sets up a plot of the discretisation
    %     update is a function_handle accepting a timestepper that updates the plot to the
    %            state of the timestepper
    function [update,figure_handle] = setupPlot(obj, type)
        figure_handle = figure();
        x = obj.grid.points();
        default_arg('type',struct);
        default_field(type,'plot_variables','uv')
        
        function S = plot_comp(comp)
            S = multiblock.Surface(obj.grid, comp);
            xlabel('$x$','interpreter','latex');
            ylabel('$y$','interpreter','latex');
            elastic.cmap;
            colorbar;
            shading interp;
        end
        
        function Surs = setup_uv(r)
            subplot(2,1,1);
            Sur1 = plot_comp(r.u);
            title(sprintf('u, t = %.2e', r.t));
            
            subplot(2,1,2);
            Sur2 = plot_comp(r.v);
            title(sprintf('v'));

            Surs = {Sur1, Sur2};
        end

        function update_uv(r, Surs)
            Surs{1}.ZData = r.u;
            Surs{1}.CData = r.u;
            subplot(2,1,1);
            title(sprintf('u, t = %.2e', r.t));
            Surs{2}.ZData = r.v;
            Surs{2}.CData = r.v;
        end
        
        function [Surs, Lines] = setup_all(r)
            % V_plot = [V_plot, obj.slipRate(r.v)];
            subplot(4,1,1);
            Sur1 = plot_comp(r.u);
            title(sprintf('u, t = %.2e', r.t));
            
            subplot(4,1,2);
            Sur2 = plot_comp(r.v);
            title(sprintf('v'));

            subplot(4,1,3);
            %%% modified by R.F.
            % line1 = plot(r.Psi, 'linewidth',  2);
            % ylabel('$\Psi$','interpreter','latex');
            % ↓
            line1 = plot(r.s, 'linewidth',  2);
            ylabel('$\tau^{\dagger}_{str}$','interpreter','latex');
            %%%

            xlabel('number along fault');

            subplot(4,1,4);
            %%% modified by R.F.
            % ylab = '$\log(|V^*|)$';
            % ↓
            ylab = '$V^*$';
            %%%
            % TODO: Make sure that r.i_t is the correct time index

            %%% modified by R.F.
            % w = [r.u; r.v; r.Psi; r.us; r.usim; r.usip];
            % ↓
            w = [r.u; r.v; r.s; r.us; r.usim; r.usip];
            %%%
            Vs = obj.V_star(r.i_t, w);
            %%% modified by R.F.           
            % line2 = plot(log10(abs(Vs)), 'linewidth',  2);
            % ↓
            line2 = plot(Vs, 'linewidth',  2);
            %%%
            xlabel('number along fault');
            ylabel(ylab,'interpreter','latex');

            Surs = {Sur1, Sur2};
            Lines = {line1, line2};
        end

        function update_all(r, Surs, Lines)
            subplot(4,1,1);
            Surs{1}.ZData = r.u;
            Surs{1}.CData = r.u;
            cmax = max(abs(r.u)) + 1e-10;
            caxis([-cmax, cmax]);
            title(sprintf('u, t = %.2e', r.t));

            subplot(4,1,2);
            Surs{2}.ZData = r.v;
            Surs{2}.CData = r.v;
            cmax = max(abs(r.v)) + 1e-10;
            caxis([-cmax, cmax]);

            subplot(4,1,3);
            %%% modified by R.F.
            % Lines{1}.YData = r.Psi; 
            % ↓
            Lines{1}.YData = r.s;
            %%%

            subplot(4,1,4);
            % TODO: Make sure that r.i_t is the correct time index
            %%% modified by R.F.
            % w = [r.u; r.v; r.Psi; r.us; r.usim; r.usip];
            % ↓
            w = [r.u; r.v; r.s; r.us; r.usim; r.usip];
            %%%
            Vs = obj.V_star(r.i_t, w);
            %%% modified by R.F.
            % Lines{2}.YData = log10(abs(Vs));
            % ↓
            Lines{2}.YData = Vs;
            %%%
        end

        function plot_with_traction_trajectory(r)
            error('Not implemented');
        end

        
        switch type.plot_variables
        case 'uv'
            r = obj.getTimeSnapshot(0);
            Surs = setup_uv(r);
            update = @(r)update_uv(r, Surs);
        case 'all'
            r = obj.getTimeSnapshot(0);
            [Surs, Lines] = setup_all(r);
            update = @(r)update_all(r, Surs, Lines);
        end
    end
    
    % Compare two functions u and v in the discrete l2 norm.
    function e = compareSolutions(obj, u, v)
        e = [];
        error('Not implemented')
    end
    
end
end