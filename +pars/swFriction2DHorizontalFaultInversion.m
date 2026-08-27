function [parameters, trueParameters] = swFriction2DHorizontalFaultInversion(opts)
    default_arg('opts',struct);
    default_field(opts,'m',251);
    default_field(opts,'m_p',26);
    default_field(opts,'k',0.003);
    default_field(opts,'misfitType','velocity');
    default_field(opts,'receiverSpacing',2);
    %%% modified by R.F.
    % default_field(opts,'inversionParameters',{'a'});
    % ↓
    default_field(opts,'inversionParameters',{'tau0'});
    %%%
    default_field(opts,'initialGuessScalings',[]);
    default_field(opts,'initialGuessValues',[]);
    default_field(opts,'doFilter',false);
    %%% added by R.F.
    default_field(opts,'swParams',[]);
    default_field(opts,'faultGeometry', []);
    %%%

    if isempty(opts.initialGuessScalings) && isempty(opts.initialGuessValues)
        error();
    end

    %% SBP operators and grid points
    opset = @sbp.D2Nonequidistant;
    % Grid points along the fault. The total number of gridpoints is
    % approximately [m,m]^2
    m = opts.m;
    m_p = opts.m_p;

    %% Domain
    T = 6;  % Final time
    xlims = [-10, 10];
    ylims = [-10, 0, 10];
    %%% added by R.F.
    if ~isempty(opts.faultGeometry)
        xlims = opts.faultGeometry.xlims;
    end
    %%%
    domain = domains.HorizontalFault(xlims, ylims);

    %% Domain material parameters
    % material.rho - cell array with function handles of (x,y) for density for each block
    % material.mu - cell array with function handles (x,y) for bulk modulus for each block
    cs = 3.464; % km/s
    rho = 2.67; % g/cm3
    mu = cs^2*rho; % GPa

    material.rho = {@(x,y) 0*x+rho, @(x,y) 0*x+rho};
    material.mu = {@(x,y) 0*x+mu, @(x,y) 0*x+mu};

    
    %% Boundary conditions
    % Store boundary ids
    b = domain.boundaryGroups.bottom;
    t = domain.boundaryGroups.surface;
    l = domain.boundaryGroups.left;
    r = domain.boundaryGroups.right;
    bc.ids = {l, r, b, t};
    
    % Set types of conditions and potential data
    % bc.type - cell array with entries 'dirichlet', 'traction', or 'outflow'
    %           for each boundary.
    % bc.data - cell array with function handles to time-dependent data. 
    %           May be specified for dirichlet or traction conditions.
    %           For zero data set entries to [].
    %           ex: bc.type = {'outflow', 'outflow', 'outflow', 'traction'};
    %               bc.data = {[],[],[],@(t) g(t)};
    bc.type = {'outflow', 'outflow', 'outflow', 'outflow'};
    bc.data = {[],[],[],[]};
    
    %% Point sources
    % Set to empty if no point sources are included. Otherwise specify
    % sources.x - cell array with coordinate vectors for receiver positions
    %             ex: sources.x = {[x1,y1], [x2,y2], ...}
    % sources.blockIds - vector with grid block ids for each source. 
    %                    See domain for detail on block ids.
    %                    ex: sources.blockIds = [blockId1, blockId2, ...] 
    % sources.funs - cell array with function handles for time dependent part of sources
    %                ex: sources.funs = {@(t) f1(t), @(t) f2(t), ...}
    sources = [];

    %%% added by R.F.
    %%% Fault sources
    fault_sources = struct;
    fault_sources.x = {};
    fault_sources.blockIds = [];
    % domain.boundaryGroups.fault
    %%%

    %% Receivers
    % receivers.x - cell array with coordinate vectors for receiver positions
    %               or character string 'surface' to use all grid points on the surface.
    %               ex: receivers.x = {[x1,y1], [x2,y2], ...}
    % receivers.blockIds - vector with grid block ids for each source. Set to []
    %                      if receivers.x = 'surface'
    %                      See domain for detail on block ids.
    %                      ex: sources.blockIds = [blockId1, blockId2, ...] 
    boundaries.outer.x = [-9, 9];
    boundaries.outer.y = [-9, 9];
    boundaries.inner.y = [-2, 2];
    if opts.receiverSpacing == 1
        boundaries.inner.x = [-6, 7];
    elseif opts.receiverSpacing == 2
        boundaries.inner.x = [-7, 7];
    end
    receivers = pars.placeReceivers(boundaries,opts.receiverSpacing);
    
    
    
    %% Friction
    %%% modified by R.F.

    % % Parameters
    % % Must set a, b, V0, f0, sigma0, D_c, Psi0
    % % These are also the variables which may be inverted for
    % % The parameters can be either scalars or function handles of
    % % (x,y).
    % 
    % if_velocity_weakening = @(x,y) (x >= - 5) & (x <= 6);
    % if_velocity_strengthening = @(x,y) 1-if_velocity_weakening(x,y);
    % 
    % % 'True' parameter values, used to generate data for misfit
    % a =  @(x,y) 0.009*if_velocity_weakening(x,y) + 0.013*if_velocity_strengthening(x,y);
    % b = 0.011;
    % sigma0 = 120;
    % f0 = 0.6;
    % V0 = 1e-6;
    % D_c = @(x,y) (0.2*if_velocity_weakening(x,y) + 1*if_velocity_strengthening(x,y));
    % Psi0 = 0.724339595021678;
    % tau0 = 72; % Note: Really sigma_yz. Will be dotted with the normal component on the fault!
    % 
    % % Set initial velocity and estimated final velocity
    % % and compute initial state, and background traction
    % A = 25;
    % Vinit = 1e-12;
    % 
    % 
    % % Set functions for friction law and state law
    % % These functions are defined in elastic.friction.inversion.
    % % See elastic.friction.inversion.generateFunctions for details
    % % on how to generate the friction law and state law of interest.
    % F = @elastic.friction.inversion.F;
    % F_V = @elastic.friction.inversion.F_V;
    % F_Psi = @elastic.friction.inversion.F_Psi;
    % G = @elastic.friction.inversion.G;
    % G_V = @elastic.friction.inversion.G_V;
    % G_Psi = @elastic.friction.inversion.G_Psi;
    % % Functions used to solve for characteristic variables
    % % Note: Must match with F above! 
    % % TODO: Should probably be generated together with above functions
    % Finv = @elastic.friction.erickson2022.asinh_inv;
    % nonlin_solve_fun = @elastic.friction.erickson2022.nonlin_solve_fun_asinh;
    % % External loading
    % x0 = 3;
    % %y0 = -8.955116271972656e-01;
    % std_dev = 2;
    % tau_L = @(t, x, y) A*gaussian1d(x,x0,std_dev);
    % 
    % % Pack parameters
    % friction.rsParams.a = a;
    % friction.rsParams.b = b;
    % friction.rsParams.sigma0 = sigma0;
    % friction.rsParams.f0 = f0;
    % friction.rsParams.V0 = V0;
    % friction.rsParams.D_c = D_c;
    % friction.rsParams.tau0 = tau0;
    % 
    % friction.loadingParams.Vinit = Vinit; 
    % friction.loadingParams.A = A;
    % 
    % % Pack functions
    % friction.funs.F = F;
    % friction.funs.F_V = F_V;
    % friction.funs.F_Psi = F_Psi;
    % friction.funs.G = G;
    % friction.funs.G_V = G_V;
    % friction.funs.G_Psi = G_Psi;
    % friction.funs.Finv = Finv;
    % friction.funs.nonlin_solve_fun = nonlin_solve_fun;
    % friction.funs.tau_L = tau_L;
        
    % ↓

    % Parameters
    % Must set taur, taup, tau0, D_c, V0, s0, sigma0, 
    % These are also the variables which may be inverted for
    % The parameters can be either scalars or function handles of
    % (x,y).

    %%% Parameter setting 1
    % if_zero_stressdrop = @(x,y) (x < -2);
    % if_positive_stressdrop = @(x,y) (x >= -2) & (x <= 5);
    % if_negative_stressdrop = @(x,y) (x > 5);
    % if_nucleation = @(x,y) (x >= -2) & (x <= -1.5);
    % % if_pos_strengthexcess = @(x,y) (x > 5);

    % % 'True' parameter values, used to generate data for misfit
    % tau0 = 72; % Note: Really sigma_yz. Will be dotted with the normal component on the fault!
    % taur = @(x,y) (0*if_zero_stressdrop(x,y) - 6*if_positive_stressdrop(x,y) + 2*if_negative_stressdrop(x,y)+tau0) ;
    % taup = @(x,y) (0*if_zero_stressdrop(x,y) - 6*if_nucleation(x,y) + 2*if_positive_stressdrop(x,y) + 10*if_negative_stressdrop(x,y)+tau0) ;
    % D_c = @(x,y) (0.2*if_zero_stressdrop(x,y) + 0.2*if_positive_stressdrop(x,y) + 3*if_negative_stressdrop(x,y));

    %%% Parameter setting 2
    % if_nucleation = @(x,y) (x >= 2.5) & (x <= 3.5);
    % if_barrier = @(x,y)  (x <= 0) | (x >= 5);
    % tau0 = 72; % Note: Really sigma_yz. Will be dotted with the normal component on the fault!
    % taur = @(x,y) (66 + 7 * if_barrier(x,y));
    % taup = @(x,y) (74 - 8 * if_nucleation(x,y));
    % D_c = @(x,y) (0.2 + 1.8 * if_barrier(x,y));

    %%% Default Parameter setting: Overwitten by opts.swParams
    if_nucleation = @(x,y) (x >= 2.5) & (x <= 3.5);
    if_barrier = @(x,y)  (x <= -2) | (x >= 5);
    tau0 = 72; 
    taur = @(x,y) (66 + 7 * if_barrier(x,y));
    taup = @(x,y) (74 - 8 * if_nucleation(x,y));
    D_c = 0.2;
    sigma0 = 120;

    % Set initial velocity and estimated final velocity
    % and compute initial state, and background traction
    A = 25;
    Vinit = 0;
    sinit = 0;
    
    % Set functions for strength
    % These functions are defined in elastic.friction.swinversion.
    % Fucntions for forward
    swFaultV = @elastic.friction.swinversion.swFaultV;
    swFaulttau = @elastic.friction.swinversion.swFaulttau;
    dsdt = @elastic.friction.swinversion.dsdt;

    % Functions for adjoint
    swF_s = @elastic.friction.swinversion.swF_s;
    defrupt = @elastic.friction.swinversion.defrupt;
    swFaultVadj = @elastic.friction.swinversion.swFaultVadj;
    swFaulttauadj = @elastic.friction.swinversion.swFaulttauadj;
    swF = @elastic.friction.swinversion.swF; % strength: not used

    % External loading
    % x0 = 3;
    % %y0 = -8.955116271972656e-01;
    % std_dev = 2;
    % tau_L = @(t, x, y) A*gaussian1d(x,x0,std_dev); 
    
    % Pack parameters
    %%% modified by R.F.
    if ~isempty(opts.swParams)
        friction.swParams = opts.swParams;
        if isfield(opts.swParams, 'tau_L')
            friction.funs.tau_L = opts.swParams.tau_L;
            friction.swParams = rmfield(friction.swParams, 'tau_L');
        end
    else
        friction.swParams.taur = taur;
        friction.swParams.taup = taup;
        friction.swParams.sigma0 = sigma0;
        friction.swParams.D_c = D_c;
        friction.swParams.tau0 = tau0;
    end
    %%%

    friction.loadingParams.Vinit = Vinit; 
    friction.loadingParams.sinit = sinit; 
    friction.loadingParams.A = A;

    % Pack functions
    friction.funs.swFaultV = swFaultV;
    friction.funs.swFaulttau = swFaulttau;
    friction.funs.dsdt = dsdt;
    
    friction.funs.defrupt = defrupt;
    friction.funs.swF = swF;
    friction.funs.swF_s = swF_s;
    friction.funs.swFaultVadj = swFaultVadj;
    friction.funs.swFaulttauadj = swFaulttauadj;

    %%% modified by R.F.
    % friction.funs.tau_L = tau_L;
    %%%

    %%%

   %% Initial conditions
    % Initial values of fields u, v as a cell array
    % of grid functions for each block.
    % Set initialconditions = [] for zero initial conditions
    initialconditions.u = {@(x,y) 0*x, @(x,y) 0*x};
    initialconditions.v = {@(x,y) 0*x - Vinit/2, @(x,y) 0*x + Vinit/2};

    %%% modified by R.F.
    % initialconditions.Psi = @(x,y) 0*x + Psi0;
    % ↓
    initialconditions.s = @(x,y) 0*x + sinit;
    %%%

    %% ERK time-stepper options.
    % For gradient computations 
    % Either use order = 4, adaptive_fwd = false;
    % or         order = 3, adaptive_fwd = true/false;
    order = 4;
    adaptive_fwd = false;
    tol = 1e-5; % % Only used if adaptive_fwd = true;
    reportRetry = false; % Only used if adaptive_fwd = true;
    k = opts.k; % Time step. Will be aligned to final time. Set to empty if default CFL conditions are used
            % (specified in the AntiplaneShear2DRS<Fwd/Adj>Discr.m)

    if opts.doFilter
        h = 20/(m-1);
        samplingFreq = 1/k;
        cornerFreq = cs/(16*h);
        filterOpts.cutoffFreq = cornerFreq/(samplingFreq/2);
        filterOpts.order = 6;
    else 
        filterOpts = [];
    end

    % Pack tsOpts struct
    tsOpts.forwardMethod.order = order;
    tsOpts.forwardMethod.adaptive = adaptive_fwd;
    tsOpts.forwardMethod.rtol = tol;
    tsOpts.forwardMethod.reportRetry = reportRetry;
    tsOpts.adjointMethod.order = order; 
    tsOpts.k = k;
    
    % Form parameter struct
    parameters = struct;
    parameters.domain = domain;
    parameters.opset = opset;
    parameters.T = T;
    parameters.material = material;
    parameters.bc = bc;
    parameters.friction = friction;
    parameters.sources = sources;
    %%% Added by R.F.
    parameters.fault_sources = fault_sources;
    %%%
    parameters.receivers = receivers;
    parameters.initialconditions = initialconditions;
    parameters.misfitType = opts.misfitType;
    parameters.tsOpts = tsOpts;
    parameters.dim = 2;
    parameters.m = m;
    parameters.m_p = m_p;
    parameters.filterOpts = filterOpts;
    
    trueParameters = parameters;
    % Perturb true parameters
    if ~isempty(opts.initialGuessScalings)
        for i = 1:numel(opts.inversionParameters)
            paramName = opts.inversionParameters{i};
            %%% modified by R.F.
            % true_param = friction.rsParams.(paramName);
            % ↓
            true_param = parameters.friction.swParams.(paramName);
            %%%
            amp = opts.initialGuessScalings(i);
            if isa(true_param,'function_handle')
                perturbed_param = @(x,y) amp*true_param(x,y);
            else
                perturbed_param = amp*true_param;
            end
            %%% modified by R.F.            
            % parameters.friction.rsParams.(paramName) = perturbed_param;
            % ↓
            parameters.friction.swParams.(paramName) = perturbed_param;
            %%%
        end
    else
        for i = 1:numel(opts.inversionParameters)
            paramName = opts.inversionParameters{i};
            %%% modified by R.F.
            % parameters.friction.rsParams.(paramName) = opts.initialGuessValues(i);
            % ↓
            if isa(opts.initialGuessValues{i},'function_handle')
                parameters.friction.swParams.(paramName) = opts.initialGuessValues{i};
            else
                parameters.friction.swParams.(paramName) = opts.initialGuessValues{i};
            end
            % parameters.friction.swParams.(paramName) = opts.initialGuessValues(i);
            %%%
        end
    end
end

% TODO: Check that gaussian matches Erics impl.
function f = gaussian1d(x, x0, d)
    f = exp(-0.5*(x-x0).^2/d^2);
end