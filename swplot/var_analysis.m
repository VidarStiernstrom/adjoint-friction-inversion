% clearvars;
% close all;
% % clf;
% addSubpaths;
% addpath('./swadjoint');

% name = 'var4';
% outpath = './swadjoint/out/';
% savedir = sprintf('%s%s', outpath, name);

function var_analysis(savedir, calx, plotOpts)

    default_arg('calx', [-2, 9]);
    default_arg('plotOpts', struct());
    default_field(plotOpts, 'xlim', [-2, 9]);
    default_field(plotOpts, 'ylim', [0, 0.5]);
    default_field(plotOpts, 'yloglim', [1e-3, 10]);
    default_field(plotOpts, 'yndlim', [1e-7, 1e-4]);
    default_field(plotOpts, 'fontsize', 14);
    default_field(plotOpts, 'varPlotIter', []);
    default_field(plotOpts, 'varPlotOneIter', []);
    default_field(plotOpts, 'RejectDeltaSensitivity', false);
    default_field(plotOpts, 'CohesiveAnalysis', false);
    
    load(sprintf('%s/vars.mat', savedir));

    fontsize = plotOpts.fontsize;
    if isempty(plotOpts.varPlotIter)
        varPlotIter = 1:length(vars.gradients);
    else
        varPlotIter = plotOpts.varPlotIter;
    end

    If2c = vars.If2c;
    H_b = vars.H_b;
    dirac_deltas = vars.dirac_deltas;

    % unit: MPa / km / s
    cs = 3.464; % km/s
    rho = 2.67; % g/cm3
    mu = cs^2*rho*1e3; % MPa
    % mu = cs^2*rho; % GPa
    % cs = 3.464*1e3; % m/s
    % rho = 2.67*1e3; % kg/m3
    % mu = cs^2*rho*1e-6; % MPa

    fault_x = full(vars.fault_x); % m
    Nx = length(fault_x);

    ParamFun = vars.swParams;
    ParamName = fieldnames(ParamFun);
    ParamVec = struct();

    % set ParamFun (function of x only)
    for j = 1:length(ParamName)
        parName = ParamName{j};
        par = ParamFun.(parName);
        % if not function handle
        if ~isa(par, 'function_handle')
            ParamFun.(parName) = @(x) par + x*0;
        else
            ParamFun.(parName) = @(x) par(x,0);
        end
        ParamVec.(parName) = ParamFun.(parName)(fault_x);
    end
    % convert m to km
    ParamFun.D_c = @(x) ParamFun.D_c(x) * 1e-3;
    ParamVec.D_c = ParamVec.D_c * 1e-3;


    %% LEFM framework

    % calx = [-0.7,9];

    ParamFun.fracE = @(x) (ParamFun.taup(x) - ParamFun.taur(x)) .* ParamFun.D_c(x) / 2; 
    ParamFun.deltau = @(x) (ParamFun.tau0(x) - ParamFun.taur(x));
    ParamFun.SIF = @(a) arrayfun(@(a) sqrt(2/pi) * integral(@(x) ParamFun.deltau(x) ./ sqrt(a - x), fault_x(1), a), a);
    % ParamFun.SIF = @(a) arrayfun(@(a) sqrt(2/pi) * integral(@(x) ParamFun.deltau(x) ./ sqrt(a - x), 0, a), a);

    % approx. version
    % ParamFun.Vr = @(x) cs*(1 - 2*mu*ParamFun.fracE(x)./(ParamFun.SIF(x).^2));
    % more correct version
    ParamFun.Vr = @(x) cs*(2 ./ ((2*mu*ParamFun.fracE(x)./ ParamFun.SIF(x).^2).^2 + 1) - 1);
    % ParamFun.Rc = @(x) 9*pi/16 * mu * ParamFun.fracE(x) ./ (ParamFun.taup(x) - ParamFun.taur(x)).^2 .* sqrt(1 - (ParamFun.Vr(x)/cs).^2);
    ParamFun.W = @(x) (ParamFun.taup(x) - ParamFun.taur(x)) ./ ParamFun.D_c(x);
    ParamFun.Rc = @(x) pi/4 * mu ./ ParamFun.W(x) .* sqrt(1 - (ParamFun.Vr(x) / cs).^2);

    % ParamFun.Trup = @(x, x0) arrayfun(@(x, x0) integral(@(xsi) 1./ParamFun.Vr(xsi), x0, x), x, x0);

    %%% Calculate Trup and Vr in LEFM
    fault_xp = fault_x(fault_x > calx(1) & fault_x < calx(2));
    Vr_vec = ParamFun.Vr(fault_xp);
    rupIters = find(Vr_vec > 0);
    fault_xp = fault_xp(rupIters); % only plot where Vr > 0
    Vr_vec = Vr_vec(rupIters);
    Lgriff = fault_xp(1);
    disp(sprintf('Griffith length: %.2f km\n', Lgriff));
    % Trup_vec = ParamFun.Trup(fault_xp(), fault_xp(1)*ones(size(fault_xp)));
    Trup_vec = zeros(size(fault_xp));
    for i = 2:length(fault_xp)
        Trup_vec(i) = (1 / Vr_vec(i-1) + 1 / Vr_vec(i)) * (fault_xp(i) - fault_xp(i-1)) / 2 + Trup_vec(i-1);
    end

    Rc_vec = ParamFun.Rc(fault_xp);

    %% Numerical 

    %%% Get Trup and Vr from numerical simulation
    Trup_Num = vars.Tdc;
    rupIterNum = find(~isnan(Trup_Num) & fault_x > calx(1) & fault_x < calx(2));

    %%% Get cohesive zone size
    % if Trup is no defined, skip
    if ~isfield(vars, 'Trup')
        disp('Trup is not defined, skip cohesive zone size calculation');
    else
        Trup0 = vars.Trup;
        % load(sprintf('%s/Trup.mat', savedir), 'Trup');
        % Trup0 = Trup;
        Rc_Num = zeros(size(rupIterNum));
        for i = 1:length(rupIterNum)
            idx = rupIterNum(i);
            Tdc_xi = Trup_Num(idx);
            xi_plus_Rc_iter = find(Trup0(:,1) >= Tdc_xi, 1, 'first');
            Rc_Num(i) = fault_x(xi_plus_Rc_iter) - fault_x(idx);
        end
    end

    % smoothing by spline function
    % smoothing_param = 0.65;
    % fitobj = fit(fault_x(rupIterNum), Trup_Num(rupIterNum), 'smoothingspline', 'SmoothingParam', smoothing_param);
    % dT_dx = differentiate(fitobj, fault_x(rupIterNum));
    % Trup_num_approx = feval(fitobj, fault_x(rupIterNum));
    % Vr_Num = 1 ./ dT_dx;

    Vr_Num = zeros(size(rupIterNum));
    for i = 1:length(rupIterNum)
        idx = rupIterNum(i);
        if i == 1
            Vr_Num(i) = (fault_x(idx+1) - fault_x(idx)) / (Trup_Num(idx+1) - Trup_Num(idx));
        elseif i == length(rupIterNum)
            Vr_Num(i) = (fault_x(idx) - fault_x(idx-1)) / (Trup_Num(idx) - Trup_Num(idx-1));
        else
            Vr_Num(i) = (fault_x(idx+1) - fault_x(idx-1)) / (Trup_Num(idx+1) - Trup_Num(idx-1));
        end
    end

    %% Vr-plot
    f = figure();
    plot(fault_xp, Vr_vec / cs, 'linewidth', 2, 'DisplayName', 'LEFM');
    hold on;
    plot(fault_x(rupIterNum), Vr_Num / cs, 'linewidth', 2, 'DisplayName', 'Num.');
    xlabel('x along fault (km)');
    ylabel('V_r/c_s');
    ylim([0, 1.1]);
    grid on;
    legend('Location', 'best');
    set(gca,'fontsize',fontsize);
    pdf_output(f, sprintf('%s/LEFM_Vr_', savedir));
    xlim([Lgriff, calx(2)]);
    pdf_output(f, sprintf('%s/LEFM_Vr', savedir));

    %% Trup plot
    f2 = figure();
    plot(fault_xp, Trup_vec, 'linewidth', 2, 'DisplayName', 'LEFM');
    hold on;
    plot(fault_x(rupIterNum), Trup_Num(rupIterNum) - Trup_Num(rupIterNum(1)), 'k--', 'linewidth', 2, 'DisplayName', 'Num');
    % plot(fault_x(rupIterNum), Trup_num_approx - Trup_num_approx(1), 'r--', 'linewidth', 2, 'DisplayName', 'Num Approx');
    xlabel('x along fault (km)');
    ylabel('T_{rup} (s)');
    xlim(calx);
    legend('Location', 'best');
    grid on;
    set(gca,'fontsize',fontsize);
    pdf_output(f2, sprintf('%s/Trup', savedir));

    %% Rc plot
    f3 = figure();
    plot(fault_xp, Rc_vec, 'linewidth', 2, 'DisplayName', 'LEFM');
    hold on;
    if isfield(vars, 'Trup')
        plot(fault_x(rupIterNum), Rc_Num, 'k--', 'linewidth', 2, 'DisplayName', 'Num');
    end
    xlabel('x along fault (km)');
    ylabel('R_c (km)');
    xlim(calx);
    grid on;
    set(gca,'fontsize',fontsize);
    pdf_output(f3, sprintf('%s/Rc', savedir));

    %% vTdc plot
    vTdc_Num = vars.vTdc; % m/s
    vTdc_pred = ParamFun.D_c(fault_xp) .* Vr_vec ./ ParamFun.Rc(fault_xp); % km/s

    f4 = figure();
    plot(fault_x(rupIterNum), vTdc_Num(rupIterNum), 'linewidth', 2, 'DisplayName', 'Num');
    hold on;
    plot(fault_xp, vTdc_pred * 1e3, 'k--', 'linewidth', 2, 'DisplayName', 'Pred');
    xlabel('x along fault (km)');
    ylabel('v_{Tdc} (m/s)');
    xlim(calx);
    grid on;
    legend('Location', 'best');
    set(gca,'fontsize',fontsize);
    pdf_output(f4, sprintf('%s/vTdc', savedir));



    %% Calculate Sensitivity for Trup in LEFM
    gradLEFM = struct();
    % approx. version
    % gradLEFM.fun.fracE = @(x) (cs - ParamFun.Vr(x)) ./ (ParamFun.Vr(x).^2) ./ ParamFun.fracE(x);
    % more correct version
    gradLEFM.fun.fracE = @(x) cs * (1 - (ParamFun.Vr(x) / cs).^2) ./ (ParamFun.Vr(x).^2) ./ ParamFun.fracE(x);
    gradLEFM.fun.fracE_ = @(x) cs * (1 - (ParamFun.Vr(x) / cs).^2) ./ (ParamFun.Vr(x).^2);
    gradLEFM.fun.fracE_cohesive = @(x) 1 / cs + x.* 0;
    gradLEFM.fun.W_cohesive = @(x) - 1 ./ ParamFun.Vr(x);
    % gradLEFM.fun.fracE_cohesive = @(x) (1 - cs ./ (2 * ParamFun.Vr(x))) / cs;
    % gradLEFM.fun.W_cohesive = @(x) - 1 ./ (2 * ParamFun.Vr(x));

    %%% Continuous to Discrete gradient
    gradLEFM.vec.fracE = H_b * If2c * gradLEFM.fun.fracE(fault_x);
    gradLEFM.vec.fracE_ = H_b * If2c * gradLEFM.fun.fracE_(fault_x);
    % gradLEFM.vec.fracE_cohesive = H_b * If2c * gradLEFM.fun.fracE_cohesive(fault_x);
    gradLEFM.vec.fracE_cohesive = ParamFun.Rc(fault_x) .* gradLEFM.fun.fracE_cohesive(fault_x);
    gradLEFM.vec.W_cohesive = ParamFun.Rc(fault_x) .* gradLEFM.fun.W_cohesive(fault_x);

    % fault_xq = fault_x(fault_x > -1 & fault_x < 9);
    % gradLEFM.vec.fracE = zeros(size(fault_xq));
    % for i = 1:length(fault_xq)-1
    %     x1 = fault_xq(i);
    %     x2 = fault_xq(i+1);
    %     gradLEFM.vec.fracE(i) = (gradLEFM.fun.fracE(x1) + gradLEFM.fun.fracE(x2)) * (x2 - x1) / 2;
    % end

    %% Calculate Sensitivity from Adjoint result
    grads = struct();
    grads.tau0 = cell(length(vars.gradients), 1);
    grads.taup = cell(length(vars.gradients), 1);
    grads.taur = cell(length(vars.gradients), 1);
    grads.D_c  = cell(length(vars.gradients), 1);
    grads.D_c_x = cell(length(vars.gradients), 1); % D_c without the direct gradient of D_c

    grads.taup_ = cell(length(vars.gradients), 1);
    grads.D_c_  = cell(length(vars.gradients), 1);
    grads.FracE = cell(length(vars.gradients), 1);
    grads.W = cell(length(vars.gradients), 1);
    grads.FracE_ = cell(length(vars.gradients), 1);
    grads.W_ = cell(length(vars.gradients), 1);


    vTdc = vars.vTdc;
    

    % for i = 1:length(vars.gradients)
    for i = varPlotIter
        grad = vars.gradients{i};
        % convert Dc m to km
        grad.D_c = grad.D_c * 1e3;

        direct_grad_D_c = H_b * dirac_deltas{i} ./ (vTdc * 1e-3);
        for j = 1:length(grads.D_c{i})
            if isnan(direct_grad_D_c(j))
                direct_grad_D_c(j) = 0;
            end
        end
        % if vars have the filed 'RupDefslip', 
        if isfield(vars, 'RupDefSlip')
            RupDefSlip = vars.RupDefSlip;
            direct_grad_D_c = RupDefSlip * direct_grad_D_c;
        end
        if ~plotOpts.RejectDeltaSensitivity
            grad.D_c = grad.D_c + direct_grad_D_c;
        end
        
        grads.tau0{i} = grad.tau0;
        grads.taup{i} = grad.taup;
        grads.taur{i} = grad.taur;
        grads.D_c{i} = grad.D_c;
        grads.D_c_x{i} = grad.D_c + direct_grad_D_c;

        grads.FracE{i} = grad.taup ./ ParamVec.D_c + grad.D_c ./ (ParamVec.taup - ParamVec.taur);
        grads.W{i} = grad.taup .* ParamVec.D_c / 2 - grad.D_c .* ParamVec.D_c.^2 ./ (ParamVec.taup - ParamVec.taur) / 2;

        grads.FracE_{i} = (ParamVec.taup - ParamVec.taur) / 2 .* grad.taup + ParamVec.D_c / 2 .* grad.D_c;
        grads.W_{i}     = (ParamVec.taup - ParamVec.taur) / 2 .* grad.taup - ParamVec.D_c / 2 .* grad.D_c;
        
        grads.taup_{i}  = (ParamVec.taup - ParamVec.taur) .* grad.taup;
        grads.D_c_{i}   = ParamVec.D_c .* grad.D_c;
    end


    %% Plot gradients

    source = vars.source;
    xplotLim = plotOpts.xlim;
    yplotLim = plotOpts.ylim;
    ylogplotLim = plotOpts.yloglim;
    yndplotLim = plotOpts.yndlim;
    % xplotLim = [-2, 9];
    % yplotLim = [0, 0.5];
    % ylogplotLim = [1e-3, 10];
    % yplotLim2 = [1e-7, 1e-4];

    % Plot grad FE with W and taur fixed 
    f = figure();
    hold on;
    % for i = 1:length(vars.gradients)
    for i = varPlotIter
        plot(fault_x, grads.FracE{i}*1e-3, 'linewidth', 2, 'DisplayName', sprintf('x_m = %.2f', full(source{i}.x{1})));
    end
    plot(fault_x, gradLEFM.vec.fracE*1e-3, 'k--', 'linewidth', 2, 'DisplayName', 'LEFM');
    legend();
    xlabel('x along fault (km)');
    title('$\delta T_{rup} (x_m) / \delta \Gamma(x) \; [\mathrm{s/MPa} \cdot \mathrm{m}]$','Interpreter','latex')
    % title('$\mathrm{Grad:} (\tau_p - \tau_r)D_c / 2 \mathrm{(MPa \cdot m)}$','Interpreter','latex')
    grid on;
    xlim(xplotLim);
    ylim(yplotLim);
    pdf_output(f, sprintf('%s/grad_dT_dFE_LEFM', savedir));
    ylim(ylogplotLim);
    yscale('log');
    grid on;
    set(gca,'fontsize',fontsize);
    pdf_output(f, sprintf('%s/grad_dT_dFE_LEFM_log', savedir));
    

    %%% Plot grad taup with D_c and taur fixed
    f2 = figure();
    hold on;
    % for i = 1:length(vars.gradients)
    for i = varPlotIter
        % grad_fracE1 = grads.taup{i} * 2 ./ ParamVec.D_c;
        % plot(fault_x, grad_fracE1*1e-3, 'linewidth', 2, 'DisplayName', sprintf('x_m = %.2f', full(source{i}.x{1})));
        plot(fault_x, grads.taup_{i}, 'linewidth', 2, 'DisplayName', sprintf('x_m = %.2f', full(source{i}.x{1})));
    end
    plot(fault_x, gradLEFM.vec.fracE_, 'k--', 'linewidth', 2, 'DisplayName', 'LEFM');
    legend();
    xlabel('x along fault (km)');
    title('$\delta T_{rup} (x_m) / \delta \tau_p (x) \times (\tau_p - \tau_r) \; [\mathrm{s/MPa} \cdot \mathrm{m}]$','Interpreter','latex')
    % title('$\mathrm{Grad:} \tau_p - \tau_r \mathrm{(MPa)}$','Interpreter','latex')
    grid on;
    set(gca,'fontsize',fontsize);
    xlim(xplotLim);
    ylim(yplotLim);
    pdf_output(f2, sprintf('%s/grad_taup', savedir));
    ylim(ylogplotLim);
    yscale('log');
    pdf_output(f2, sprintf('%s/grad_taup_log', savedir));

    %%% Plot grad D_c with taup and taur fixed
    f3 = figure();
    hold on;
    % for i = 1:length(vars.gradients)
    for i = varPlotIter
        % grad_fracE2 = grads.D_c{i} * 2 ./ (ParamVec.taup - ParamVec.taur);
        % plot(fault_x, grad_fracE2*1e-3, 'linewidth', 2, 'DisplayName', sprintf('x_m = %.2f', full(source{i}.x{1})));
        plot(fault_x, grads.D_c_{i}, 'linewidth', 2, 'DisplayName', sprintf('x_m = %.2f', full(source{i}.x{1})));
    end
    plot(fault_x, gradLEFM.vec.fracE_, 'k--', 'linewidth', 2, 'DisplayName', 'LEFM');
    legend();
    xlabel('x along fault (km)');
    title('$\delta T_{rup} (x_m) / \delta D_c (x) \times D_c \; [\mathrm{s/MPa} \cdot \mathrm{m}]$','Interpreter','latex')
    % title('$\mathrm{Grad:} \tau_p - \tau_r \mathrm{(MPa)}$','Interpreter','latex')
    % grid on;
    xlim(xplotLim);
    ylim(yplotLim);
    grid on;
    set(gca,'fontsize',fontsize);
    pdf_output(f3, sprintf('%s/grad_Dc', savedir));
    ylim(ylogplotLim);
    yscale('log');
    pdf_output(f3, sprintf('%s/grad_Dc_log', savedir));

    %%% Plot grad non-dimensionalised FE with W and taur fixed
    f4 = figure();
    hold on;
    % for i = 1:length(vars.gradients)
    for i = varPlotIter
        plot(fault_x, grads.FracE_{i}, 'linewidth', 2, 'DisplayName', sprintf('x_m = %.2f', full(source{i}.x{1})));
    end
    plot(fault_x, gradLEFM.vec.fracE_, 'k--', 'linewidth', 2, 'DisplayName', 'LEFM');
    if plotOpts.CohesiveAnalysis
        plot(fault_x, gradLEFM.vec.fracE_ + gradLEFM.vec.fracE_cohesive, 'k--', 'linewidth', 2, 'DisplayName', 'coh.');
    end
    xlim(xplotLim);
    % yscale('log');
    ylim(yndplotLim);
    grid on;
    legend();
    xlabel('x along fault (km)');
    title('$\delta T_{rup} (x_m) / \delta \Gamma(x) \times \Gamma (x) \; [\mathrm{s}]$','Interpreter','latex')
    set(gca,'fontsize',fontsize);
    pdf_output(f4, sprintf('%s/grad_dT_dFE_nonDim', savedir));
    ylim(ylogplotLim);
    yscale('log');
    pdf_output(f4, sprintf('%s/grad_dT_dFE_nonDim_log', savedir));

    if plotOpts.CohesiveAnalysis
        f6 = figure();
        hold on;
        for i = varPlotIter
            plot(fault_x, grads.FracE_{i} - gradLEFM.vec.fracE_, 'linewidth', 2, 'DisplayName', sprintf('x_m = %.2f', full(source{i}.x{1})));
        end
        plot(fault_x, gradLEFM.vec.fracE_cohesive, 'k--', 'linewidth', 2, 'DisplayName', 'coh.');
        xlim(xplotLim);
        % ylim([-0.1, 0.1]);
        ylim([-1e-2, max(gradLEFM.vec.fracE_cohesive)*1.2]);
        grid on;
        legend();
        xlabel('x along fault (km)');
        title('$\delta T_{rup} (x_m) / \delta \Gamma(x) \times \Gamma (x) - \mathrm{LEFM} \; [\mathrm{s}]$','Interpreter','latex')
        set(gca,'fontsize',fontsize);
        pdf_output(f6, sprintf('%s/grad_dT_dFE_nonDim_coh', savedir));
        yscale('log');
        ylim([1e-4, max(gradLEFM.vec.fracE_cohesive)*1.2]);
        pdf_output(f6, sprintf('%s/grad_dT_dFE_nonDim_coh_log', savedir));
    end


    %%% Plot grad non-dimensionalised W with FE and taur fixed
    f5 = figure();
    hold on;
    % for i = 1:length(vars.gradients)
    for i = varPlotIter
        plot(fault_x, abs(grads.W_{i}), 'linewidth', 2, 'DisplayName', sprintf('x_m = %.2f', full(source{i}.x{1})));
    end
    % plot(fault_x, gradLEFM.vec.fracE_, 'k--', 'linewidth', 2, 'DisplayName', 'LEFM');
    if plotOpts.CohesiveAnalysis
        plot(fault_x, abs(gradLEFM.vec.W_cohesive), 'k--', 'linewidth', 2, 'DisplayName', 'coh');
    end
    xlim(xplotLim);
    grid on;
    ylim(yndplotLim);
    legend();
    % yscale('log');
    xlabel('x along fault (km)');
    title('$|\delta T_{rup} (x_m) / \delta W(x)| \times W(x) \; [\mathrm{s}]$','Interpreter','latex')
    set(gca,'fontsize',fontsize);
    pdf_output(f5, sprintf('%s/grad_dT_dW_nonDim', savedir));
    ylim(ylogplotLim);
    yscale('log');
    pdf_output(f5, sprintf('%s/grad_dT_dW_nonDim_log', savedir));

    if ~ isempty(plotOpts.varPlotOneIter)
        i = plotOpts.varPlotOneIter; 
        f7 = figure();
        hold on;

        plot(fault_x, abs(grads.FracE_{i}), 'linewidth', 3, 'DisplayName', 'FE');
        plot(fault_x, abs(grads.W_{i}), 'linewidth', 3, 'DisplayName', 'W', 'linestyle', ':');
        plot(fault_x, abs(gradLEFM.vec.fracE_) .* (fault_x <= source{i}.x{1}), 'k--', 'linewidth', 3, 'DisplayName', 'FE by LEFM');
        legend();
        xlabel('x along fault (km)');
        title('$|\delta T_{rup} (x_m) / \delta \bar{p}(x)| \; [\mathrm{s}]$','Interpreter','latex')
    % title('$\mathrm{Grad:} (\tau_p - \tau_r)D_c / 2 \mathrm{(MPa \cdot m)}$','Interpreter','latex')
        grid on;
        xlim(xplotLim);
        ylim(yplotLim);
        set(gca,'fontsize',fontsize);
        pdf_output(f7, sprintf('%s/grad_dT_dp_xm=%.1f', savedir, full(source{i}.x{1})));
    end



%%% Botu
% approx_threshold_x = 2;
% approx_threshold_iter = find(fault_x(rupIterNum) > approx_threshold_x, 1);
% rupIterNum1 = rupIterNum(1:approx_threshold_iter);
% rupIterNum2 = rupIterNum(approx_threshold_iter+1:end);

% smoothparam1 = 0.8;
% smoothparam2 = 0.3;
% fitobj1 = fit(fault_x(rupIterNum1), Trup_Num(rupIterNum1), 'smoothingspline', 'SmoothingParam', smoothparam1);
% fitobj2 = fit(fault_x(rupIterNum2), Trup_Num(rupIterNum2), 'smoothingspline', 'SmoothingParam', smoothparam2);

% dT_dx1 = differentiate(fitobj1, fault_x(rupIterNum1));
% dT_dx2 = differentiate(fitobj2, fault_x(rupIterNum2));

% Trup_num_approx = [feval(fitobj1, fault_x(rupIterNum1)); feval(fitobj2, fault_x(rupIterNum2))];
% Vr_Num = [1 ./ dT_dx1; 1 ./ dT_dx2];


% approximation by quadrature
% vr_num = a3 x^2 + a2 x + a1
% Ndim = 4;
% Kernel = zeros(length(rupIterNum), Ndim);
% Kernel_div = zeros(length(rupIterNum), Ndim);
% for i = 1:length(rupIterNum)
%     idx = rupIterNum(i);
%     for j = 1:Ndim
%         Kernel(i, j) = fault_x(idx)^(j-1);
%         Kernel_div(i, j) = (j-1) * fault_x(idx)^(j-2);
%     end
% end
% a = lsqr(Kernel, Trup_Num(rupIterNum));

% Trup_num_approx = Kernel * a;
% Vr_Num = 1 ./ (Kernel_div * a);

% Vr_Num = NaN * ones(length(rupIterNum), 1);
% for ir = 1:length(rupIterNum)
%     i = rupIterNum(ir);
%     Vr_Num(ir) = (fault_x(i+4) - fault_x(i-4)) / (Trup_Num(i+4) - Trup_Num(i-4));
% end
% Vr_num = min(Vr_Num, cs);

