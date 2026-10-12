
% adjoint simulation to calculate the sensitivity of rupture time w.r.t the slip-weakening parameters
% c.f. variational analysis of LEFM

% nucleation with initial slow-propagating rupture
% taup-taur: 4 MPa Dc = 0.2 m / taup-taur: 5 MPa Dc = 0.16 m / taup-taur: 3.2 MPa Dc = 0.25 m

clearvars;
close all;
clf;
addSubpaths;
addpath('./swplot');

% D_c0 = 0.16;
% D_c0 = 0.2;
D_c0 = 0.25;
name = sprintf('var_Dc%.2f', D_c0);

outpath = './swout/';
savedir = sprintf('%s%s', outpath, name);

optObjFunc = @AntiplaneShear2DSWFrictionOpt;
parSetFunc = @pars.swFriction2DHorizontalFaultInversion;

invPars = {'tau0','taup','taur','D_c'}; % parameter to invert for

parSetOpts.m = 251; % Grid points in x-direction
parSetOpts.m_p = 251; % Grid points on parameter grid
% % for debugging
% parSetOpts.m = 51; 
% parSetOpts.m_p = 51;

parSetOpts.scalingStudy = true; 
parSetOpts.misfitType = 'velocity';
parSetOpts.initialGuessScalings = 1; %[];
parSetOpts.initialGuessValues = [];

%%% fault range
fault_xlims = [-8, 12];
parSetOpts.faultGeometry.xlims = fault_xlims;
axlims = [fault_xlims(1), fault_xlims(2), -10, 10];


order = 8; % error occurs if order = 4

maxIter = 3;
loadDataPath = '';
loadOptPath = '';
lbfgs = true;

fontsize = 18;

T = 9.5;

% Location of the adjoint source (i.e., the location where the rupture time is measured) 
adj_fault_source_locs = {[2], [3], [4], [5], [6], [7], [8], [9], [10.5]};
% adj_fault_source_locs = {[7]}; % for debugging
Nadj = length(adj_fault_source_locs);

%%% Run with parameters, save data and create optimization object
% load the paramter of true model
[parset, ~] = parSetFunc(parSetOpts);

x = linspace(fault_xlims(1),fault_xlims(2),parSetOpts.m_p); % Parameter grid

%%% Parameter settings: semi-infinite crack
% Define parameters for forward simulation

if_zero_stressdrop = @(x,y) (x < 0);
if_positive_stressdrop = @(x,y) (x >= 0) & (x <= 11);
if_negative_stressdrop = @(x,y) (x > 11);
if_nucleation = @(x,y) (x >= 0) & (x <= 1);

FracEnegy0 = 0.4; % MPa m
taup0 = 2 * FracEnegy0 / D_c0 + 72; % MPa
disp(sprintf('name: %s, taup - taur = %.2f MPa, W = %.2f MPa / m', name, taup0 - 72, (taup0 - 72) / D_c0));

tau0 = @(x,y) 72 * if_zero_stressdrop(x,y) + 74*if_positive_stressdrop(x,y) + 69*if_negative_stressdrop(x,y);
taup = @(x,y) 72 * if_zero_stressdrop(x,y) + taup0*(if_positive_stressdrop(x,y) + if_negative_stressdrop(x,y));
taur = @(x,y) 72 + x.* 0;
D_c = D_c0;

sigma0 = 120;


% slowly propagating nucleation with external stress loading
vr = 0.5; % km/s
x0 = -0.2;
xnuc0 = 0;

slope = (tau0(xnuc0,0) - taur(xnuc0,0)) / (xnuc0 - x0);
tau_str_L = @(t, x, y) 72 + (x > x0+vr*t).*(x-x0-vr*t)*slope;
parset.friction.funs.tau_str_L = tau_str_L;


parset.friction.swParams.taur = taur;
parset.friction.swParams.taup = taup;
parset.friction.swParams.sigma0 = sigma0;
parset.friction.swParams.D_c = D_c;
parset.friction.swParams.tau0 = tau0;

% plot parameters
mkdir(savedir);
pars = {taup, taur, tau0};
parlabels = {'$\tau_p$', '$\tau_r$','$\tau_0$'};
parstyles = {'-', '-', '-'};
plot_pars(pars, parlabels, parstyles, '$Stress$ (MPa)', parSetOpts.m, sprintf('%s/pars.pdf', savedir), fault_xlims);

optObj = optObjFunc(parset, invPars, order); %, k);

% optObj.setSyntheticReceiverData(trueReceiverData, trueTimeIntegrationData);
initialPars = optObj.pars;

x0 = optObj.parametersNativeToVector(optObj.pars);

%% Run forward simulation

%%% Calculate the gradient
parsNative = optObj.parametersVectorToNative(x0);
optObj.setParameters(parsNative);
optObj.updateForwardDiscr();

fwdSavedir = sprintf('%s/fwd', savedir);
if ~exist(fwdSavedir, 'dir')
    mkdir(fwdSavedir);
end
saveOpts.saveDir = fwdSavedir;
saveOpts.frameSpacing = 10;


tic
optObj.runForward(true, T, saveOpts);
toc
disp('Forward simulation finish')

% adjoint with residual = 0
ReceiverData = optObj.forwardReceiverRecordings;
TimeIntegrationData = optObj.forwardTimeIntegrationData;
optObj.setSyntheticReceiverData(ReceiverData, TimeIntegrationData);

% optObj.updateAdjointDiscr();
T_fwd = optObj.forwardTimeIntegrationData.T;

res_files = dir(sprintf('%s/*.mat', fwdSavedir));
Nres = numel(res_files);

paths.fwd = fwdSavedir;
uvmaxs.fwd = [1,1.5];
make_movie(paths, 'vVS', uvmaxs, Nres, {'fwd'}, axlims, fault_xlims); %make movie for true simulation
plotFlag = struct();
plotFlag.true = false;
plotFlag.fwd = true;
plotFlag.adj = false;
contFlag = struct();
contFlag.Trup = false;
contFlag.Trups = true;
contFlag.Tdc = true;
if isa(D_c, 'function_handle')
    x = optObj.adjointDiscr.flt_pos(:,1);
    contFlag.Dc = D_c(x);
else
    contFlag.Dc = D_c;
end
close all;
make_fig(paths, 0, Nres, fontsize, {'fwd'}, plotFlag, contFlag, fault_xlims);
[Tdc, vTdc, Trup] = optObj.getRuptureTime();

vars = struct();
vars.gradients = cell(Nadj, 1);
vars.source = cell(Nadj, 1);

fprintf('----------------------------\n');
for i = 1:Nadj
    adj_fault_source_loc = adj_fault_source_locs{i};
    fprintf("Iter %d Adjoint:\n", i)
    adj_name = sprintf('adj_x=%.1f', adj_fault_source_loc(1));
    adjSavedir = sprintf('%s%s/%s/adj/', outpath, name, adj_name);

    if ~exist(adjSavedir, 'dir')
        mkdir(adjSavedir);
    end
    saveOpts.saveDir = adjSavedir;

    optObj.updateAdjointVarDiscr(adj_fault_source_loc);
    % get fault source 
    fault_source = optObj.adjointDiscr.fault_sources;
    fprintf('Adjoint source at x = %.3f\n', full(fault_source.x{1}));

    tic
    optObj.runAdjoint(true, T, saveOpts);
    toc
    disp('Adjoint simulation finish')

    grad = optObj.gradientFormula();

    vars.gradients{i} = grad;
    vars.source{i} = fault_source;
    vars.dirac_deltas{i} = optObj.adjointDiscr.fault_dirac_deltas;
    save(sprintf('%s/grad.mat', adjSavedir), 'x', 'grad');

    % plot
    paths.adj = adjSavedir;
    uvmaxs.adj = [0.8, 0.8];
    make_movie(paths, 'vVS', uvmaxs, Nres, {'adj'}, axlims, fault_xlims);
    adjcmax = 0.8;
    plotFlag = struct();
    plotFlag.true = false;
    plotFlag.fwd = false;
    plotFlag.adj = true;
    make_fig(paths, adjcmax, Nres, fontsize,{'fwd', 'adj'}, plotFlag, contFlag, fault_xlims);

    %%% delete results
    delete([adjSavedir, '/*.mat']);
    rmdir(adjSavedir, 's');
    
    % fprintf("%d iter adjoint has finished!\n", i)
    fprintf('----------------------------\n');
    close all;
end

% delete true results
delete([fwdSavedir, '/*.mat']);
rmdir(fwdSavedir, 's');

vars.H_b = optObj.H_b;
vars.If2c = optObj.If2c;

fault_x = optObj.adjointDiscr.flt_pos(:,1);
swParams = parset.friction.swParams;
vars.swParams = swParams;
vars.fault_x = fault_x;
vars.Tdc = Tdc;
vars.vTdc = vTdc;
vars.Trup = Trup;
save(sprintf('%s/vars.mat', savedir), 'vars');

disp("Finish!")


