% adjoint simulation from the estimation results 

% clearvars;
% close all;
% clf;
% addSubpaths;
% addpath('./swadjoint');

% outpath = './swadjoint/out/inv4/';
% sim_iters = [20];

function inv_sim(loadOpts, sim_iters, T, fontsize)
    default_arg('T', 6);
    default_arg('fontsize', 20);

    load(sprintf('%s%s', loadOpts.savepath, loadOpts.filename), 'optObj', 'parSetFun', 'optPars', 'parSetOpts', 'order', 'history');
    outpath = loadOpts.savepath;

    optObjFunc = optObj;
    parSetFunc = parSetFun;
    invPars = optPars;

    [parset, trueParset] = parSetFunc(parSetOpts);

    %% True forward simulation

    disp('Generating synthetic data')
    trueOptObj = optObjFunc(trueParset, invPars, order);

    trueSavedir = sprintf('%s/true', outpath);
    if ~exist(trueSavedir, 'dir')
        mkdir(trueSavedir);
    end
    saveOpts.saveDir = trueSavedir;
    saveOpts.frameSpacing = 30;

    tic
    trueOptObj.runForward(true, T, saveOpts);
    toc
    disp('True forward simulation finish')

    trueReceiverData = trueOptObj.forwardReceiverRecordings;
    trueTimeIntegrationData = trueOptObj.forwardTimeIntegrationData;
    k = trueOptObj.k;
    truePars = trueOptObj.pars;

    % for seismograms plot
    trueObs = struct();
    trueObs.loc = trueParset.receivers.x;
    nStages = trueTimeIntegrationData.nStages;
    trueObs.data = trueOptObj.removeStagedData(trueOptObj.forwardReceiverRecordings, nStages);
    trueObs.t = trueOptObj.removeStagedData(trueOptObj.forwardTimeIntegrationData.T, nStages);
    clear trueOptObj;

    %% plot true seismograms
    obsSaveOpts.savepath = outpath;
    obsSaveOpts.filename = 'true_obs';
    % plot_obs(trueObs, obsSaveOpts, true, [-0.35, 0.35]); % fault plot / caxis
    % timeseries plot
    tvOpts = struct();
    tvOpts.timeseries = true;
    tvOpts.plot_receivers = [46, 48, 50, 53, 54, 60, 62, 66, 74, 80];
    tvOpts.delv = 0.4;
    tvOpts.truePlot = false;
    tvOpts.fontsize = fontsize;
    plot_obs(trueObs, obsSaveOpts, true, [-0.35, 0.35], tvOpts); % fault plot / caxis

    res_files = dir(sprintf('%s/*.mat', trueSavedir));
    Nres = numel(res_files);

    paths.true = trueSavedir;
    uvmaxs.true = [1,1.5];
    make_movie(paths, 'vVS', uvmaxs, Nres, {'true'}); %make movie for true simulation

    for iter = sim_iters
        fprintf("%d iter adjoint start!\n", iter)
        iter_name = sprintf('iter%d', iter);
        savedir = sprintf('%s%s/', outpath, iter_name);
        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end

        %%% Set adjoint optObj
        optObj = optObjFunc(parset, invPars, order, k);
        optObj.setSyntheticReceiverData(trueReceiverData, trueTimeIntegrationData);

        x0 = history.x(:, iter+2); % parameter vector at iteration: Vector form
        % x0 = optObj.parametersNativeToVector(optObj.pars);
        parsNative = optObj.parametersVectorToNative(x0);

        optObj.setParameters(parsNative);
        optObj.updateForwardDiscr();

        %%%% Plot parameters
        taur = trueParset.friction.swParams.taur;
        taup = trueParset.friction.swParams.taup;
        tau0 = trueParset.friction.swParams.tau0;
        tau0_fwd = parsNative.tau0;

        x = linspace(-10,10,parSetOpts.m_p); % Parameter grid
        pars = {taup, taur, tau0, tau0_fwd};
        parlabels = {'$\tau_p$', '$\tau_r$','$\tau_0$ true', '$\tau_0$ fwd'};
        parstyles = {'-', '-', '-', 'k--'};
        plot_pars(pars, parlabels, parstyles, '$Stress$ (MPa)', parSetOpts.m_p, sprintf('%s%s/pars.pdf', outpath, iter_name));

        %%% Run Forward 
        fwdSavedir = sprintf('%s%s/fwd', outpath, iter_name);
        if ~exist(fwdSavedir, 'dir')
            mkdir(fwdSavedir);
        end
        saveOpts.saveDir = fwdSavedir;

        tic
        optObj.runForward(true, T, saveOpts);
        toc
        disp('Forward simulation finish')

        % plot fowrward seismograms
        fwdObs = trueObs;
        fwdObs.data = optObj.removeStagedData(optObj.forwardReceiverRecordings, nStages);
        fwdObsSaveOpts.savepath = savedir;
        fwdObsSaveOpts.filename = 'fwd_obs';
        % plot_obs(fwdObs, fwdObsSaveOpts, false, [-0.35, 0.35]);
        tvOpts.truePlot = true;
        tvOpts.trueObs = trueObs;
        tvOpts.timeseries = true;
        plot_obs(fwdObs, fwdObsSaveOpts, false, [-0.35, 0.35], tvOpts);
        resObs = trueObs;
        resObs.data = fwdObs.data - trueObs.data;
        resObsSaveOpts.savepath = savedir;
        resObsSaveOpts.filename = 'res_obs';
        clim = max(abs(resObs.data(:))) + 1e-10;
        % plot_obs(resObs, resObsSaveOpts, false, [-clim, clim]);
        tvOpts.timeseries = false;
        plot_obs(resObs, resObsSaveOpts, false, [-clim, clim], tvOpts);

        %%% Run adjoint
        optObj.updateAdjointDiscr();

        adjSavedir = sprintf('%s%s/adj', outpath, iter_name);
        if ~exist(adjSavedir, 'dir')
            mkdir(adjSavedir);
        end
        saveOpts.saveDir = adjSavedir;

        tic
        optObj.runAdjoint(true, T, saveOpts);
        toc
        disp('Adjoint simulation finish')

        grad = optObj.gradientFormula();

        save(sprintf('%s/grad.mat', adjSavedir), 'x', 'grad');

        % plot
        paths.fwd = fwdSavedir; % paths for result files
        paths.adj = adjSavedir;
        uvmaxs.fwd = [1,1.5];
        % uvmaxs.adj = [0.05, 0.05];
        uvmaxs.adj = [0.1, 0.1];
        make_movie(paths, 'vVS', uvmaxs, Nres, {'fwd', 'adj'});
        adjcmax = 0.1;
        make_fig(paths, adjcmax, Nres, fontsize);

        %%% delete results
        delete([fwdSavedir, '/*.mat']);
        delete([adjSavedir, '/*.mat']);
        rmdir(fwdSavedir, 's');
        rmdir(adjSavedir, 's');
        
        fprintf("%d iter adjoint has finished!\n", iter)
        close all;
    end

    % delete true results
    delete([trueSavedir, '/*.mat']);
    rmdir(trueSavedir, 's');
end
