function [parsVec, misFit, flag, output, history, initialPars, truePars, parSet] = optimize_antiplaneshear_sw(maxIter, m, m_p, misfitType, optPars, val, scalingStudy, swParams, saveOpts, loadData, horizontalFault, loadOptPath, lbfgs, algorithm)
    
    addSubpaths;
    default_arg('swParams',[]);
    default_arg('scalingStudy',true);
    default_arg('loadData',false);
    default_arg('horizontalFault',false);
    default_arg('loadOptPath','');
    default_arg('lbfgs',false);
    % added by R.F.
    default_arg('algorithm', 'bfgs');
    %
    default_arg('saveOpts',struct);
    % default_arg('filename','swadjoint/inversion/test.mat');
    order = 8;
    % optObj = @AntiplaneShear2DRSFrictionOpt;
    optObj = @AntiplaneShear2DSWFrictionOpt;
    
    % if loadData
    %      % Generate this data by running fractal_fault_generate_data.m with the 
    %      % pars.rsFriction2DFractalFaultForward parameter set.
    %     assert(isequal(misfitType,'velocity'));
    %     data_str = 'highres';
    %     loadDataPath = 'mat/fractal_fault_m1001_mp1001/receiverData.mat';
    % else
    %     loadDataPath = '';
    %     data_str = 'synth';
    % end
    loadDataPath = '';
    data_str = 'synth';

    % if horizontalFault
    %     fault_str = 'horiz';
    %     parSetFun = @pars.rsFriction2DHorizontalFaultInversion;
    % else
    %     fault_str = 'rough';
    %     parSetFun = @pars.rsFriction2DFractalFaultInversion;
    % end
    fault_str = 'horiz';
    parSetFun = @pars.swFriction2DHorizontalFaultInversion;

    if scalingStudy
        parSetOpts.initialGuessScalings = val;
        parSetOpts.initialGuessValues = [];
    else
        parSetOpts.initialGuessScalings = [];
        parSetOpts.initialGuessValues = val;
    end
    parSetOpts.misfitType = misfitType;
    parSetOpts.inversionParameters = optPars;
    parSetOpts.m = m;
    parSetOpts.m_p = m_p;
    %%% added by R.F.
    parSetOpts.swParams = swParams;
    if lbfgs
        algorithm = 'lbfgs';
        disp('Running L-BFGS optimization');
    else
        disp(['Running ', algorithm, ' optimization']);
    end
    %%%
    % [parsVec, misFit, flag, output, history, initialPars, truePars, parSet] = adjopt.optimize_lbfgs(optObj, parSetFun, optPars, parSetOpts, order, maxIter, loadDataPath, loadOptPath, lbfgs);
    [parsVec, misFit, flag, output, history, initialPars, truePars, parSet, hessian] = adjopt.optimize_func(optObj, parSetFun, optPars, parSetOpts, order, maxIter, loadDataPath, loadOptPath, algorithm);


    %%% deleted by R.F.
    % if lbfgs
    %     opt_str = 'lbfgs';
    % else
    %     opt_str = 'bfgs';
    % end

    % if scalingStudy
    %     scalings = parSetOpts.initialGuessScalings;
    %     str_ending = 'scaled';
    %     for i = 1:length(scalings)
    %         parname = optPars{i};
    %         s = scalings(i);
    %         str_ending = erase(sprintf('%s_%s%1.2f',str_ending,parname,s),'.');
    %     end
    % else
    %     vals = parSetOpts.initialGuessValues;
    %     str_ending = 'vals';
    %     for i = 1:length(vals)
    %         parname = optPars{i};
    %         p = vals(i);
    %         % str_ending = erase(sprintf('%s_%s%1.4f',str_ending,parname,p),'.');
    %         str_ending = erase(sprintf('%s_%s%1.1f',str_ending,parname,p),'.');
    %     end
    % end
    % filename = sprintf('swadjoint/inversion/%s_%s_maxiter%d_m%d_mp%d_%s_%s_%s.mat',fault_str,opt_str,maxIter,m,m_p,misfitType,data_str,str_ending);
    %%%

    if isfield(saveOpts, 'filename')
        filename = saveOpts.filename;
    else
        if horizontalFault
            fault_str = 'horiz';
        else
            fault_str = 'rough';
        end
        if lbfgs
            opt_str = 'lbfgs';
        else
            % opt_str = 'bfgs';
            opt_str = algorithm;
        end
        filename = sprintf('%s_%s_maxiter%d_m%d_mp%d_%s_%s.mat', fault_str, opt_str, maxIter, m, m_p, misfitType, data_str);
    end

    if isfield(saveOpts, 'savepath')
        savepath = saveOpts.savepath;
    else
        savepath = 'swadjoint/inversion/';
    end

    if isfile(filename)
        filecounter = 1;
        split_str = split(filename,'.');
        newfilename = [split_str{1},'_',num2str(filecounter),'.',split_str{2}];
        while isfile(newfilename)
            filecounter = filecounter+1;
            newfilename = [split_str{1},'_',num2str(filecounter),'.',split_str{2}];
        end
        filename = newfilename;
    end
    % mkdir('swadjoint/inversion/');
    % save(filename, '-v7.3','-nocompression');
    if ~exist(savepath, 'dir')
        mkdir(savepath);
    end
    save(sprintf('%s%s', savepath, filename), '-v7.3', '-nocompression');
end