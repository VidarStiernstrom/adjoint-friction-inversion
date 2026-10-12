clearvars;
close all;
addSubpaths;
% addpath('./swplot');
addpath('./swadjoint_paper/swplot');

% name = 'var0'; 
outpath = './swadjoint_paper/swout/';

for Dc = [0.16, 0.2, 0.25]
    name = sprintf('var_Dc%.2f', Dc);

    savedir = sprintf('%s%s', outpath, name);
    % calx = [-1,9];
    calx = [1, 11];

    plotOpts = struct();
    % plotOpts.xlim = [0.5, 9];
    plotOpts.xlim = [2.5, 11];
    plotOpts.ylim = [0, 0.5];
    plotOpts.yloglim = [1e-3, 10];
    plotOpts.yndlim = [0, 0.5];
    plotOpts.fontsize = 16;
    plotOpts.varPlotOneIter = 1;
    if Dc == 0.2
        plotOpts.varPlotIter = 2:8;
        plotOpts.varPlotOneIter = 6;
    end
    plotOpts.RejectDeltaSensitivity = false;
    plotOpts.CohesiveAnalysis = false;
    var_analysis(savedir, calx, plotOpts);
end

