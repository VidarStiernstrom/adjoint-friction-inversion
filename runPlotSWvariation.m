clearvars;
close all;
addSubpaths;
addpath('./swplot');

name = 'var0'; 
outpath = './swout/';
savedir = sprintf('%s%s', outpath, name);
calx = [-1,9];

plotOpts = struct();
plotOpts.xlim = [0.5, 9];
plotOpts.ylim = [0, 0.5];
plotOpts.yloglim = [1e-3, 10];
plotOpts.yndlim = [0, 0.5];
plotOpts.fontsize = 16;
% plotOpts.varPlotIter = 2:8;
% plotOpts.varPlotOneIter = 6;
plotOpts.RejectDeltaSensitivity = false;
plotOpts.CohesiveAnalysis = false;
var_analysis(savedir, calx, plotOpts);

