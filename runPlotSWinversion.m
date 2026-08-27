clearvars;
close all;
addSubpaths;
addpath('./swplot');

loadOpts.savepath = './swout/inv0/';
loadOpts.filename = 'result.mat';
plot_iters = 1:10:50;
sim_iters = [1,4,50];
% % for debugging
% plot_iters = 1; 
% sim_iters = 1;
T = 6;
lim = [71.5, 72.7];
figsize = [50, 50, 600, 400];
finalplot = true; % if true, plot the final result of the inversion
inv_plot(loadOpts, plot_iters, lim, figsize, 'result', finalplot);
inv_sim(loadOpts, sim_iters, T);

