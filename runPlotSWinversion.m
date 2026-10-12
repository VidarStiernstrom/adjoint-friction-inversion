clearvars;
close all;
addSubpaths;
addpath('./swplot');

%% Parameter setup
loadOpts.savepath = './swout/inv0/';
loadOpts.filename = 'result.mat';
fontsize = 18;

%%% Figure 2 / 3 / 4

%%% Simulation with inverted parameters
T = 6;
sim_iters = [3,4,50];
inv_sim(loadOpts, sim_iters, T, fontsize);


%%% plot learning curve
plot_loss(loadOpts, fontsize);


%%% Inverted parameter plot
plot_iters = [2, 5, 15, 20, 40];
lim = [];
figsize = [50, 50, 600, 400];
finalplot = true; % if true, plot the final result of the inversion
plot_var = 'deltau_diff'; % 'tau0' or 'deltau'
inv_plot(loadOpts, plot_iters , lim, figsize, 'result', finalplot, fontsize, plot_var);
inv_plot(loadOpts, plot_iters , lim, figsize, 'result_deltau', finalplot, fontsize, 'deltau');

plot_marker_type = {':square', ':x', ':x'};
inv_plot(loadOpts, [3,4,5] , lim, figsize, 'result_iter345', false, fontsize, plot_var, plot_marker_type);

inv_plot(loadOpts, [] , [-0.3, 0.7], figsize, 'result_final', finalplot, fontsize, plot_var);

%%% Hessian analysis
plotOpts = struct();
plotOpts.test_iters = [3,4];
plotOpts.eigenvalue_plot = false;
plotOpts.no_Hessian_plot = true;
plotOpts.fontsize = fontsize;
inv_hessian_analysis(loadOpts, plotOpts);

