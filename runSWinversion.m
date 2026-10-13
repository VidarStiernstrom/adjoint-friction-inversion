% inversion to estimate the initial stress using L-BFGS algorithm
% inversion from gradual nucleation by external loading
% increasing Dc

clearvars;
close all;
addSubpaths;
addpath('./swplot');

maxIter = 50; % max number of optimization iterations
m = 251; % Grid points in x-direction
% % for debugging
% maxIter = 1; 
% m = 51;

m_p = 26; % Grid points on parameter grid
scalingStudy = false; 
misfitType = 'velocity';
optPars = {'tau0'};

saveOpts.savepath = './swout/inv0/';
saveOpts.filename = 'result.mat';

if ~exist(saveOpts.savepath, 'dir')
    mkdir(saveOpts.savepath);
end

x = linspace(-10,10,m_p); % Parameter grid
faultx_lim = [-10, 10];

%%% define swParams for true model
x0 = 0;
tau0 = 72;
taur = @(x,y) 68 + 6 * ((6 < abs(x) & abs(x) < 7) .* (abs(x) - 6) + (abs(x) >= 7));
taup = @(x,y) taur(x,y) + 6;
D_c = @(x,y) (0.2 + 0.2 * (abs(x) >= 1.2).*(abs(x) - 1.2));
sigma0 = 120;

plot_pars({D_c}, {'$D_c$'}, {'--x'}, '$D_c$ (m)', m_p, sprintf('%sDc', saveOpts.savepath), faultx_lim, [50, 50, 400, 300]);
plot_pars({taur, tau0}, {'$\tau_r$', '$\tau_0$'}, {'--x', '--x'}, '$\tau$ (MPa)', m_p, sprintf('%stau', saveOpts.savepath), faultx_lim, [50, 50, 400, 300]);

% nucleation by external stress loading
A = 2.1;
std_dev = 1;
tau_L = @(t, x, y) A*exp(-0.5*(x-x0).^2/std_dev^2); % Gaussian function

% pack parameters
swParams = struct();
swParams.taur = taur;
swParams.taup = taup;
swParams.sigma0 = sigma0;
swParams.D_c = D_c;
swParams.tau0 = tau0;
swParams.tau_L = tau_L;


%%% initial value
val = {72.5};
% tau0_ini = @(x,y) (72 + 2.7 * if_nucleation(x,y)); % initial guess for tau0
% val = {tau0_ini}; % ! initial values should be a function handle

loadData = false; 
horizontalFault = false;
loadOptPath = '';
lbfgs = true;
algorithm = 'lbfgs';

[parsVec, misFit, flag, output, history, initialPars, truePars, parSet] = optimize_antiplaneshear_sw(maxIter, m, m_p, misfitType, optPars, val, scalingStudy, swParams, saveOpts, loadData, horizontalFault, loadOptPath, lbfgs);

figure();
tau0_inverted = history.x(:,end); % Inverted results at final iteration
plot(x,truePars.tau0, x,initialPars.tau0,':x', x,tau0_inverted,'k--x', 'linewidth',2);
xlabel('$x$ (km)','interpreter','latex');
legend('$\mathbf{\tau_0}$, true','$\mathbf{\tau_0}$, initial','$\mathbf{\tau_0}$, current','interpreter','latex','Location','southwest');
grid on;
set(gca,'fontsize',16);


