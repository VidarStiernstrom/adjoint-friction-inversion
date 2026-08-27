% clearvars;
% close all;
% addSubpaths;
% addpath('./swadjoint');

function inv_plot(loadOpts, plot_iter, lim, figsize, fname, finalplot)
    default_arg('plot_iter', []); % iterations to plot
    default_arg('lim', []); % y-axis limits
    default_arg('figsize', []); % figure size
    default_arg('fname', 'result');
    default_arg('finalplot', false);

    load(sprintf('%s%s', loadOpts.savepath, loadOpts.filename), 'history', 'm_p', 'initialPars', 'truePars');
    % load(sprintf('%s%s', loadOpts.savepath, loadOpts.filename), 'optObj', 'parSetFun', 'optPars', 'parSetOpts', 'order', 'history', 'm_p', 'initialPars', 'truePars');

    path = loadOpts.savepath;
    % [~, fname, ~] = fileparts(loadOpts.filename);

    % optimization results
    if isempty(figsize)
        f1 = figure();
    else
        f1 = figure('Position', figsize);
    end
    x = linspace(-10,10,m_p); % Parameter grid
    tau0_inverted = history.x(:,end);

    % [parset, ~] = parSetFun(parSetOpts);
    % k = 6;
    % optObj = optObj(parset, optPars, order, k);
    % x0 = history.x(:, end); % parameter vector at iteration: Vector form
    % parsNative = optObj.parametersVectorToNative(x0);
    % tau0_inverted = parsNative.tau0;
    % if isa(tau0_inverted, 'function_handle')
    %     disp('tau0_inverted is a function handle');
    %     tau0_inverted = tau0_inverted(x, 0);
    % end

    if finalplot
        plot(x, truePars.tau0, x, initialPars.tau0, ':x', x, tau0_inverted, 'k--x', 'linewidth', 2);
        legend('$\mathbf{\tau_0}$, true','$\mathbf{\tau_0}$, initial','$\mathbf{\tau_0}$, final','interpreter','latex','Location','southwest');
    else
        plot(x, truePars.tau0, x, initialPars.tau0, ':x', 'linewidth', 2);
        legend('$\mathbf{\tau_0}$, true','$\mathbf{\tau_0}$, initial','interpreter','latex','Location','southwest');
    end
    % plot(x,truePars.tau0, x,initialPars.tau0,':x', x,tau0_inverted,'k--x', 'linewidth',2);
    % legend('$\mathbf{\tau_0}$, true','$\mathbf{\tau_0}$, initial','$\mathbf{\tau_0}$, final','interpreter','latex','Location','southwest');
    hold on;
    for i = plot_iter
        plot(x,history.x(:,i+2),':x','linewidth',2, 'DisplayName', sprintf('$\\mathbf{\\tau_0}$, iter %d', i));
    end
    lgd = legend('show');
    lgd.Location = 'southwest';
    lgd.NumColumns = 2;
    if ~isempty(lim)
        ylim(lim);
    end
    xlabel('$x$ (km)','interpreter','latex');
    grid on;
    set(gca,'fontsize',16);
    pdf_output(f1, sprintf('%s%s', path, fname));
end
