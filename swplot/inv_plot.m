% clearvars;
% close all;
% addSubpaths;
% addpath('./swadjoint');

function inv_plot(loadOpts, plot_iter, lim, figsize, fname, finalplot, fontsize, plot_var, plot_marker_type)
    default_arg('plot_iter', []); % iterations to plot
    default_arg('lim', []); % y-axis limits
    default_arg('figsize', []); % figure size
    default_arg('fname', 'result');
    default_arg('finalplot', false);
    default_arg('fontsize', 16);
    default_arg('plot_var', 'tau0');
    default_arg('plot_marker_type', {});

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

    switch plot_var
        case 'tau0'
            if finalplot
                plot(x, truePars.tau0, x, initialPars.tau0, ':x', x, tau0_inverted, 'k--x', 'linewidth', 2);
                legend('$\mathbf{\tau_0}$, true','$\mathbf{\tau_0}$, initial','$\mathbf{\tau_0}$, final','interpreter','latex','Location','southwest');
            else
                plot(x, truePars.tau0, x, initialPars.tau0, ':x', 'linewidth', 2);
                legend('$\mathbf{\tau_0}$, true','$\mathbf{\tau_0}$, initial','interpreter','latex','Location','southwest');
            end
            hold on;
            for i = plot_iter
                plot(x,history.x(:,i+2),':x','linewidth',2, 'DisplayName', sprintf('$\\mathbf{\\tau_0}$, iter %d', i));
            end
            title('Inverted $\mathbf{\tau_0}$ [MPa]', 'interpreter','latex');

        case 'deltau'
            truePars.deltau = truePars.tau0 - truePars.taur;
            deltau_inverted = tau0_inverted - truePars.taur;
            initialPars.deltau = initialPars.tau0 - truePars.taur;

            if finalplot
                plot(x, truePars.deltau, x, initialPars.deltau, ':x', x, deltau_inverted, 'k--x', 'linewidth', 2);
                legend('true','initial','final','Location','southwest');
            else
                plot(x, truePars.deltau, x, initialPars.deltau, ':x', 'linewidth', 2);
                legend('true','initial','Location','southwest');
            end
            hold on;
            for i = plot_iter
                deltau_iter = history.x(:,i+2) - truePars.taur;
                plot(x, deltau_iter, ':x', 'linewidth', 2, 'DisplayName', sprintf('Iter %d', i));
            end
            title('Inverted $\mathbf{\Delta \tau} [MPa]$', 'interpreter','latex');

        case 'deltau_diff'
            truePars.deltau = truePars.tau0 - truePars.taur;
            deltau_inverted = tau0_inverted - truePars.taur;
            initialPars.deltau = initialPars.tau0 - truePars.taur;

            deltau_inverted_diff = deltau_inverted - truePars.deltau;
            initialPars.deltau_diff = initialPars.deltau - truePars.deltau;

            if finalplot
                plot(x, zeros(size(x)), x, initialPars.deltau_diff, ':x', x, deltau_inverted_diff, 'k--x', 'linewidth', 2);
                legend('','initial','final','Location','southwest');
            else
                plot(x, zeros(size(x)), x, initialPars.deltau_diff, ':x', 'linewidth', 2);
                legend('','initial','Location','southwest');
            end
            hold on;

            kount = 1;
            for i = plot_iter
                deltau_iter = history.x(:,i+2) - truePars.taur;
                deltau_iter_diff = deltau_iter - truePars.deltau;
                if isempty(plot_marker_type)
                    plot(x, deltau_iter_diff, ':x', 'linewidth', 2, 'DisplayName', sprintf('Iter %d', i));
                else
                    plot(x, deltau_iter_diff, plot_marker_type{kount}, 'linewidth', 2, 'DisplayName', sprintf('Iter %d', i));
                end
                kount = kount + 1;
            end
            title('Dynamic stress drop: Est. - Tru. [MPa]');
    end
    lgd = legend('show');
    lgd.Location = 'southwest';
    lgd.NumColumns = 2;
    if ~isempty(lim)
        ylim(lim);
    end
    xlabel('$x$ (km)','interpreter','latex');
    grid on;
    set(gca,'fontsize',fontsize);
    pdf_output(f1, sprintf('%s%s', path, fname));
end
