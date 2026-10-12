% clearvars;
% close all;
% addSubpaths;
% addpath('./swadjoint');

% loadOpts.path = 'swadjoint/out/inv7/';
% loadOpts.filename = 'result.mat';

% test_iters = 3:4;
% eigenvalue_plot = false;
% no_Hessian_plot = false;

function inv_hessian_analysis(loadOpts, plotOpts)
    default_arg('plotOpts', struct());
    default_field(plotOpts, 'test_iters', []);
    default_field(plotOpts, 'eigenvalue_plot', false);
    default_field(plotOpts, 'no_Hessian_plot', true);
    default_field(plotOpts, 'fontsize', 18);

    test_iters = plotOpts.test_iters;
    eigenvalue_plot = plotOpts.eigenvalue_plot;
    no_Hessian_plot = plotOpts.no_Hessian_plot;
    fontsize = plotOpts.fontsize;

    load(sprintf('%s%s', loadOpts.savepath, loadOpts.filename), 'history', 'm_p', 'initialPars', 'truePars');

    savepath = sprintf('%s%s', loadOpts.savepath, 'hessian_analysis/');
    % create directory if it does not exist
    if ~exist(savepath, 'dir')
        mkdir(savepath);
    end

    % delete fist and last colums of history.x, history.fval, and history.grad since they are the initial and final values, which are not plotted in the optimization process
    history.x = history.x(:, 2:end-1);
    history.fval = history.fval(2:end-1);
    history.grad = history.grad(:, 2:end-1);

    s_list = history.x(:, 2:end) - history.x(:, 1:end-1);
    y_list = history.grad(:, 2:end) - history.grad(:, 1:end-1);

    m_lbfgs= 10;
    tol = 1e-12;

    final_iter = 50;
    n = m_p;
    I = eye(n);

    history.B = cell(1, final_iter+1);
    update_pred_list = zeros(n, final_iter+1);
    update_pred_LBFGS_list = zeros(n, final_iter+1);
    update_actual_list = zeros(n, final_iter+1);
    grad_normalized_list = zeros(n, final_iter+1);

    first_valid_index = true;

    for k = 0:final_iter
        % at iteration k, we know g_k, and s_i for i = k-1, k-2, ..., k-m
        g_k = history.grad(:, k+1);

        alpha_list = zeros(1, final_iter+1); % store alpha_i for LBFGS
        
        %%% for LBFGS
        q = g_k;
        for i = k-1:-1:max(k-m_lbfgs, 0)
            s_i = s_list(:, i+1); % s_i = x_{i+1} - x_i 
            y_i = y_list(:, i+1); % y_i = g_{i+1} - g_i
            sy_i = s_i' * y_i;
            if sy_i > tol
                rho_i = 1 / sy_i;
                alpha_i = rho_i * (s_i' * q);
                alpha_list(i+1) = alpha_i;
                q = q - alpha_i * y_i;
            end
        end
        %%%
        
        gamma = 1.0;

        % Choose gamma from the most recent valid pair
        for i = k-1:-1:max(k-m_lbfgs, 0)
            s_i = s_list(:, i+1); % s_i = x_{i+1} - x_i 
            y_i = y_list(:, i+1); % y_i = g_{i+1} - g_i
            sy_i = s_i' * y_i;
            if sy_i > tol
                if first_valid_index
                    fprintf('First valid pair found at iteration %d: sy = %.6e; use gamma = 1\n', i, sy_i);
                    first_valid_index = false;
                    break;
                else
                    yy_i = y_i' * y_i;
                    gamma = sy_i / yy_i;
                    break;
                end
            end
        end

        B = gamma * I;

        z = B * q; %%% for LBFGS

        for i = max(k-m_lbfgs, 0):k-1
            s_i = s_list(:, i+1);
            y_i = y_list(:, i+1);
            sy_i = s_i' * y_i;
            if sy_i > tol
                rho_i = 1 / sy_i;
                V = (I - rho_i * (s_i * y_i'));
                B = V * B * V' + rho_i * (s_i * s_i');
                %%% for LBFGS
                alpha_i = alpha_list(i+1);
                beta_i = rho_i * (y_i' * z);
                z = z + s_i * (alpha_i - beta_i);
                %%%
            end
        end

        % check update direction
        update_pred = - B * g_k;
        update_pred = update_pred / norm(update_pred);

        update_pred_LBFGS = - z / norm(z); %%% for LBFGS

        if k == final_iter
            update_actual = NaN(n, 1); % not update at iteration 50
        else
        update_actual = s_list(:, k+1);
            update_actual = update_actual / norm(update_actual);
        end

        grad_normalized = - g_k / norm(g_k);

        history.B{k+1} = B;
        update_pred_list(:, k+1) = update_pred;
        update_pred_LBFGS_list(:, k+1) = update_pred_LBFGS;
        update_actual_list(:, k+1) = update_actual;
        grad_normalized_list(:, k+1) = grad_normalized;
    end



    %% plot

    final_x = history.x(:, end);
    % test_iters = 3:4;
    % eigenvalue_plot = false;
    % no_Hessian_plot = false;

    %%% plot hessian at test_iters
    for i = 1:length(test_iters)
        iter = test_iters(i);
        InvHessian = history.B{iter+1}; %B{1} -> iter0
        figure;
        imagesc(InvHessian);
        % logscale
        set(gca, 'ColorScale', 'log');
        set(gca, 'fontsize', fontsize);
        colorbar;
        title(sprintf('inv Hessian at iteration %d', iter));
        pdf_output(gcf, sprintf('%s%s%d', savepath, 'invHessian_iter', iter));
        % saveas(gcf, sprintf('%sHessian_iter%d.png', loadOpts.savepath, iter));
    end

    %%% plot final inv hessian
    final_invhessian = history.B{end};
    f = figure;
    imagesc(final_invhessian);
    % logscale
    set(gca, 'ColorScale', 'log');
    set(gca, 'fontsize', fontsize);
    colorbar;
    title(sprintf('inv Hessian at final iteration'));
    pdf_output(f, sprintf('%sinvHessian_final', savepath));

    %%% plot uncertainty from diagonal of final inv hessian
    cov_ = diag(final_invhessian); % Inverse of the covariance matrix is the Hessian at
    % disp(covinv);
    sigv = history.fval(30);
    disp(strcat('normalize loss by final loss: ', num2str(sigv)));
    cov = cov_ * sigv;

    faultx = linspace(-10, 10, m_p);
    f = figure;
    hold on;
    plot(faultx, truePars.tau0, 'linewidth', 2)
    plot(faultx, initialPars.tau0, ':x', 'linewidth', 2)
    plot(faultx, final_x, 'k--', 'linewidth', 2)
    plot(faultx, final_x + 2 * cov, 'b-', 'linewidth', 1);
    plot(faultx, final_x - 2 * cov, 'b-', 'linewidth', 1);
    xlabel('$x$ (km)', 'interpreter', 'latex');
    legend('$\mathbf{\tau_0}$, true', '$\mathbf{\tau_0}$, initial', '$\mathbf{\tau_0}$, final', 'interpreter', 'latex', 'Location', 'southwest');
    grid on;
    set(gca, 'fontsize', fontsize);
    pdf_output(f, sprintf('%s%s', savepath, 'result_final'));


    %%%% plot update direction (actual vs predicted) at test_iters
    for LBFGS_plot = [true, false]
        for i = 1:length(test_iters)
            test_iter = test_iters(i);
            update_pred = update_pred_list(:, test_iter+1);
            update_pred_LBFGS = update_pred_LBFGS_list(:, test_iter+1);
            update = update_actual_list(:, test_iter+1);
            grad_normalized = grad_normalized_list(:, test_iter+1);

            % plot update_pred vs update
            f = figure;
            hold on;
            if ~ no_Hessian_plot
                plot(faultx, update_pred, 'o--', 'linewidth', 2, 'DisplayName', 'predicted update');
            end
            plot(faultx, update, 'x--', 'linewidth', 2, 'DisplayName', 'actual update');
            plot(faultx, grad_normalized, 'k--', 'linewidth', 2, 'DisplayName', 'gradient');
            if LBFGS_plot
                plot(faultx, update_pred_LBFGS, 's--', 'linewidth', 2, 'DisplayName', 'predicted update (LBFGS)');
            end
            legend('show', 'Location', 'best');
            grid on;
            xlabel('$x$ (km)', 'interpreter', 'latex');
            ylabel('update', 'interpreter', 'latex');
            title(sprintf('Update at iteration %d -> %d', test_iter, test_iter+1));
            set(gca, 'fontsize', fontsize);
            if LBFGS_plot
                pdf_output(f, sprintf('%s%s', savepath, 'update_iter' + string(test_iter) + '_LBFGS'));
            elseif no_Hessian_plot
                pdf_output(f, sprintf('%s%s', savepath, 'update_iter' + string(test_iter) + '_noHessian'));
            else
                pdf_output(f, sprintf('%s%s', savepath, 'update_iter' + string(test_iter)));
            end
        end
    end

    %% SVD　for final inv hessian
    if eigenvalue_plot
        [U, S, V] = svd(final_invhessian);
        sig = diag(S);
        f1 = figure;
        semilogy(sig, 'o-', 'linewidth', 2);
        xlabel('Index', 'interpreter', 'latex');
        ylabel('Singular value', 'interpreter', 'latex');
        grid on;
        set(gca, 'fontsize', fontsize);
        pdf_output(f1, sprintf('%s%s', savepath, 'svd_values'));

        f2 = figure();
        % figsize = [50, 50, 400, 200];
        set(gcf, 'Position', [50, 50, 1000, 600]);
        for i = 1:26
            subplot(4, 7, i);
            iter = length(sig) - (i-1);
            plot(faultx, U(:, iter), 'linewidth', 2);
            title(sprintf('SV:%.2e', sig(iter)));
            set(gca, 'fontsize', fontsize);
            grid on;
        end

        pdf_output(f2, sprintf('%s%s', savepath, 'svd_modes'));
    end
end


