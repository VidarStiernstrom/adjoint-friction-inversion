function plot_obs(obs, saveOpts, plotReceiver, clims, plotOpts)
    default_arg('plotReceiver', false);
    default_arg('clims', []);
    default_arg('plotOpts', struct());
    default_field(plotOpts, 'timeseries', false);
    default_field(plotOpts, 'plot_receivers', []);
    default_field(plotOpts, 'delv', 0.3);
    default_field(plotOpts, 'truePlot', false);
    default_field(plotOpts, 'fontsize', 20);

    fontsize = plotOpts.fontsize;

    Nobs = numel(obs.loc);
    obs_x = zeros(1, Nobs);
    obs_y = zeros(1, Nobs);
    obs_r = zeros(1, Nobs);
    for i = 1:Nobs
        obs_x(i) = obs.loc{i}(1);
        obs_y(i) = obs.loc{i}(2);
        obs_r(i) = sqrt(obs_x(i)^2 + obs_y(i)^2) * ((obs_y(i) >= 0)*2 - 1);
    end
    [~, sorted_iters] = sort(obs_r); % Sort indices based on obs_r
    data = obs.data(sorted_iters,:)'; % Reorder data based on obs_r

    if plotOpts.truePlot
        trueObs = plotOpts.trueObs;
        true_data = trueObs.data(sorted_iters,:)';
    end

    % [R, T] = meshgrid(obs_r, obs.t);
    [Iter, T] = meshgrid(1:88, obs.t); % 2001 x 88

    close all;
    f = figure();
    subplot(1,1,1)
    s = pcolor(Iter, T, data);
    set(s, 'EdgeColor', 'none'); % Remove grid lines
    if ~isempty(clims)
        caxis(clims); % Set color axis limits
    end
    elastic.cmap;
    colorbar;
    xlabel('obs iteration');
    ylabel('t');
    grid on;
    set(gca, 'fontsize', fontsize);
    pdf_output(f, sprintf('%s%s', saveOpts.savepath, saveOpts.filename));

    % velocity time series
    % plot_receivers = [46, 48, 50, 53, 54, 60, 62, 66, 74, 80];
    if plotOpts.timeseries
        plot_receivers = plotOpts.plot_receivers;
        delv = plotOpts.delv;
        f = figure();
        subplot(1,1,1)
        hold on;
        for i = 1:length(plot_receivers)
            idx = plot_receivers(i);
            if plotOpts.truePlot
                plot(trueObs.t, true_data(:,idx) + i*delv, 'b', 'linewidth', 3);
                plot(obs.t, data(:,idx) + i*delv, 'r--', 'linewidth', 2);
            else
                plot(obs.t, data(:,idx) + i*delv, 'linewidth', 3, 'DisplayName', sprintf('Index %d', idx));
            end
        end
        xlabel('t');
        grid on;
        if ~plotOpts.truePlot
            legend('show', 'Location', 'best');
        end
        xlim([0, obs.t(end)]);
        set(gca, 'fontsize', fontsize);
        pdf_output(f, sprintf('%s%s_tV', saveOpts.savepath, saveOpts.filename));
    end


    if plotReceiver
        f2 = figure();
        subplot(1,1,1)
        hold on;
        plot([-10,10], [0,0], 'linewidth', 2);
        scatter(obs_x, obs_y, 50, 'k', 'x');
        axis equal;
        xlabel('x (km)');
        ylabel('y (km)');
        xlim([-10, 10]);
        ylim([-10, 10]);
        box on;
        ax = gca;
        ax.LineWidth = 1.5;
        set(gca, 'fontsize', fontsize);
        pdf_output(f2, sprintf('%sreceiver', saveOpts.savepath));

        f3 = figure();
        subplot(1,1,1)
        hold on;
        plot([-10,10], [0,0], 'linewidth', 2);
        textscatter(obs_x(sorted_iters), obs_y(sorted_iters), string(1:Nobs)); %, 'FontSize', 12, 'Color', 'k', 'Marker', 'x');
        axis equal;
        xlabel('x (km)');
        ylabel('y (km)');
        xlim([-10, 10]);
        ylim([-10, 10]);
        box on;
        ax = gca;
        ax.LineWidth = 1.5;
        set(gca, 'fontsize', fontsize);
        pdf_output(f3, sprintf('%sreceiver_iter', saveOpts.savepath));
    end
end