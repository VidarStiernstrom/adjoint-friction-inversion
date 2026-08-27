function make_fig(paths, adjcmax, Nres, sim_types, plotFlag, contFlag)
    default_arg('Nres', 69); % number of results to plot
    default_arg('sim_types', {'true', 'fwd', 'adj'});
    default_arg('plotFlag', struct());
    default_field(plotFlag, 'true', true);
    default_field(plotFlag, 'fwd', true);
    default_field(plotFlag, 'adj', true);
    default_arg('contFlag', struct());
    default_field(contFlag, 'Trup', false);
    default_field(contFlag, 'Trups', true);
    default_field(contFlag, 'Tdc', true);
    default_field(contFlag, 'Dc', 0.2);

    % sim_types = {'true','fwd', 'adj'};
    s_true = [];
    V_true = [];

    s_fwd = [];
    V_fwd = [];

    taustr_adj = [];
    V_adj = [];

    t_true = [];
    t = [];
    t_adj = [];

    % for i = 1:3
    for i = 1:length(sim_types)
        sim_type = sim_types{i};
        % savedir =  sprintf('%s%s/%s', outpath, name, sim_type);
        savedir = paths.(sim_type);
        for it = 1:Nres
            fpath = sprintf('%s/%d.mat',savedir,it);
            load(fpath);
            r = saveData;
            if it == 1
                faultx = full(r.fault(:,1));
            end
            if sim_type == "true"
                s_true = [[s_true], [r.s]];
                V_true = [[V_true], [r.V]];
                t_true = [[t_true], [r.t]];
            end
            if sim_type == "fwd"
                s_fwd = [[s_fwd], [r.s]];
                V_fwd = [[V_fwd], [r.V]];
                t = [[t], [r.t]];
            elseif sim_type == "adj"
                taustr_adj = [[taustr_adj], [r.s]];
                V_adj = [[V_adj], [r.V]];
                t_adj = [[t_adj], [r.t]];
            end
        end
        if sim_type == "adj"
            t_adj = t_adj(1,end) - t_adj;
            [T_adj, X_adj] = meshgrid(t_adj, faultx);
        end
        if sim_type == "fwd"
            Tend = t(1, end);
            [T,X] = meshgrid(t, faultx);
        end
    end

    % t_adj = t_adj(1,end) - t_adj;
    % [T_adj, X_adj] = meshgrid(t_adj, faultx);

    % [T,X] = meshgrid(t, faultx);
    % Tend = t_true(1, end);

    if plotFlag.fwd == true
        % Fault properties
        f1 = figure;
        subplot(1,1,1);
        hold on;
        s1 = pcolor(X, T, s_fwd);
        % s.FaceColor = 'interp';
        s1.EdgeColor = 'none';

        % contour(X, T, s_fwd,  [1e-5, 0.2], 'k', 'linewidth', 1.5);
        if contFlag.Trup == true
            contour(X, T, V_fwd > 0, 'k', 'linewidth', 0.1);
        elseif contFlag.Trups == true
            contour(X, T, s_fwd > 0, 'k', 'linewidth', 1.5);
        end
        if contFlag.Tdc == true
            contour(X, T, s_fwd > contFlag.Dc, 'k', 'linewidth', 1.5);
        end

        xlabel('x along fault');
        ylabel('t');
        title('Slip: forward simulation');
        colorbar;
        % pdf_output(f1, sprintf('%s%s/fwd_s', outpath, name));
        pdf_output(f1, sprintf('%s/../fwd_s', paths.fwd));

        % foreward slip velocity
        f2 = figure;
        subplot(1,1,1);
        hold on;
        s2 = pcolor(X, T, V_fwd);
        s2.EdgeColor = 'none';
        axis([-10, 10, 0, Tend]);
        % contour(X, T, V_fwd, [1e-5, 1], 'k', 'linewidth', 1.5);
        xlabel('x along fault');
        ylabel('t');
        title("V: forward");
        % cmax = max(abs(V_fwd(:))) + 1e-10;
        caxis([-2, 2]);
        elastic.cmap;
        colorbar;
        % pdf_output(f2, sprintf('%s%s/fwd_V', outpath, name));
        pdf_output(f2, sprintf('%s/../fwd_V', paths.fwd));
    end

    if plotFlag.true == true
        % true slip velocity
        f4 = figure;
        subplot(1,1,1);
        hold on;
        s4 = pcolor(X, T, V_true);
        s4.EdgeColor = 'none';
        % contour(X, T, V_true, [1e-5, 1], 'k', 'linewidth', 1.5);
        axis([-10, 10, 0, Tend]);
        xlabel('x along fault');
        ylabel('t');
        title("V: true");
        caxis([-2, 2]);
        elastic.cmap;
        colorbar;
        % pdf_output(f4, sprintf('%s%s/true_V', outpath, name));
        pdf_output(f4, sprintf('%s/../true_V', paths.true));

        % residual slip
        f5 = figure;
        subplot(1,1,1);
        hold on;
        s5 = pcolor(X, T, V_fwd - V_true);
        s5.EdgeColor = 'none';
        % contour(X, T, s_fwd,  [1e-5, 0.2], 'k', 'linewidth', 1.5);
        if contFlag.Trup == true
            contour(X, T, V_fwd > 0, 'k', 'linewidth', 0.1);
        elseif contFlag.Trups == true
            contour(X, T, s_fwd > 0, 'k', 'linewidth', 1.5);
        end
        if contFlag.Tdc == true
            contour(X, T, s_fwd > contFlag.Dc, 'k', 'linewidth', 1.5);
        end
        axis([-10, 10, 0, Tend]);
        xlabel('x along fault');
        ylabel('t');
        title("residual: V_{fwd} - V_{true}");
        cmax = max(abs(V_fwd(:) - V_true(:))) + 1e-10;
        caxis([-cmax, cmax]);
        elastic.cmap;
        colorbar;
        % pdf_output(f5, sprintf('%s%s/res_V', outpath, name));
        pdf_output(f5, sprintf('%s/../res_V', paths.fwd));
    end

    if plotFlag.adj == true
        % adjoint strength
        f6 = figure;
        subplot(1,1,1);
        hold on;
        s6 = pcolor(X_adj, T_adj, taustr_adj);
        % s.FaceColor = 'interp';
        s6.EdgeColor = 'none';
        % contour(X, T, s_fwd,  [1e-5, 0.2], 'k', 'linewidth', 1.5);
        if contFlag.Trup == true
            contour(X, T, V_fwd > 0, 'k', 'linewidth', 0.1);
        elseif contFlag.Trups == true
            contour(X, T, s_fwd > 0, 'k', 'linewidth', 1.5);
        end
        if contFlag.Tdc == true
            contour(X, T, s_fwd > contFlag.Dc, 'k', 'linewidth', 1.5);
        end
        xlabel('x along fault');
        ylabel('t');
        title('Adjoint strength');
        cmax = max(abs(taustr_adj(:)));
        caxis([-cmax, cmax]);
        elastic.cmap;
        colorbar;
        % pdf_output(f6, sprintf('%s%s/adj_str', outpath, name));
        pdf_output(f6, sprintf('%s/../adj_str', paths.adj));

        % adjoint slip velocity
        f7 = figure;
        subplot(1,1,1);
        hold on;
        s7 = pcolor(X_adj, T_adj, V_adj);
        % s.FaceColor = 'interp';
        s7.EdgeColor = 'none';
        % contour(X, T, s_fwd,  [1e-5, 0.2], 'k', 'linewidth', 1.5);
        if contFlag.Trup == true
            contour(X, T, V_fwd > 0, 'k', 'linewidth', 0.1);
        elseif contFlag.Trups == true
            contour(X, T, s_fwd > 0, 'k', 'linewidth', 1.5);
        end
        if contFlag.Tdc == true
            contour(X, T, s_fwd > contFlag.Dc, 'k', 'linewidth', 1.5);
        end
        xlabel('x along fault');
        ylabel('t');
        title('Adjoint slip velocity');
        if isempty(adjcmax)
            cmax = max(abs(V_adj(:))) * adjcmax + 1e-10;
        else
            cmax = adjcmax;
        end
        caxis([-cmax, cmax]);
        elastic.cmap;
        colorbar;
        % pdf_output(f7, sprintf('%s%s/adj_V', outpath, name));
        pdf_output(f7, sprintf('%s/../adj_V', paths.adj));

        % gradient
        % load(sprintf('%s%s/adj/grad.mat', outpath, name), 'x', 'grad');
        load(sprintf('%s/grad.mat', paths.adj), 'x', 'grad');
        faultpx = x;

        f8 = figure();
        hold on;
        optvars = fieldnames(grad);
        for j = 1:length(optvars)
            plot(faultpx, grad.(optvars{j}), 'linewidth', 2, 'DisplayName', optvars{j});
        end
        % plot(faultpx,grad.tau0,'linewidth',2, 'DisplayName', '$\tau_{0}$');
        % plot(faultpx,grad.taup,'linewidth',2, 'DisplayName', '$\tau_{p}$');
        % plot(faultpx,grad.taur,'linewidth',2, 'DisplayName', '$\tau_{r}$');
        % plot(faultpx,grad.D_c,'linewidth',2, 'DisplayName', '$D_c$');
        xlabel('x along fault (km)');
        % legend('show', 'interpreter', 'latex', 'Location', 'best');
        legend('show', 'Location', 'best');
        grid on;
        set(gca,'fontsize',16);
        title('Gradient');
        % pdf_output(f8, sprintf('%s%s/grad', outpath, name));
        pdf_output(f8, sprintf('%s/../grad', paths.adj));

        % tadj-Vadj plot
        plotx = linspace(-3, 3, 7);
        f9 = figure();
        hold on;
        for i = 1:length(plotx)
            xiter = find(faultx > plotx(i), 1, 'first');
            if isempty(xiter)
                continue;
            end
            plot(T_adj(xiter,:), V_adj(xiter,:), 'DisplayName', sprintf('x = %.1f', plotx(i)), 'linewidth', 2);
        %     plot(T(xiter,:), V_fwd(xiter,:), 'DisplayName', sprintf('x = %.2f', faultx(xiter)), 'linewidth', 2, 'linestyle', ':');
        %     plot(T(xiter,:), V_true(xiter,:), 'DisplayName', sprintf('x = %.2f', faultx(xiter)), 'linewidth', 2);
        end
        xlabel('$t$ (s)', 'interpreter', 'latex');
        ylabel('$V^{\dagger}$ (m/s)', 'interpreter', 'latex');
        legend('show', 'interpreter', 'latex', 'Location', 'best');
        grid on;
        set(gca,'fontsize',16);
        xlim([0, Tend]);
        % pdf_output(f9, sprintf('%s%s/t_V', outpath, name));
        pdf_output(f9, sprintf('%s/../adj_t_V', paths.adj));
    end
end
