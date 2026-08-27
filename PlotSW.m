classdef PlotSW < handle
    properties
        r
        r0
        type
        grid
        faultx
        update
        fh
    end

    methods
        function obj = PlotSW(ini_path, type)
            default_arg('type',struct);
            % default_field(type,'plot_variables','all');
            default_field(type,'axlims', []);
            default_field(type,'umax', []);
            default_field(type,'vmax', []);
            default_field(type,'sim', 'fwd');
            default_field(type,'plot_variables', 'all'); % or 'vVS'
            obj.r = [];
            load(ini_path);
            obj.r0 = saveData; 
            obj.grid = obj.r0.grid; % Assuming grid is stored in the loaded file
            obj.faultx = full(obj.r0.fault(:,1)); % Assuming faultx is the first column of the fault matrix
            obj.type = type;
            % [obj.update, obj.fh] = obj.setupPlot(type);
            obj.update = [];
            obj.fh = []; 
        end

        function plot_result(obj, path)
            load(path);
            obj.r = saveData;
            if isempty(obj.r)
                error('Results not loaded. Call load_result first.');
            end
            % 初回のみfigureを作成
            if isempty(obj.fh)
                [obj.update, obj.fh] = obj.setupPlot(obj.type);
            end
            % 既存のfigureを使用
            figure(obj.fh);
            obj.update(obj.r);
            drawnow;
        end

        function [update,figure_handle] = setupPlot(obj, type)
            f = figure();

            if type.plot_variables == "all"
                f.Position = [100, 100, 300, 1000]; % Set figure size
                tiledlayout(6,1);
                [Surs, Lines] = setup_all(obj.r0);
                update = @(r) update_all(r, Surs, Lines, type);
            elseif type.plot_variables == "vVS"
                f.Position = [100, 100, 360, 500]; % Set figure size
                tiledlayout(5,1);
                [Surs, Lines] = setup_vVS(obj.r0);
                update = @(r) update_vVS(r, Surs, Lines, type);
            end
            figure_handle = f;

            function [Surs, Lines] = setup_vVS(r)
                nexttile(1,[3,1])
                Sur1 = plot_comp(r.v);
                if type.sim == "adj"
                    title(sprintf('$du^{\\dagger} / dt^{\\dagger}, t^{\\dagger} = %.2e$', r.t), "Interpreter", "latex");
                else
                    title(sprintf('$v, t = %.2e$', r.t), "Interpreter", "latex");
                end
                pbaspect([1 1 1]);

                nexttile(5,[1,1])
                line1 = plot(obj.faultx,r.s, 'linewidth',  2);
                if type.sim == "adj"
                    ylabel('$\tau_{str}^{\dagger}$','interpreter','latex');
                else
                    ylabel('$s^*$','interpreter','latex');
                end
                % ylabel('$s^*$','interpreter','latex');
                xlabel('x along fault');
                ytickformat('%.2f');

                % subplot(4,1,4);
                nexttile(4,[1,1])
                if type.sim == "adj"
                    ylab = '$V^{\dagger}$';
                else
                    ylab = '$V$';
                end
                % ylab = '$V^*$';                
                line2 = plot(obj.faultx, r.V, 'linewidth',  2);
                xlabel('x along fault');
                ylabel(ylab,'interpreter','latex');
                ytickformat('%.2f');

                Surs = {Sur1};
                Lines = {line1, line2};
            end

            function update_vVS(r, Surs, Lines, type)
                nexttile(1,[3,1]);
                Surs{1}.ZData = r.v;
                Surs{1}.CData = r.v;
                if ~isempty(type.vmax)
                    cmax = type.vmax;
                else
                    cmax = max(abs(r.v)) + 1e-10;
                end
                caxis([-cmax, cmax]);
                if type.sim == "adj"
                    title(sprintf('$du^{\\dagger} / dt^{\\dagger}, t^{\\dagger} = %.2e$', r.t), "Interpreter", "latex");
                else
                    title(sprintf('$v, t = %.2e$', r.t), "Interpreter", "latex");
                end
                if ~isempty(type.axlims)
                    axis(type.axlims)
                end

                nexttile(5,[1,1]);
                Lines{1}.YData = r.s; 

                nexttile(4,[1,1]);
                Lines{2}.YData = r.V;
            end

            function [Surs, Lines] = setup_all(r)
                % V_plot = [V_plot, obj.slipRate(r.v)];
                % subplot(4,1,1);
                % tiledlayout(6,1)
                nexttile(1,[2,1])
                Sur1 = plot_comp(r.u);
                if type.sim == "adj"
                    title(sprintf('$u^{\\dagger}, t^{\\dagger} = %.2e$', r.t), "Interpreter", "latex");
                else
                    title(sprintf('$u, t = %.2e$', r.t), "Interpreter", "latex");
                end
                % title(sprintf('u, t = %.2e', r.t));
                pbaspect([2 2 1]);
                % expand the size of y-axis
                
                % subplot(4,1,2);
                nexttile(3,[2,1])
                Sur2 = plot_comp(r.v);
                if type.sim == "adj"
                    title('$du^{\dagger} / dt^{\dagger} $', "Interpreter", "latex");
                else
                    title('$v$', "Interpreter", "latex");
                end
                % title(sprintf('$v$', "Interpreter", "latex"));
                pbaspect([1 1 1]);

                % subplot(4,1,3);
                nexttile(5,[1,1])
                %%% modified by R.F.
                % line1 = plot(r.Psi, 'linewidth',  2);
                % ylabel('$\Psi$','interpreter','latex');
                % ↓
                line1 = plot(obj.faultx,r.s, 'linewidth',  2);
                if type.sim == "adj"
                    ylabel('$\tau_{str}^{\dagger}$','interpreter','latex');
                else
                    ylabel('$s^*$','interpreter','latex');
                end
                % ylabel('$s^*$','interpreter','latex');
                %%%
                xlabel('x along fault');
                ytickformat('%.2f');

                % subplot(4,1,4);
                nexttile(6,[1,1])
                if type.sim == "adj"
                    ylab = '$V^{\dagger}$';
                else
                    ylab = '$V$';
                end
                % ylab = '$V^*$';
                %ylab = '$\log(|V^*|)$';
                %%% modified by R.F.
                % w = [r.u; r.v; r.Psi; r.us; r.usim; r.usip];% r.Vs];
                % ↓
                % w = [r.u; r.v; r.s; r.us; r.usim; r.usip];% r.Vs];
                %%%
                
                % Vs = obj.V_star(r.t, w);
                %line2 = plot(log10(abs(Vs)), 'linewidth',  2);
                line2 = plot(obj.faultx, r.V, 'linewidth',  2);
                %line2 = plot(log10(abs(r.Vs)), 'linewidth',  2);
                xlabel('x along fault');
                ylabel(ylab,'interpreter','latex');
                ytickformat('%.2f');

                Surs = {Sur1, Sur2};
                Lines = {line1, line2};
            end

            function S = plot_comp(comp)
                S = multiblock.Surface(obj.grid, comp);
                xlabel('$x$','interpreter','latex');
                ylabel('$y$','interpreter','latex');
                elastic.cmap;
                colorbar;
                shading interp;
            end

            function update_all(r, Surs, Lines, type)
                % nexttileを使用して、setupPlotと一致させる
                nexttile(1,[2,1]);
                Surs{1}.ZData = r.u;
                Surs{1}.CData = r.u;
                % cmax = max(abs(r.u)) + 1e-10;
                if ~isempty(type.umax)
                    cmax = type.umax;
                else
                    cmax = max(abs(r.u)) + 1e-10;
                end
                caxis([-cmax, cmax]);
                if type.sim == "adj"
                    title(sprintf('$u^{\\dagger}, t^{\\dagger} = %.2e$', r.t), "Interpreter", "latex");
                else
                    title(sprintf('$u, t = %.2e$', r.t), "Interpreter", "latex");
                end
                % title(sprintf('u, t = %.2e', r.t));
                if ~isempty(type.axlims)
                    axis(type.axlims)
                end

                nexttile(3,[2,1]);
                Surs{2}.ZData = r.v;
                Surs{2}.CData = r.v;
                % cmax = max(abs(r.v)) + 1e-10;
                if ~isempty(type.vmax)
                    cmax = type.vmax;
                else
                    cmax = max(abs(r.v)) + 1e-10;
                end
                caxis([-cmax, cmax]);
                if ~isempty(type.axlims)
                    axis(type.axlims)
                end

                nexttile(5,[1,1]);
                %%% modified by R.F.
                % Lines{1}.YData = r.Psi; 
                % ↓
                Lines{1}.YData = r.s; 
                %%%           

                nexttile(6,[1,1]);
                %%% modified by R.F.
                % w = [r.u; r.v; r.Psi; r.us; r.usim; r.usip];% r.Vs];
                % ↓
                % w = [r.u; r.v; r.s; r.us; r.usim; r.usip];% r.Vs];
                %%%            
                % Vs = obj.V_star(r.t, w);
                %Lines{2}.YData = log10(abs(Vs));
                Lines{2}.YData = r.V;
                % Lines{2}.YData = Vs;
                % Lines{2}.YData = log10(abs(r.Vs));
            end



        end
    end
end