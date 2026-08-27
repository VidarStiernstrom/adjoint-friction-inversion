function plot_pars(pars, parlabels, parstyles, yaxname, m_p, savepath, figsize)
    default_arg('yaxname', '$Stress$ (MPa)')
    default_arg('savepath', '');
    default_arg('figsize', [])
    x = linspace(-10, 10, m_p); % Parameter grid
    % pars = {@tau0, @taur, @taup};

    if isempty(figsize)
        f = figure();
    else
        f = figure('Position', figsize);
    end
    legend_labels = cell(length(pars), 1);
    for i = 1:length(pars)
        if isa(pars{i}, 'function_handle')
            tau_vec_i = pars{i}(x, 0); % Assuming y=0 for 2D horizontal fault
        else
            if isscalar(pars{i})
                tau_vec_i = pars{i} * ones(size(x)); % Constant value
            else
                % If pars{i} is a vector, use it directly
                % Assuming pars{i} is already in the correct format for plotting against x
                assert(length(pars{i}) == m_p, 'Parameter vector length does not match m_p');
                tau_vec_i = pars{i};
            end
        end
        % tau_vec_i = pars{i}(x, 0);
        plot(x, tau_vec_i, parstyles{i}, 'linewidth', 2);
        hold on;
        legend_labels{i} = sprintf('%s', parlabels{i});
    end
    legend(legend_labels, 'interpreter', 'latex', 'Location', 'southwest');
    xlabel('$x$ (km)', 'interpreter', 'latex');
    ylabel(yaxname, 'interpreter', 'latex');
    set(gca, 'fontsize', 16);
    grid on;
    if ~isempty(savepath)
        pdf_output(f, savepath);
    end
end
