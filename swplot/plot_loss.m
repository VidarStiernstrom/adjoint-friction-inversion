function plot_loss(loadOpts, fontsize)

    load(sprintf('%s%s', loadOpts.savepath, loadOpts.filename), 'history');
    loss = history.fval(2:end-1);
    iter = 0:length(loss)-1;

    f = figure;
    plot(iter, loss, '-x', 'linewidth', 2);
    xlabel('Iteration');
    ylabel('Loss');
    set(gca, 'YScale', 'log');
    set(gca, 'fontsize', fontsize);
    grid on;
    pdf_output(f, sprintf('%s%s', loadOpts.savepath, 'loss.pdf'));

end
