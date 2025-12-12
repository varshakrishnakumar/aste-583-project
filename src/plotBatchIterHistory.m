function plotBatchIterHistory(iterHist, titleStr, figName)
    if nargin < 2 || isempty(titleStr)
        titleStr = 'Batch iteration history';
    end
    if nargin < 3 || isempty(figName)
        figName = 'Batch_IterHist';
    end
    it  = iterHist(:,1);
    dX  = iterHist(:,2);
    rr  = iterHist(:,4);
    rho = iterHist(:,5);
    f = figure('Name', figName);
    tlo = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
    nexttile(tlo);
    semilogy(it, dX, '-o', 'LineWidth', 1.5, 'MarkerSize', 6);
    grid on;
    xlabel('Iteration');
    ylabel('||\Delta X||');
    title(sprintf('%s: state update norm', titleStr));
    nexttile(tlo);
    yyaxis left;
    plot(it, rr, '-o', 'LineWidth', 1.5, 'MarkerSize', 6);
    ylabel('RR RMS [mm/s]');
    yyaxis right;
    plot(it, rho, '-s', 'LineWidth', 1.5, 'MarkerSize', 6);
    ylabel('Range RMS [m]');
    grid on;
    xlabel('Iteration');
    title(sprintf('%s: residual RMS', titleStr));
    legend('RR RMS','Range RMS','Location','best');
end