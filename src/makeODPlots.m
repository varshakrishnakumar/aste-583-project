function figs = makeODPlots(out, opts)
    if nargin < 2, opts = struct(); end
    opts = setDef(opts,'zoomPct',99);
    opts = setDef(opts,'save',false);
    opts = setDef(opts,'outDir',"");
    opts = setDef(opts,'prefix',"");
    opts = setDef(opts,'showHist',true);
    opts = setDef(opts,'showCorr',true);
    isEKF = isfield(out,'x_hist') && isfield(out,'P_hist') && isfield(out,'nis');
    isBatch = isfield(out,'res_prefit') && isfield(out,'iterHist');
    if ~isEKF && ~isBatch
        error('makeODPlots:UnknownOut', ...
            'Could not detect EKF or Batch output struct. Missing expected fields.');
    end
    if isEKF
        t_days   = (out.tk(:) - out.ET0) / 86400;
        stIDs    = out.station(:);
        pre_rr   = out.res_pre_rr(:);
        post_rr  = out.res_post_rr(:);
        pre_rho  = out.res_pre_rho(:);
        post_rho = out.res_post_rho(:);
        nis      = out.nis(:);
        if isfield(out,'dof')
            dof = out.dof(:);
        else
            dof = ones(size(nis));
            dof(isfinite(pre_rho)) = 2;
        end
        x_hist = out.x_hist;
        P_hist = out.P_hist;
        x_final = out.x_final(:);
        P_final = out.P_final;
        kind = "EKF";
    else
        t_days   = (out.tMeas(:) - out.ET0) / 86400;
        stIDs    = out.stationIDs(:);
        pre_rr   = out.res_prefit.rr(:);
        post_rr  = out.res_postfit.rr(:);
        pre_rho  = out.res_prefit.rho(:);
        post_rho = out.res_postfit.rho(:);
        x_final = out.X_est(:);
        P_final = out.P_est;
        kind = "BATCH";
    end
    uSt = unique(stIDs(:)).';
    C   = lines(numel(uSt));
    figs = gobjects(0);
    f = figure('Color','w','Name', kind + " RR residuals by station");
    figs(end+1) = f;
    if hasTiled()
        tlo = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
        ax1 = nexttile; ax2 = nexttile;
    else
        ax1 = subplot(2,1,1);
        ax2 = subplot(2,1,2);
    end
    plotByStation(ax1, t_days, stIDs, pre_rr*1e3, uSt, C);
    styleResidualAxes(ax1, 'Days since detection', 'Prefit range-rate resid [m/s]', opts.zoomPct);
    title(ax1, kind + " prefit range-rate residuals", 'Interpreter','none');
    plotByStation(ax2, t_days, stIDs, post_rr*1e3, uSt, C);
    styleResidualAxes(ax2, 'Days since detection', 'Postfit range-rate resid [m/s]', opts.zoomPct);
    title(ax2, kind + " postfit range-rate residuals", 'Interpreter','none');
    addLegendOutside(ax2, uSt);
    if any(isfinite(pre_rho)) || any(isfinite(post_rho))
        f = figure('Color','w','Name', kind + " Range residuals by station");
        figs(end+1) = f;
        if hasTiled()
            tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
            ax1 = nexttile; ax2 = nexttile;
        else
            ax1 = subplot(2,1,1);
            ax2 = subplot(2,1,2);
        end
        plotByStation(ax1, t_days, stIDs, pre_rho*1e3, uSt, C);
        styleResidualAxes(ax1, 'Days since detection', 'Prefit range resid [m]', opts.zoomPct);
        title(ax1, kind + " prefit range residuals", 'Interpreter','none');
        plotByStation(ax2, t_days, stIDs, post_rho*1e3, uSt, C);
        styleResidualAxes(ax2, 'Days since detection', 'Postfit range resid [m]', opts.zoomPct);
        title(ax2, kind + " postfit range residuals", 'Interpreter','none');
        addLegendOutside(ax2, uSt);
    end
    [sig_rho, sig_rr] = defaultSigmas(kind, stIDs);
    e_rr  = post_rr ./ sig_rr;
    e_rho = post_rho ./ sig_rho;
    f = figure('Color','w','Name', kind + " Standardized postfit residuals");
    figs(end+1) = f;
    if hasTiled()
        tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
        ax1 = nexttile; ax2 = nexttile; ax3 = nexttile; ax4 = nexttile;
    else
        ax1 = subplot(2,2,1); ax2 = subplot(2,2,2);
        ax3 = subplot(2,2,3); ax4 = subplot(2,2,4);
    end
    plotByStation(ax1, t_days, stIDs, e_rr, uSt, C);
    styleStdAxes(ax1, 'Days since detection', 'Std postfit range-rate resid [σ]');
    title(ax1, 'Standardized RR (postfit)', 'Interpreter','none');
    histogram(ax2, e_rr(isfinite(e_rr)), 60); grid(ax2,'on');
    xline(ax2, 0,'k-'); xline(ax2, 3,'k--'); xline(ax2, -3,'k--');
    xlabel(ax2, 'Std postfit RR resid [σ]'); ylabel(ax2,'count');
    title(ax2, 'Histogram (RR)', 'Interpreter','none');
    if any(isfinite(e_rho))
        plotByStation(ax3, t_days, stIDs, e_rho, uSt, C);
        styleStdAxes(ax3, 'Days since detection', 'Std postfit range resid [σ]');
        title(ax3, 'Standardized range (postfit)', 'Interpreter','none');
        histogram(ax4, e_rho(isfinite(e_rho)), 60); grid(ax4,'on');
        xline(ax4, 0,'k-'); xline(ax4, 3,'k--'); xline(ax4, -3,'k--');
        xlabel(ax4, 'Std postfit range resid [σ]'); ylabel(ax4,'count');
        title(ax4, 'Histogram (range)', 'Interpreter','none');
    else
        axis(ax3,'off'); axis(ax4,'off');
        text(ax3, 0.1, 0.5, 'No range measurements in this dataset.', 'FontSize', 12);
    end
    f = figure('Color','w','Name', kind + " RMS by station");
    figs(end+1) = f;
    rr_pre_rms   = nan(numel(uSt),1);
    rr_post_rms  = nan(numel(uSt),1);
    rho_pre_rms  = nan(numel(uSt),1);
    rho_post_rms = nan(numel(uSt),1);
    for i = 1:numel(uSt)
        ii = (stIDs == uSt(i));
        rr_pre_rms(i)   = rms_omitnan(pre_rr(ii)*1e3);
        rr_post_rms(i)  = rms_omitnan(post_rr(ii)*1e3);
        rho_pre_rms(i)  = rms_omitnan(pre_rho(ii)*1e3);
        rho_post_rms(i) = rms_omitnan(post_rho(ii)*1e3);
    end
    if hasTiled()
        tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
        ax1 = nexttile; ax2 = nexttile;
    else
        ax1 = subplot(1,2,1); ax2 = subplot(1,2,2);
    end
    bar(ax1, [rr_pre_rms, rr_post_rms]); grid(ax1,'on');
    set(ax1,'XTick',1:numel(uSt),'XTickLabel',stationLegend(uSt), ...
        'XTickLabelRotation',25);
    ylabel(ax1,'RMS RR [m/s]');
    legend(ax1, {'prefit','postfit'}, 'Location','best');
    title(ax1,'Range-rate RMS by station','Interpreter','none');
    set(ax1,'YScale','log');
    yl_rr = [max(min(rr_post_rms)/2, 1e-6), max(rr_pre_rms)*2];
    ylim(ax1, yl_rr);
    if any(isfinite(rho_pre_rms)) || any(isfinite(rho_post_rms))
        bar(ax2, [rho_pre_rms, rho_post_rms]); grid(ax2,'on');
        set(ax2,'XTick',1:numel(uSt),'XTickLabel',stationLegend(uSt), ...
            'XTickLabelRotation',25);
        ylabel(ax2,'RMS range [m]');
        legend(ax2, {'prefit','postfit'}, 'Location','best');
        title(ax2,'Range RMS by station','Interpreter','none');
        set(ax2,'YScale','log');
        yl_rho = [max(min(rho_post_rms)/2, 1e-2), max(rho_pre_rms)*2];
        ylim(ax2, yl_rho);
    else
        axis(ax2,'off');
        text(ax2, 0.1, 0.5, 'No range measurements in this dataset.', 'FontSize', 12);
    end
    if isEKF
        f = figure('Color','w','Name','NIS');
        figs(end+1) = f;
        semilogy(t_days, nis, '.'); grid on; hold on;
        alpha = 0.05;
        lo = arrayfun(@(k) chi2inv_local(alpha/2, max(1,round(dof(k)))), 1:numel(dof)).';
        hi = arrayfun(@(k) chi2inv_local(1-alpha/2, max(1,round(dof(k)))), 1:numel(dof)).';
        semilogy(t_days, lo, 'k--');
        semilogy(t_days, hi, 'k--');
        xlabel('Days since detection');
        ylabel('NIS');
        title('NIS with 95% chi-square bounds','Interpreter','none');
        legend({'NIS','lower bound','upper bound'},'Location','best');
    end
    if isBatch
        IH = out.iterHist;
        it = IH(:,1);
        dx = IH(:,2);
        cost = IH(:,3);
        rrmm = IH(:,4);
        rhom = IH(:,5);
        f = figure('Color','w','Name','Batch convergence');
        figs(end+1) = f;
        if hasTiled()
            tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
            ax1 = nexttile; ax2 = nexttile; ax3 = nexttile; ax4 = nexttile;
        else
            ax1 = subplot(2,2,1); ax2 = subplot(2,2,2);
            ax3 = subplot(2,2,3); ax4 = subplot(2,2,4);
        end
        semilogy(ax1, it, cost, '-o'); grid(ax1,'on'); xlabel(ax1,'iter'); ylabel(ax1,'cost');
        title(ax1,'Cost','Interpreter','none');
        semilogy(ax2, it, dx, '-o'); grid(ax2,'on'); xlabel(ax2,'iter'); ylabel(ax2,'||dX||');
        title(ax2,'Step norm','Interpreter','none');
        plot(ax3, it, rrmm, '-o'); grid(ax3,'on'); xlabel(ax3,'iter'); ylabel(ax3,'RR RMS [mm/s]');
        title(ax3,'RR RMS','Interpreter','none');
        plot(ax4, it, rhom, '-o'); grid(ax4,'on'); xlabel(ax4,'iter'); ylabel(ax4,'Range RMS [m]');
        title(ax4,'Range RMS','Interpreter','none');
    end
    if opts.showCorr && ~isempty(P_final)
        f = figure('Color','w','Name', kind + " Covariance correlation");
        figs(end+1) = f;
        corr = cov2corr(P_final);
        imagesc(corr); axis equal tight; colorbar; grid on;
        title('Posterior correlation matrix','Interpreter','none');
        labels = {'x','y','z','vx','vy','vz','k','lat4','lon4','bias'};
        set(gca,'XTick',1:10,'YTick',1:10,'XTickLabel',labels,'YTickLabel',labels);
        xtickangle(35);
    end
    if isfield(out,'X0') && ~isempty(out.X0)
        X0 = out.X0(:);
        dr = x_final(1:3) - X0(1:3);
        f = figure('Color','w','Name', kind + " Δr_xy + ellipse");
        figs(end+1) = f;
        hold on; grid on;
        plot(0,0,'ko','MarkerFaceColor','k','DisplayName','prior origin');
        scatter(dr(1), dr(2), 40, 'filled', 'DisplayName','\Delta r at ET0');
        xlabel('\Delta r_x [km]'); ylabel('\Delta r_y [km]');
        title('\Delta r_{XY} at ET0 with 3\sigma ellipse','Interpreter','tex');
        Pxy = P_final(1:2,1:2);
        plotCovEllipse(dr(1:2), Pxy, 3);
        sig_x = sqrt(max(P_final(1,1),0));
        sig_y = sqrt(max(P_final(2,2),0));
        xlim([dr(1)-4*sig_x, dr(1)+4*sig_x]);
        ylim([dr(2)-4*sig_y, dr(2)+4*sig_y]);
        axis equal;
        legend('Location','best');
    end
    if opts.save
        if ~isfield(opts,'outDir') || strlength(string(opts.outDir))==0
            opts.outDir = fullfile(pwd, "figs_" + kind);
        end
        if ~exist(opts.outDir,'dir'), mkdir(opts.outDir); end
        if ~isfield(opts,'prefix') || strlength(string(opts.prefix))==0
            opts.prefix = kind;
        end
        saveFigureSet(figs, opts.outDir, char(opts.prefix));
        fprintf('Saved %d figures to: %s\n', numel(figs), opts.outDir);
    end
end
function tf = hasTiled()
    tf = exist('tiledlayout','file') == 2;
end
function plotByStation(ax, t_days, stIDs, y, uSt, C)
    axes(ax);
    hold(ax,'on'); grid(ax,'on');
    for i = 1:numel(uSt)
        ii = (stIDs == uSt(i)) & isfinite(y);
        h = scatter(ax, t_days(ii), y(ii), 10, C(i,:), 'filled');
        try
            set(h,'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.6);
        catch
        end
    end
    yline(ax,0,'k-');
end
function styleResidualAxes(ax, xlab, ylab, zoomPct)
    xlabel(ax, xlab, 'Interpreter','none');
    ylabel(ax, ylab, 'Interpreter','none');
    kids = ax.Children;
    yy = [];
    for k = 1:numel(kids)
        if isprop(kids(k),'YData')
            yy = [yy; kids(k).YData(:)];
        end
    end
    yl = robustSymLim(yy, zoomPct);
    ylim(ax, yl);
end
function styleStdAxes(ax, xlab, ylab)
    xlabel(ax, xlab, 'Interpreter','none');
    ylabel(ax, ylab, 'Interpreter','none');
    yline(ax, 3,'k--'); yline(ax, -3,'k--');
    ylim(ax, [-6 6]);
end
function yl = robustSymLim(y, pct)
    y = y(isfinite(y));
    if isempty(y)
        yl = [-1 1];
        return;
    end
    a = sort(abs(y(:)));
    n = numel(a);
    idx = max(1, min(n, ceil((pct/100)*n)));
    m = a(idx);
    if m == 0
        m = max(a);
        if m == 0, m = 1; end
    end
    yl = [-m m];
end
function addLegendOutside(ax, uSt)
    legend(ax, stationLegend(uSt), 'Location','bestoutside');
end
function s = stationLegend(uSt)
    s = cell(1,numel(uSt));
    for i=1:numel(uSt)
        s{i} = sprintf('#%d %s', uSt(i), stationName(uSt(i)));
    end
end
function name = stationName(id)
    switch id
        case 1, name = 'Goldstone';
        case 2, name = 'Canberra';
        case 3, name = 'Madrid';
        case 4, name = 'Antarctica';
        otherwise, name = 'Unknown';
    end
end
function [sig_rho, sig_rr] = defaultSigmas(kind, stIDs)
    sigma_rho_dsn = 1e-3;
    sigma_rho_ant = 1e-2;
    sigma_rr_dsn  = 1e-7;
    sigma_rr_ant  = 1e-6;
    if kind == "EKF"
        Rscale_rr_dsn  = 5;
        Rscale_rr_ant  = 10;
        Rscale_rho_dsn = 5;
        Rscale_rho_ant = 5;
    else
        Rscale_rr_dsn  = 5;
        Rscale_rr_ant  = 5;
        Rscale_rho_dsn = 5;
        Rscale_rho_ant = 5;
    end
    sig_rho = nan(size(stIDs));
    sig_rr  = nan(size(stIDs));
    for k = 1:numel(stIDs)
        if stIDs(k) == 4
            sig_rho(k) = Rscale_rho_ant * sigma_rho_ant;
            sig_rr(k)  = Rscale_rr_ant  * sigma_rr_ant;
        else
            sig_rho(k) = Rscale_rho_dsn * sigma_rho_dsn;
            sig_rr(k)  = Rscale_rr_dsn  * sigma_rr_dsn;
        end
    end
end
function x = chi2inv_local(p, k)
    x = 2 * gammaincinv(p, k/2);
end
function v = rms_omitnan(x)
    x = x(isfinite(x));
    if isempty(x), v = NaN; else, v = sqrt(mean(x.^2)); end
end
function C = cov2corr(P)
    P = 0.5*(P+P.');
    d = sqrt(max(diag(P),0));
    D = d*d.';
    C = P ./ D;
    C(~isfinite(C)) = 0;
    C = max(min(C,1),-1);
end
function plotCovEllipse(mu, P2, nsig)
    P2 = 0.5*(P2+P2.');
    [V,D] = eig(P2);
    D = max(D, 0);
    A = V*sqrt(D);
    th   = linspace(0, 2*pi, 200);
    circ = [cos(th); sin(th)];
    ell  = mu(:) + nsig * (A * circ);
    plot(ell(1,:), ell(2,:), 'k--', 'LineWidth',1.2, ...
        'DisplayName', sprintf('%g-sigma ellipse', nsig));
end
function saveFigureSet(figs, outDir, prefix)
    for i = 1:numel(figs)
        f = figs(i);
        if ~isvalid(f), continue; end
        nm = f.Name;
        if isempty(nm), nm = sprintf('Figure%d', f.Number); end
        safe = regexprep(nm, '[^a-zA-Z0-9_-]', '_');
        base = sprintf('%s_%02d_%s', prefix, i, safe);
        pngPath = fullfile(outDir, [base '.png']);
        pdfPath = fullfile(outDir, [base '.pdf']);
        set(f,'Color','w');
        if exist('exportgraphics','file') == 2
            exportgraphics(f, pngPath, 'Resolution', 300);
            exportgraphics(f, pdfPath, 'ContentType', 'vector');
        else
            print(f, fullfile(outDir, base), '-dpng', '-r300');
            print(f, fullfile(outDir, base), '-dpdf', '-r300');
        end
    end
end
function opts = setDef(opts, field, val)
    if ~isfield(opts, field) || isempty(opts.(field))
        opts.(field) = val;
    end
end