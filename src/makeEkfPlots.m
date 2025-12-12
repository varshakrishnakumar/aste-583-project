function figs = makeEkfPlots(span, t_days, stationID, ...
                      res_pre_rr, res_post_rr, ...
                      res_pre_rho, res_post_rho, ...
                      x_hist, P_hist, varargin)
    spanStr   = char(string(span));
    t_days    = t_days(:);
    stationID = stationID(:);
    nSteps = size(x_hist,2);
    if nSteps == numel(t_days)
        t_state = t_days;
    else
        t_state = linspace(min(t_days), max(t_days), nSteps).';
        warning('makeEkfPlots:TimeMismatch', ...
            'size(x_hist,2) ~= length(t_days). Using linspace timebase for state/cov plots.');
    end
    nis = [];
    dof = [];
    v = varargin;
    if ~isempty(v) && isnumeric(v{1})
        nis = v{1}; v(1) = [];
        if ~isempty(v) && isnumeric(v{1})
            dof = v{1}; v(1) = [];
        end
    end
    optsStruct = struct();
    if numel(v) == 1 && isstruct(v{1})
        optsStruct = v{1};
        v = {};
    end
    p = inputParser;
    p.CaseSensitive = false;
    p.addParameter('fitMask', []);
    p.addParameter('dcoDay', []);
    p.addParameter('ltmDay', 15);
    p.addParameter('covPlane', 'RT');
    p.addParameter('stateLabels', {});
    p.addParameter('saveDir', '');
    p.addParameter('prefix', '');
    p.addParameter('maxEllipses', 6);
    p.parse(v{:});
    opt = p.Results;
    opt = mergeStruct(opt, optsStruct);
    if isempty(opt.dcoDay)
        s = lower(spanStr);
        if contains(s,'0to6') || contains(s,'0-6') || contains(s,'0_6')
            opt.dcoDay = 6;
        else
            opt.dcoDay = 14;
        end
    end
    if isempty(opt.fitMask)
        opt.fitMask = true(size(t_days));
    else
        opt.fitMask = logical(opt.fitMask(:));
        if numel(opt.fitMask) ~= numel(t_days)
            error('fitMask must be the same length as t_days.');
        end
    end
    uSt = unique(stationID(isfinite(stationID)));
    [stNames, stColors, stMarkers, stIsAnt] = stationStyles(uSt);
    sigmaRR_DSN_mmps = 0.1;
    sigmaRR_ANT_mmps = 1.0;
    sigmaR_DSN_m     = 1.0;
    sigmaR_ANT_m     = 10.0;
    rr_pre_mmps  = res_pre_rr(:)  * 1e6;
    rr_post_mmps = res_post_rr(:) * 1e6;
    figs = struct();
    figs.rr = plotResidualStack( ...
        sprintf('EKF_%s_RR', spanStr), spanStr, ...
        t_days, stationID, opt.fitMask, opt.dcoDay, ...
        rr_pre_mmps, rr_post_mmps, ...
        'Range-rate residual [mm/s]', ...
        {'Prefit RR residuals', 'Postfit RR residuals'}, ...
        stNames, stColors, stMarkers, ...
        sigmaRR_DSN_mmps, sigmaRR_ANT_mmps);
    hasRange = any(isfinite(res_pre_rho(:))) || any(isfinite(res_post_rho(:)));
    if hasRange
        rho_pre_m  = res_pre_rho(:)  * 1e3;
        rho_post_m = res_post_rho(:) * 1e3;
        figs.rho = plotResidualStack( ...
            sprintf('EKF_%s_Range', spanStr), spanStr, ...
            t_days, stationID, opt.fitMask, opt.dcoDay, ...
            rho_pre_m, rho_post_m, ...
            'Range residual [m]', ...
            {'Prefit range residuals', 'Postfit range residuals'}, ...
            stNames, stColors, stMarkers, ...
            sigmaR_DSN_m, sigmaR_ANT_m);
    else
        figs.rho = figure('Name', sprintf('EKF_%s_Range', spanStr), 'Color','w');
        ax = axes(figs.rho); axis(ax,'off');
        text(ax, 0.5, 0.5, 'No range measurements in this dataset.', ...
            'HorizontalAlignment','center', 'FontSize',12);
    end
    figs.uncertainty = plotUncertainty( ...
        sprintf('EKF_%s_Uncertainty', spanStr), spanStr, ...
        t_state, x_hist, P_hist, opt.dcoDay, opt.ltmDay, opt.covPlane, opt.maxEllipses);
    figs.params = plotParamsHistory( ...
        sprintf('EKF_%s_Params', spanStr), spanStr, ...
        t_state, x_hist, P_hist, opt.stateLabels, opt.dcoDay);
    figs.corr = plotCorrelationMatrix( ...
        sprintf('EKF_%s_Corr', spanStr), spanStr, P_hist(:,:,end), opt.stateLabels);
    if ~isempty(nis) && ~isempty(dof)
        figs.nis = plotNIS( ...
            sprintf('EKF_%s_NIS', spanStr), spanStr, ...
            t_days, nis, dof, opt.dcoDay);
    end
    if ~isempty(opt.saveDir)
        saveDir = char(string(opt.saveDir));
        if ~exist(saveDir, 'dir'); mkdir(saveDir); end
        prefix = char(string(opt.prefix));
        saveAllFigures(figs, saveDir, prefix, spanStr);
    end
end
function out = mergeStruct(out, in)
    if isempty(in); return; end
    fn = fieldnames(in);
    for k = 1:numel(fn)
        out.(fn{k}) = in.(fn{k});
    end
end
function [names, colors, markers, isAnt] = stationStyles(uSt)
    baseNames = containers.Map('KeyType','double','ValueType','char');
    baseNames(1) = 'Goldstone';
    baseNames(2) = 'Canberra';
    baseNames(3) = 'Madrid';
    baseNames(4) = 'Antarctica';
    baseColors = containers.Map('KeyType','double','ValueType','any');
    baseColors(1) = [0.0000 0.4470 0.7410];
    baseColors(2) = [0.8500 0.3250 0.0980];
    baseColors(3) = [0.9290 0.6940 0.1250];
    baseColors(4) = [0.4940 0.1840 0.5560];
    baseMarkers = containers.Map('KeyType','double','ValueType','char');
    baseMarkers(1) = 'o';
    baseMarkers(2) = 'o';
    baseMarkers(3) = 'o';
    baseMarkers(4) = 's';
    names   = cell(size(uSt));
    colors  = zeros(numel(uSt),3);
    markers = cell(size(uSt));
    isAnt   = false(size(uSt));
    extra = lines(max(1,numel(uSt)));
    for i = 1:numel(uSt)
        id = uSt(i);
        if isKey(baseNames,id)
            names{i} = sprintf('#%d %s', id, baseNames(id));
            colors(i,:) = baseColors(id);
            markers{i} = baseMarkers(id);
        else
            names{i} = sprintf('#%d', id);
            colors(i,:) = extra(i,:);
            markers{i} = 'o';
        end
        isAnt(i) = (id == 4);
    end
end
function f = plotResidualStack(figName, spanStr, t, stID, fitMask, dcoDay, yPre, yPost, yLabel, tileTitles, stNames, stColors, stMarkers, sigmaDSN, sigmaANT)
    f = figure('Name', figName, 'Color','w');
    tl = tiledlayout(f, 2, 1, 'TileSpacing','compact', 'Padding','compact');
    ax1 = nexttile(tl, 1);
    plotResidualTile(ax1, t, stID, fitMask, dcoDay, yPre, yLabel, tileTitles{1}, stNames, stColors, stMarkers, sigmaDSN, sigmaANT, true);
    ax2 = nexttile(tl, 2);
    plotResidualTile(ax2, t, stID, fitMask, dcoDay, yPost, yLabel, tileTitles{2}, stNames, stColors, stMarkers, sigmaDSN, sigmaANT, false);
    linkaxes([ax1 ax2], 'x');
    sgtitle(tl, sprintf('EKF %s: %s', spanStr, erase(figName, "EKF_"+spanStr+"_")), 'FontWeight','bold');
end
function plotResidualTile(ax, t, stID, fitMask, dcoDay, y, yLabel, ttl, stNames, stColors, stMarkers, sigmaDSN, sigmaANT, showLegend)
    axes(ax);
    hold(ax,'on');
    grid(ax,'on'); grid(ax,'minor');
    ax.Layer = 'top';
    ax.Box = 'on';
    hasPass = any(~fitMask & isfinite(y) & isfinite(t));
    if hasPass
        x0 = dcoDay;
        x1 = max(t(isfinite(t)));
        if x1 > x0
        end
    end
    yline(ax, 0, '-', 'Color',[0 0 0], 'LineWidth',0.8, 'HandleVisibility','off');
    yline(ax, +3*sigmaDSN, '--', 'Color',[0.2 0.2 0.2], 'LineWidth',0.8, 'HandleVisibility','off');
    yline(ax, -3*sigmaDSN, '--', 'Color',[0.2 0.2 0.2], 'LineWidth',0.8, 'HandleVisibility','off');
    yline(ax, +3*sigmaANT, ':',  'Color',[0.2 0.2 0.2], 'LineWidth',0.9, 'HandleVisibility','off');
    yline(ax, -3*sigmaANT, ':',  'Color',[0.2 0.2 0.2], 'LineWidth',0.9, 'HandleVisibility','off');
    xline(ax, dcoDay, ':', 'DCO', 'Color',[0 0 0], 'LineWidth',1.0, ...
        'LabelVerticalAlignment','middle', 'LabelHorizontalAlignment','left');
    uSt = unique(stID(isfinite(stID)));
    for i = 1:numel(uSt)
        id = uSt(i);
        ii = findStationIndex(id, uSt, stNames);
        idxStation = (stID == id) & isfinite(y) & isfinite(t);
        if ~any(idxStation); continue; end
        k = find(uSt == id, 1, 'first');
        c = stColors(k,:);
        mk = stMarkers{k};
        idxFit  = idxStation & fitMask;
        idxPass = idxStation & ~fitMask;
        if any(idxFit)
            scatter(ax, t(idxFit), y(idxFit), 18, ...
                'Marker', mk, ...
                'MarkerFaceColor', c, ...
                'MarkerEdgeColor', c, ...
                'MarkerFaceAlpha', 0.70, ...
                'MarkerEdgeAlpha', 0.70, ...
                'DisplayName', stNames{k});
        end
        if any(idxPass)
            scatter(ax, t(idxPass), y(idxPass), 18, ...
                'Marker', mk, ...
                'MarkerFaceColor', 'none', ...
                'MarkerEdgeColor', c, ...
                'LineWidth', 0.9, ...
                'MarkerEdgeAlpha', 0.90, ...
                'HandleVisibility', 'off');
        end
    end
    yl = robustSymLim(y, 1.15, [3*sigmaDSN, 3*sigmaANT]);
    ylim(ax, yl);
    if hasPass
        x0 = dcoDay;
        x1 = max(t(isfinite(t)));
        if x1 > x0
            y0 = yl(1); y1 = yl(2);
            patch(ax, [x0 x1 x1 x0], [y0 y0 y1 y1], [0 0 0], ...
                'FaceAlpha', 0.035, 'EdgeColor','none', 'HandleVisibility','off');
            ch = ax.Children;
            ax.Children = [ch(2:end); ch(1)];
        end
    end
    xlabel(ax, 'Days since detection');
    ylabel(ax, yLabel);
    title(ax, ttl, 'FontWeight','bold');
    txt = sprintf('\\pm3\\sigma guides: DSN %.3g, Ant %.3g', 3*sigmaDSN, 3*sigmaANT);
    text(ax, 0.01, 0.97, txt, 'Units','normalized', 'VerticalAlignment','top', ...
        'FontSize', 10, 'Color',[0.15 0.15 0.15]);
    if showLegend
        lg = legend(ax, 'Location','eastoutside');
        lg.Box = 'off';
        lg.Interpreter = 'tex';
    end
end
function yl = robustSymLim(y, pad, mustInclude)
    y = y(isfinite(y));
    if isempty(y)
        m = 1;
    else
        ay = sort(abs(y));
        k  = max(1, round(0.99*numel(ay)));
        m  = ay(k);
        if m == 0
            m = max(ay(end), 1);
        end
    end
    if nargin >= 3 && ~isempty(mustInclude)
        m = max(m, max(abs(mustInclude)));
    end
    yl = pad * [-m m];
end
function k = findStationIndex(id, uSt, stNames)
    k = find(uSt == id, 1, 'first');
end
function f = plotUncertainty(figName, spanStr, t, x_hist, P_hist, dcoDay, ltmDay, covPlane, maxEllipses)
    f = figure('Name', figName, 'Color','w');
    tl = tiledlayout(f, 1, 2, 'TileSpacing','compact', 'Padding','compact');
    t   = t(:);
    n   = size(x_hist, 2);
    if numel(t) ~= n
        t = linspace(min(t), max(t), n).';
    end
    [~, idxDCO] = min(abs(t - dcoDay));
    rD = x_hist(1:3, idxDCO);
    vD = x_hist(4:6, idxDCO);
    Rhat = rD / max(norm(rD), eps);
    h    = cross(rD, vD);
    if norm(h) < 1e-12
        tmp = [1;0;0];
        if abs(dot(tmp,Rhat)) > 0.9; tmp = [0;1;0]; end
        Nhat = cross(Rhat, tmp); Nhat = Nhat / max(norm(Nhat), eps);
    else
        Nhat = h / norm(h);
    end
    That = cross(Nhat, Rhat); That = That / max(norm(That), eps);
    A_rtn = [Rhat.'; That.'; Nhat.']; % inertial -> RTN at DCO
    ax1 = nexttile(tl, 1);
    hold(ax1,'on'); grid(ax1,'on'); grid(ax1,'minor');
    ax1.Layer = 'top'; ax1.Box = 'on';
    r  = x_hist(1:3, :);
    dr = r - rD;
    switch upper(covPlane)
        case 'XY'
            p1 = dr(1, :);
            p2 = dr(2, :);
            xlabel(ax1,'\Deltax [km]');
            ylabel(ax1,'\Deltay [km]');
            planeLabel = 'XY (relative to DCO)';
            isRT = false;
        otherwise
            prtn = A_rtn * dr;
            p1   = prtn(1, :);
            p2   = prtn(2, :);
            xlabel(ax1,'\DeltaR [km]');
            ylabel(ax1,'\DeltaT (along-track) [km]');
            planeLabel = 'RT (relative to DCO)';
            isRT = true;
    end
    plot(ax1, p1, p2, '-', 'LineWidth', 1.2, 'Color',[0 0 0], ...
        'DisplayName','Mean trajectory (relative)');
    scatter(ax1, 0, 0, 35, 'filled', 'MarkerFaceColor',[0 0 0], ...
        'DisplayName','DCO');
    nE   = max(1, min(maxEllipses, n));
    idxE = unique(round(linspace(1, n, nE)));
    for ii = 1:numel(idxE)
        kk   = idxE(ii);
        mu   = [p1(kk); p2(kk)];
        Ppos = P_hist(1:3,1:3,kk);
        if isRT
            P_rtn = A_rtn * Ppos * A_rtn.';   % 3x3
            P2    = P_rtn(1:2,1:2);
        else
            P2    = Ppos(1:2,1:2);
        end
        alpha = 0.10 + 0.22*(ii/numel(idxE));
        plotCovEllipsePatch(ax1, mu, P2, 3, [0.1 0.4 0.8], alpha);
    end
    axis(ax1,'equal');
    title(ax1, sprintf('3\\sigma position covariance in %s', planeLabel), ...
        'FontWeight','bold');
    ax2 = nexttile(tl, 2);
    hold(ax2,'on'); grid(ax2,'on'); grid(ax2,'minor');
    ax2.Layer = 'top'; ax2.Box = 'on';
    sigmaR = nan(n,1);
    sigmaT = nan(n,1);
    sigmaN = nan(n,1);
    for kk = 1:n
        Ppos  = P_hist(1:3,1:3,kk);
        P_rtn = A_rtn * Ppos * A_rtn.';
        s     = sqrt(max(diag(P_rtn), 0));
        sigmaR(kk) = s(1);
        sigmaT(kk) = s(2);
        sigmaN(kk) = s(3);
    end
    plot(ax2, t, sigmaR, '-', 'LineWidth',1.3, 'DisplayName','\sigma_R');
    plot(ax2, t, sigmaT, '-', 'LineWidth',1.3, 'DisplayName','\sigma_T (along-track)');
    plot(ax2, t, sigmaN, '-', 'LineWidth',1.3, 'DisplayName','\sigma_N');
    xline(ax2, dcoDay, ':', 'DCO', 'Color',[0 0 0], 'LineWidth',1.0);
    xline(ax2, ltmDay, '--', 'LTM', 'Color',[0 0 0], 'LineWidth',1.0);
    yline(ax2, 100, '--', '100 km along-track req', ...
        'Color',[0.2 0.2 0.2], 'LineWidth',1.0);
    [~, idxL] = min(abs(t - ltmDay));
    sigT_L = sigmaT(idxL);
    text(ax2, 0.02, 0.98, ...
        sprintf('\\sigma_T at t=%.2f d: %.2f km', t(idxL), sigT_L), ...
        'Units','normalized', 'VerticalAlignment','top', 'FontSize',10);
    xlabel(ax2,'Days since detection');
    ylabel(ax2,'1\sigma position uncertainty [km]');
    title(ax2,'Position uncertainty in RTN (DCO-based RTN frame)', ...
        'FontWeight','bold');
    lg = legend(ax2,'Location','best'); lg.Box = 'off';
    sgtitle(tl, sprintf('EKF %s uncertainty diagnostics', spanStr), ...
        'FontWeight','bold');
end
function plotCovEllipsePatch(ax, mu, P2, nsig, rgb, faceAlpha)
    P2 = 0.5*(P2 + P2.');
    [V,D] = eig(P2);
    D = max(D, 0);
    A = V * sqrt(D);
    th = linspace(0, 2*pi, 240);
    circ = [cos(th); sin(th)];
    ell = mu(:) + nsig * (A * circ);
    patch(ax, ell(1,:), ell(2,:), rgb, ...
        'FaceAlpha', faceAlpha, 'EdgeColor', rgb, 'LineWidth', 1.0, ...
        'HandleVisibility','off');
end
function f = plotParamsHistory(figName, spanStr, t, x_hist, P_hist, stateLabels, dcoDay)
    f = figure('Name', figName, 'Color','w');
    n = size(x_hist,1);
    if isempty(stateLabels) || numel(stateLabels) ~= n
        stateLabels = defaultStateLabels(n);
    end
    idxParams = 7:n;
    if isempty(idxParams)
        ax = axes(f); axis(ax,'off');
        text(ax, 0.5, 0.5, 'No augmented parameters beyond [r;v] to plot.', ...
            'HorizontalAlignment','center', 'FontSize',12);
        return;
    end
    nP = numel(idxParams);
    nCols = 2;
    nRows = ceil(nP / nCols);
    tl = tiledlayout(f, nRows, nCols, 'TileSpacing','compact', 'Padding','compact');
    for k = 1:nP
        i = idxParams(k);
        ax = nexttile(tl, k);
        hold(ax,'on'); grid(ax,'on'); grid(ax,'minor'); ax.Layer='top'; ax.Box='on';
        x = x_hist(i,:).';
        s = sqrt(max(squeeze(P_hist(i,i,:)), 0));
        s3 = 3*s;
        [xPlot, yLab] = convertParamUnits(x, stateLabels{i});
        fill(ax, [t; flipud(t)], [xPlot + 3*convertScale(stateLabels{i})*s; flipud(xPlot - 3*convertScale(stateLabels{i})*s)], ...
            [0.2 0.6 0.2], 'FaceAlpha',0.10, 'EdgeColor','none', 'HandleVisibility','off');
        plot(ax, t, xPlot, '-', 'LineWidth',1.3, 'Color',[0.1 0.4 0.1], 'DisplayName','Estimate');
        xline(ax, dcoDay, ':', 'DCO', 'Color',[0 0 0], 'LineWidth',1.0, 'HandleVisibility','off');
        title(ax, stateLabels{i}, 'FontWeight','bold');
        xlabel(ax, 'Days since detection');
        ylabel(ax, yLab);
        txt = sprintf('final: %.4g \\pm %.3g (3\\sigma)', xPlot(end), 3*convertScale(stateLabels{i})*s(end));
        text(ax, 0.02, 0.98, txt, 'Units','normalized', 'VerticalAlignment','top', 'FontSize',9);
    end
    sgtitle(tl, sprintf('EKF %s augmented-parameter histories (with 3\\sigma)', spanStr), 'FontWeight','bold');
end
function labels = defaultStateLabels(n)
    if n == 10
        labels = {'x','y','z','v_x','v_y','v_z', ...
                  'k_{SRP}','lat_4','lon_4','rho_bias_dot'};
    else
        labels = arrayfun(@(i) sprintf('x_%d', i), 1:n, 'UniformOutput', false);
    end
end
function [xPlot, yLab] = convertParamUnits(x, label)
    labelLower = lower(label);
    xPlot = x;
    yLab = 'value';
    if contains(labelLower,'lat') || contains(labelLower,'lon')
        if max(abs(x)) <= 2*pi + 0.25
            xPlot = rad2deg(x);
            yLab = '[deg]';
        else
            yLab = '[deg]';
        end
        return;
    end
    if contains(labelLower,'bias') || contains(labelLower,'rho') || contains(labelLower,'rr')
        xPlot = x * 1e6;
        yLab = '[mm/s]';
        return;
    end
    if contains(labelLower,'srp') || contains(labelLower,'k')
        yLab = '[–]';
        return;
    end
end
function sc = convertScale(label)
    labelLower = lower(label);
    if contains(labelLower,'bias') || contains(labelLower,'rho') || contains(labelLower,'rr')
        sc = 1e6;
    else
        sc = 1;
    end
end
function f = plotCorrelationMatrix(figName, spanStr, P, stateLabels)
    f = figure('Name', figName, 'Color','w');
    ax = axes(f);
    n = size(P,1);
    if isempty(stateLabels) || numel(stateLabels) ~= n
        stateLabels = defaultStateLabels(n);
    end
    sig = sqrt(max(diag(P), 0));
    denom = sig * sig.';
    C = P ./ denom;
    C(~isfinite(C)) = 0;
    imagesc(ax, C, [-1 1]);
    axis(ax,'equal'); axis(ax,'tight');
    colormap(ax, redBlueMap(256));
    cb = colorbar(ax); cb.Label.String = 'corr';
    title(ax, sprintf('Posterior correlation matrix (%s)', spanStr), 'FontWeight','bold');
    set(ax, 'XTick',1:n, 'XTickLabel',stateLabels, ...
            'YTick',1:n, 'YTickLabel',stateLabels, ...
            'XTickLabelRotation',45);
    hold(ax,'on');
    for k = 0.5:1:(n+0.5)
        plot(ax, [0.5 n+0.5], [k k], '-', 'Color',[0 0 0 0.08], 'LineWidth',0.7, 'HandleVisibility','off');
        plot(ax, [k k], [0.5 n+0.5], '-', 'Color',[0 0 0 0.08], 'LineWidth',0.7, 'HandleVisibility','off');
    end
    for i = 1:n
        for j = 1:n
            v = C(i,j);
            if abs(v) > 0.55
                txtColor = [1 1 1];
            else
                txtColor = [0 0 0];
            end
            text(ax, j, i, sprintf('%.2f', v), ...
                'HorizontalAlignment','center', 'FontSize',8, 'Color',txtColor);
        end
    end
end
function cmap = redBlueMap(m)
    if nargin < 1; m = 256; end
    m2 = floor(m/2);
    blueToWhite = [linspace(0,1,m2)', linspace(0,1,m2)', ones(m2,1)];
    whiteToRed  = [ones(m-m2,1), linspace(1,0,m-m2)', linspace(1,0,m-m2)'];
    cmap = [blueToWhite; whiteToRed];
end
function f = plotNIS(figName, spanStr, t_days, nis, dof, dcoDay)
    f = figure('Name', figName, 'Color','w');
    ax = axes(f);
    hold(ax,'on'); grid(ax,'on'); grid(ax,'minor'); ax.Layer='top'; ax.Box='on';
    nis = nis(:);
    dof = dof(:);
    ok = isfinite(nis) & isfinite(dof) & (dof > 0) & isfinite(t_days);
    if ~any(ok)
        title(ax, sprintf('NIS (%s): no valid data', spanStr), 'FontWeight','bold');
        return;
    end
    t = t_days(ok);
    n = nis(ok);
    k = dof(ok);
    nu  = n ./ k;
    lo  = chi2inv_local(0.025, k) ./ k;
    hi  = chi2inv_local(0.975, k) ./ k;
    scatter(ax, t, nu, 18, 'filled', 'MarkerFaceAlpha',0.7, 'MarkerEdgeAlpha',0.7, ...
        'DisplayName','NIS/DOF');
    plot(ax, t, lo, '--', 'LineWidth',1.1, 'Color',[0.2 0.2 0.2], 'DisplayName','95% lower');
    plot(ax, t, hi, '-.', 'LineWidth',1.1, 'Color',[0.2 0.2 0.2], 'DisplayName','95% upper');
    yline(ax, 1, ':', 'Expected', 'Color',[0 0 0], 'LineWidth',1.0, 'HandleVisibility','off');
    xline(ax, dcoDay, ':', 'DCO', 'Color',[0 0 0], 'LineWidth',1.0, 'HandleVisibility','off');
    fracOut = mean((nu < lo) | (nu > hi));
    text(ax, 0.02, 0.98, sprintf('Fraction outside 95%% bounds: %.1f%%%%', 100*fracOut), ...
        'Units','normalized', 'VerticalAlignment','top', 'FontSize',10);
    xlabel(ax,'Days since detection');
    ylabel(ax,'NIS/DOF');
    title(ax, sprintf('NIS consistency (normalized) with 95%% \\chi^2 bounds (%s)', spanStr), 'FontWeight','bold');
    lg = legend(ax,'Location','best'); lg.Box='off';
end
function x = chi2inv_local(p, k)
    x = 2 * gammaincinv(p, k/2);
end
function saveAllFigures(figs, saveDir, prefix, spanStr)
    names = fieldnames(figs);
    for i = 1:numel(names)
        f = figs.(names{i});
        if isempty(f) || ~isgraphics(f); continue; end
        base = sprintf('%sEKF_%s_%s', prefix, spanStr, names{i});
        pngPath = fullfile(saveDir, [base '.png']);
        pdfPath = fullfile(saveDir, [base '.pdf']);
        try
            exportgraphics(f, pngPath, 'Resolution', 300);
            exportgraphics(f, pdfPath, 'ContentType','vector');
        catch
            print(f, pngPath, '-dpng', '-r300');
            print(f, pdfPath, '-dpdf', '-painters');
        end
    end
end