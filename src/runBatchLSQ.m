function out = runBatchLSQ(span, mode, opts)
    if nargin < 1 || isempty(span), span = '0to14'; end
    if nargin < 2 || isempty(mode), mode = 'handout'; end
    if nargin < 3, opts = struct(); end
    span = string(span);
    mode = lower(string(mode));
    if ~isfield(opts,'maxIter'),     opts.maxIter    = 8; end
    if ~isfield(opts,'tol_dx'),      opts.tol_dx     = 1e-6; end
    if ~isfield(opts,'useLM'),       opts.useLM      = false; end
    if ~isfield(opts,'lmLambda'),    opts.lmLambda   = 1e-6; end
    if ~isfield(opts,'robust'),      opts.robust     = false; end
    if ~isfield(opts,'huberC'),      opts.huberC     = 3.0; end
    if ~isfield(opts,'earlyMode'),   opts.earlyMode  = "use"; end
    if ~isfield(opts,'earlyDays'),   opts.earlyDays  = 2; end
    if ~isfield(opts,'earlyScale'),  opts.earlyScale = 10; end
    if ~isfield(opts,'makePlots'),   opts.makePlots   = true;  end
    if ~isfield(opts,'saveFigures'), opts.saveFigures = false; end
    if ~isfield(opts,'outDir')
        opts.outDir = fullfile(pwd, sprintf("figs_Batch_
    end
    if ~isfield(opts,'figPrefix')
        opts.figPrefix = sprintf("Batch_
    end
    useLogK     = (mode == "handout_logk" || mode == "ekf_prior_logk");
    freezeK     = (mode == "ekf_k");
    useEkfPrior = (mode == "ekf_k" || mode == "ekf_prior" || mode == "ekf_prior_logk");
    [params, sc, st, X0_phys, P0_phys] = projectInit();
    [ET0, ~, ~] = initTime();
    params.t0_et       = ET0;
    params.bias_end_et = ET0 + 6*params.tday;
    if useEkfPrior
        try
            fprintf('runBatchLSQ: getting k_SRP from EKF (recommend 6to14)...\n');
            ekf   = runPrelimEKF([], '0to6');
            k_ekf = ekf.x_final(7);
            fprintf('  EKF k_SRP = %.6f\n', k_ekf);
            X0_phys(7) = k_ekf;
            P0_phys(7,7) = max(P0_phys(7,7), (0.3)^2);
        catch ME
            warning('runBatchLSQ:EKFfail','EKF k fetch failed: %s. Using handout prior k instead. (%s)', ME.message);
            if freezeK
                freezeK = false;
            end
        end
    end
    if useLogK
        if X0_phys(7) <= 0
            error('runBatchLSQ:BadPriorK','Need positive prior k to use log-k mode.');
        end
        X_apri = X0_phys;
        X_apri(7) = log(X0_phys(7));
        P_apri = P0_phys;
        P_apri(7,7) = P0_phys(7,7) / (X0_phys(7)^2);
    else
        X_apri = X0_phys;
        P_apri = P0_phys;
    end
    [W0, P_apri] = safeWhiten(P_apri);
    files = resolveFiles(span);
    measAll.t = []; measAll.stationID = []; measAll.range = []; measAll.rr = [];
    for f = files
        m = parseMeasurementData(f, ET0);
        measAll.t         = [measAll.t;         m.t(:)];
        measAll.stationID = [measAll.stationID; m.stationID(:)];
        measAll.range     = [measAll.range;     m.range(:)];
        measAll.rr        = [measAll.rr;        m.rr(:)];
    end
    [tMeas, ord] = sort(measAll.t(:));
    stationIDs   = measAll.stationID(ord);
    obs_rho      = measAll.range(ord);
    obs_rr       = measAll.rr(ord);
    N            = numel(tMeas);
    earlyCut = ET0 + opts.earlyDays*86400;
    if opts.earlyMode == "drop"
        keep      = (tMeas > earlyCut);
        tMeas     = tMeas(keep);
        stationIDs= stationIDs(keep);
        obs_rho   = obs_rho(keep);
        obs_rr    = obs_rr(keep);
        N         = numel(tMeas);
    end
    [tU, ~, idxU] = unique(tMeas, 'stable');
    sig.rho_dsn = 1e-3;
    sig.rho_ant = 1e-2;
    sig.rr_dsn  = 1e-7;
    sig.rr_ant  = 1e-6;
    scale.rho_dsn = 5;
    scale.rho_ant = 5;
    scale.rr_dsn  = 5;
    scale.rr_ant  = 5;
    maxIter = opts.maxIter;
    tol_dx  = opts.tol_dx;
    odeOpts = odeset('RelTol',1e-11,'AbsTol',1e-13);
    odeFun = @(t,X) projectDynamics(t, X, params, sc);
    X_est = X_apri;
    iterHist = zeros(maxIter,5);
    fprintf('Batch LSQ using %s (%d measurements), mode=%s\n', strjoin(files,", "), N, mode);
    fprintf('%4s  %12s  %12s  %14s  %14s\n','it','||dX||','cost','RR_RMS [mm/s]','Rho_RMS [m]');
    H_all = []; r_all = []; sigma_all = [];
    for it = 1:maxIter
        X_prop0 = X_est(:);
        if useLogK
            X_prop0(7) = exp(X_est(7));
        end
        X_aug0 = [X_prop0; reshape(eye(10),100,1)];
        if abs(tU(1) - ET0) < 1e-6
            tspan     = tU;
            dropFirst = false;
        else
            tspan     = [ET0; tU];
            dropFirst = true;
        end
        [~, Xaug] = ode113(odeFun, tspan, X_aug0, odeOpts);
        if dropFirst, Xaug = Xaug(2:end,:); end
        XU   = Xaug(:,1:10);
        PhiU = Xaug(:,11:end);
        k0_phys = X_prop0(7);
        nRange = sum(isfinite(obs_rho));
        M      = N + nRange;
        H_all     = zeros(M,10);
        r_all     = zeros(M,1);
        sigma_all = zeros(M,1);
        rowType   = zeros(M,1);
        row = 0;
        for k = 1:N
            ui  = idxU(k);
            xk_phys = XU(ui,:).';
            Phi     = reshape(PhiU(ui,:).',10,10);
            stID = stationIDs(k);
            biasActive = (tMeas(k) <= params.bias_end_et);
            x_meas = xk_phys;
            if ~biasActive
                x_meas(10) = 0;
            end
            [y_pred, Hloc] = measurementModel(tMeas(k), x_meas, stID, params, st);
            if ~biasActive
                Hloc(2,10) = 0;
            end
            isEarly = (tMeas(k) <= earlyCut);
            if isfinite(obs_rho(k))
                row         = row + 1;
                r_all(row)  = obs_rho(k) - y_pred(1);
                H_row = Hloc(1,:) * Phi;
                if useLogK
                    H_row(7) = H_row(7) * k0_phys;
                end
                if freezeK
                    H_row(7) = 0;
                end
                H_all(row,:)   = H_row;
                sigma_all(row) = measSigma(stID,"rho",sig,scale,isEarly,opts);
                rowType(row)   = 1;
            end
            row        = row + 1;
            r_all(row) = obs_rr(k) - y_pred(2);
            H_row = Hloc(2,:) * Phi;
            if useLogK
                H_row(7) = H_row(7) * k0_phys;
            end
            if freezeK
                H_row(7) = 0;
            end
            H_all(row,:)   = H_row;
            sigma_all(row) = measSigma(stID,"rr",sig,scale,isEarly,opts);
            rowType(row)   = 2;
        end
        H_all     = H_all(1:row,:);
        r_all     = r_all(1:row);
        sigma_all = sigma_all(1:row);
        rowType   = rowType(1:row);
        w = 1 ./ sigma_all;
        if opts.robust
            e = r_all .* w;
            w = w .* huberSqrtWeight(e, opts.huberC);
        end
        Hw = H_all .* w;
        rw = r_all .* w;
        Hw_tot = [Hw; W0];
        rw_tot = [rw; W0*(X_apri - X_est)];
        if opts.useLM
            Nmat = Hw_tot.'*Hw_tot;
            g    = Hw_tot.'*rw_tot;
            D    = diag(max(diag(Nmat),1e-12));
            dX   = (Nmat + opts.lmLambda*D) \ g;
        else
            [Q,R] = qr(Hw_tot,0);
            dX    = R \ (Q' * rw_tot);
        end
        X_est = X_est + dX;
        cost = rw_tot.' * rw_tot;
        rr_rms_mmps  = sqrt(mean(r_all(rowType==2).^2)) * 1e6;
        if any(rowType==1)
            rho_rms_m = sqrt(mean(r_all(rowType==1).^2)) * 1e3;
        else
            rho_rms_m = NaN;
        end
        iterHist(it,:) = [it, norm(dX), cost, rr_rms_mmps, rho_rms_m];
        fprintf('%4d  %12.3e  %12.3e  %14.3f  %14.3f\n', it, norm(dX), cost, rr_rms_mmps, rho_rms_m);
        if norm(dX) < tol_dx
            iterHist = iterHist(1:it,:);
            break;
        end
    end
    if opts.useLM
        w = 1 ./ sigma_all;
        if opts.robust
            e = r_all .* w;
            w = w .* huberSqrtWeight(e, opts.huberC);
        end
        Hw   = H_all .* w;
        Nmat = Hw.'*Hw + (W0.'*W0);
        P_int = invSymPD(Nmat);
    else
        Ri    = R \ eye(10);
        P_int = Ri * Ri.';
        P_int = 0.5*(P_int + P_int.');
    end
    X_final = X_est;
    P_final = P_int;
    if useLogK
        eta    = X_est(7);
        k_phys = exp(eta);
        X_final(7) = k_phys;
        J       = eye(10);
        J(7,7)  = k_phys;
        P_final = J * P_int * J.';
        P_final = 0.5*(P_final + P_final.');
    end
    res_prefit  = computeResiduals(X_apri);
    res_postfit = computeResiduals(X_est);
    rr_pre  = res_prefit.rr(isfinite(res_prefit.rr));
    rr_post = res_postfit.rr(isfinite(res_postfit.rr));
    rr_rms_pre  = sqrt(mean(rr_pre.^2));
    rr_rms_post = sqrt(mean(rr_post.^2));
    rho_rms_pre  = NaN;
    rho_rms_post = NaN;
    if any(isfinite(res_prefit.rho))
        rho_pre  = res_prefit.rho(isfinite(res_prefit.rho));
        rho_post = res_postfit.rho(isfinite(res_postfit.rho));
        rho_rms_pre  = sqrt(mean(rho_pre.^2));
        rho_rms_post = sqrt(mean(rho_post.^2));
    end
    fprintf('\nFinal Batch LSQ residual stats (mode=%s):\n', mode);
    fprintf('RR RMS prefit  = %.6g km/s (%.3f mm/s)\n', rr_rms_pre,  rr_rms_pre*1e6);
    fprintf('RR RMS postfit = %.6g km/s (%.3f mm/s)\n', rr_rms_post, rr_rms_post*1e6);
    if ~isnan(rho_rms_pre)
        fprintf('Rho RMS prefit  = %.6g km (%.3f m)\n', rho_rms_pre,  rho_rms_pre*1e3);
        fprintf('Rho RMS postfit = %.6g km (%.3f m)\n', rho_rms_post, rho_rms_post*1e3);
    end
    sig3 = 3*sqrt(max(diag(P_final),0));
    labels = {'x [km]','y [km]','z [km]','vx [km/s]','vy [km/s]','vz [km/s]',...
              'k_SRP','lat_4 [rad]','lon_4 [rad]','bias_rr [km/s]'};
    fprintf('\nFinal Batch LSQ estimate at ET0 (mode=%s):\n', mode);
    for i=1:10
        fprintf('%-12s  %.6e   (± %.6e 3σ)\n', labels{i}, X_final(i), sig3(i));
    end
    if opts.makePlots
        try
            makeBatchPlots();
        catch ME
            warning('runBatchLSQ:PlottingFailed', ...
                'Batch plotting failed: %s', ME.message);
        end
    end
    if opts.saveFigures
        if ~exist(opts.outDir,'dir')
            mkdir(opts.outDir);
        end
        saveAllFigures(opts.outDir, char(opts.figPrefix));
    end
    out.X_est        = X_final;
    out.P_est        = P_final;
    out.iterHist     = iterHist;
    out.files        = files;
    out.ET0          = ET0;
    out.res_prefit   = res_prefit;
    out.res_postfit  = res_postfit;
    out.tMeas        = tMeas;
    out.stationIDs   = stationIDs;
    out.X0           = X0_phys;
    out.P0           = P0_phys;
    out.span         = span;
    out.rr_rms_pre   = rr_rms_pre;
    out.rr_rms_post  = rr_rms_post;
    out.rho_rms_pre  = rho_rms_pre;
    out.rho_rms_post = rho_rms_post;
    function files = resolveFiles(spec)
        spec = string(spec);
        if numel(spec)>1, files = spec; return; end
        s = lower(strtrim(spec));
        if s=="0to6" || s=="0-6"
            files = "ASTE583_Project_LTB_Measurements_0-6D_Truth.csv";
        elseif s=="6to14" || s=="6-14"
            files = "ASTE583_Project_LTB_Measurements_6D-14D_Truth.csv";
        elseif s=="0to14" || s=="0-14" || s=="all"
            files = ["ASTE583_Project_LTB_Measurements_0-6D_Truth.csv", ...
                     "ASTE583_Project_LTB_Measurements_6D-14D_Truth.csv"];
        else
            files = spec;
        end
    end
    function plotBatchIterHistory(iterHist, titleStr, figName)
        it    = iterHist(:,1);
        dX    = iterHist(:,2);
        rrRMS = iterHist(:,4);
        rhoRMS= iterHist(:,5);
        f = figure('Name', figName, 'Color','w');
        tlo = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
        nexttile(tlo);
        semilogy(it, dX, '-o', 'LineWidth', 1.5, 'MarkerSize', 6);
        grid on; box on;
        xlabel('Iteration');
        ylabel('||\Delta X||');
        title([titleStr ' – state update norm']);
        set(gca,'FontSize',12,'LineWidth',1);
        nexttile(tlo);
        hold on;
        plot(it, rrRMS, '-o', 'LineWidth', 1.5, 'MarkerSize', 6, ...
             'DisplayName','RR RMS [mm/s]');
        if any(~isnan(rhoRMS))
            plot(it, rhoRMS, '-s', 'LineWidth', 1.5, 'MarkerSize', 6, ...
                 'DisplayName','Range RMS [m]');
        end
        grid on; box on;
        xlabel('Iteration');
        ylabel('RMS');
        title([titleStr ' – residual RMS']);
        legend('Location','best');
        set(gca,'FontSize',12,'LineWidth',1);
    end
    function res = computeResiduals(X_epoch_int)
        X0prop = X_epoch_int(:);
        if useLogK
            X0prop(7) = exp(X_epoch_int(7));
        end
        X_aug0 = [X0prop; reshape(eye(10),100,1)];
        if abs(tU(1)-ET0)<1e-6
            tspanR     = tU;
            dropFirstR = false;
        else
            tspanR     = [ET0; tU];
            dropFirstR = true;
        end
        [~, XaugR] = ode113(odeFun, tspanR, X_aug0, odeOpts);
        if dropFirstR, XaugR = XaugR(2:end,:); end
        XU_r = XaugR(:,1:10);
        res.rho = nan(N,1);
        res.rr  = nan(N,1);
        for kk=1:N
            xk = XU_r(idxU(kk),:).';
            biasActive = (tMeas(kk) <= params.bias_end_et);
            if ~biasActive
                xk(10) = 0;
            end
            y_pred = measurementModel(tMeas(kk), xk, stationIDs(kk), params, st);
            if isfinite(obs_rho(kk))
                res.rho(kk) = obs_rho(kk) - y_pred(1);
            end
            res.rr(kk) = obs_rr(kk) - y_pred(2);
        end
    end
    function makeBatchPlots()
    spanChar = char(span);
    modeChar = char(mode);
    tDays = (tMeas - ET0) / 86400;
    dcoDay = inferDcoDay(spanChar);
    ltmDay = 15;
    earlyDay = opts.earlyDays;
    uSt = unique(stationIDs(:));
    [stLbl, stClr, stMrk] = stationStyles(uSt);
    rr_pre  = res_prefit.rr(:)  * 1e6;
    rr_post = res_postfit.rr(:) * 1e6;
    sigRR_DSN  = 0.1;
    sigRR_ANT  = 1.0;
    plotResidualStack( ...
        sprintf('Batch_%s_%s_RR', spanChar, modeChar), ...
        sprintf('Batch %s (%s): Range-rate residuals', spanChar, modeChar), ...
        tDays, stationIDs, rr_pre, rr_post, ...
        'Range-rate residual [mm/s]', ...
        dcoDay, earlyDay, ...
        stLbl, stClr, stMrk, ...
        sigRR_DSN, sigRR_ANT);
    hasRange = any(isfinite(res_prefit.rho(:))) || any(isfinite(res_postfit.rho(:)));
    if hasRange
        rho_pre  = res_prefit.rho(:)  * 1e3;
        rho_post = res_postfit.rho(:) * 1e3;
        sigR_DSN = 1.0;
        sigR_ANT = 10.0;
        plotResidualStack( ...
            sprintf('Batch_%s_%s_Range', spanChar, modeChar), ...
            sprintf('Batch %s (%s): Range residuals', spanChar, modeChar), ...
            tDays, stationIDs, rho_pre, rho_post, ...
            'Range residual [m]', ...
            dcoDay, earlyDay, ...
            stLbl, stClr, stMrk, ...
            sigR_DSN, sigR_ANT);
    else
        f = figure('Name', sprintf('Batch_%s_%s_Range', spanChar, modeChar), 'Color','w');
        ax = axes(f); axis(ax,'off');
        text(ax, 0.5, 0.5, 'No range measurements in this dataset.', ...
            'HorizontalAlignment','center', 'FontSize',12);
    end
    plotIterHistoryBetter(iterHist, ...
        sprintf('Batch %s (%s)', spanChar, modeChar), ...
        sprintf('Batch_%s_%s_IterHist', spanChar, modeChar));
    plotUncertaintyAndRequirement( ...
        sprintf('Batch_%s_%s_Uncertainty', spanChar, modeChar), ...
        sprintf('Batch %s (%s): Uncertainty diagnostics', spanChar, modeChar), ...
        X_final, P_final, ...
        dcoDay, ltmDay);
    plotParamSummary( ...
        sprintf('Batch_%s_%s_Params', spanChar, modeChar), ...
        sprintf('Batch %s (%s): Estimated parameters @ ET0', spanChar, modeChar), ...
        X_final, P_final, X0_phys);
    function dco = inferDcoDay(spanText)
        s = lower(string(spanText));
        if contains(s,"0to6") || contains(s,"0-6")
            dco = 6;
        else
            dco = 14;
        end
    end
    function [labels, colors, markers] = stationStyles(uStLocal)
        colorsMap = containers.Map('KeyType','double','ValueType','any');
        colorsMap(1) = [0.0000 0.4470 0.7410];
        colorsMap(2) = [0.8500 0.3250 0.0980];
        colorsMap(3) = [0.9290 0.6940 0.1250];
        colorsMap(4) = [0.4940 0.1840 0.5560];
        labels  = cell(size(uStLocal));
        colors  = zeros(numel(uStLocal),3);
        markers = cell(size(uStLocal));
        for iS = 1:numel(uStLocal)
            id = uStLocal(iS);
            labels{iS} = sprintf('#%d %s', id, stationName(id));
            if isKey(colorsMap,id)
                colors(iS,:) = colorsMap(id);
            else
                colors(iS,:) = lines(1);
            end
            if id == 4
                markers{iS} = 's';
            else
                markers{iS} = 'o';
            end
        end
    end
    function plotResidualStack(figName, superTitle, t, stID, yPre, yPost, yLabel, dco, early, labels, colors, markers, sigDSN, sigANT)
        f = figure('Name', figName, 'Color','w');
        tl = tiledlayout(f, 2, 1, 'TileSpacing','compact', 'Padding','compact');
        ax1 = nexttile(tl, 1);
        plotResidualTile(ax1, t, stID, yPre, yLabel, 'Prefit', dco, early, labels, colors, markers, sigDSN, sigANT, true);
        ax2 = nexttile(tl, 2);
        plotResidualTile(ax2, t, stID, yPost, yLabel, 'Postfit', dco, early, labels, colors, markers, sigDSN, sigANT, false);
        linkaxes([ax1 ax2], 'x');
        sgtitle(tl, superTitle, 'FontWeight','bold');
    end
    function plotResidualTile(ax, t, stID, y, yLabelLocal, whichStr, dco, early, labels, colors, markers, sigDSN, sigANT, showLegend)
        hold(ax,'on');
        grid(ax,'on'); grid(ax,'minor');
        ax.Box = 'on';
        ax.Layer = 'top';
        if isfinite(early) && early > min(t) && early < max(t)
            ylTmp = [-1 1];
            patch(ax, [min(t) early early min(t)], [ylTmp(1) ylTmp(1) ylTmp(2) ylTmp(2)], [0 0 0], ...
                'FaceAlpha', 0.03, 'EdgeColor','none', 'HandleVisibility','off');
        end
        xline(ax, dco, ':', 'DCO', 'Color',[0 0 0], 'LineWidth',1.0, ...
        'LabelVerticalAlignment','middle', 'LabelHorizontalAlignment','left', ...
        'HandleVisibility','off');
        yline(ax, 0, '-', 'Color',[0 0 0], 'LineWidth',0.9, 'HandleVisibility','off');
        yline(ax, +3*sigDSN, '--', 'Color',[0.25 0.25 0.25], 'LineWidth',0.9, 'HandleVisibility','off');
        yline(ax, -3*sigDSN, '--', 'Color',[0.25 0.25 0.25], 'LineWidth',0.9, 'HandleVisibility','off');
        yline(ax, +3*sigANT,  ':', 'Color',[0.25 0.25 0.25], 'LineWidth',1.0, 'HandleVisibility','off');
        yline(ax, -3*sigANT,  ':', 'Color',[0.25 0.25 0.25], 'LineWidth',1.0, 'HandleVisibility','off');
        maxPtsPerStation = 25000;
        for iS = 1:numel(uSt)
            id = uSt(iS);
            idx = (stID == id) & isfinite(t) & isfinite(y);
            if ~any(idx), continue; end
            tt = t(idx);
            yy = y(idx);
            if numel(yy) > maxPtsPerStation
                stride = ceil(numel(yy)/maxPtsPerStation);
                tt = tt(1:stride:end);
                yy = yy(1:stride:end);
            end
            scatter(ax, tt, yy, 14, ...
                'Marker', markers{iS}, ...
                'MarkerFaceColor', colors(iS,:), ...
                'MarkerEdgeColor', colors(iS,:), ...
                'MarkerFaceAlpha', 0.55, ...
                'MarkerEdgeAlpha', 0.55, ...
                'DisplayName', labels{iS});
        end
        yl = robustSymLim(y, max(3*sigDSN, 3*sigANT));
        ylim(ax, 1.15*yl);
        ch = ax.Children;
        if ~isempty(ch) && isa(ch(end),'matlab.graphics.primitive.Patch')
            p = ch(end);
            p.YData = [ylim(ax) fliplr(ylim(ax))];
        end
        title(ax, sprintf('%s %s', whichStr, whichMeasurementLabel(yLabelLocal)), 'FontWeight','bold');
        xlabel(ax, 'Days since detection');
        ylabel(ax, yLabelLocal);
        yyAll = y(isfinite(y));
        if ~isempty(yyAll)
            rmsVal = sqrt(mean(yyAll.^2));
            text(ax, 0.01, 0.98, sprintf('RMS(%s): %.3g', whichStr, rmsVal), ...
                'Units','normalized', 'VerticalAlignment','top', 'FontSize',10);
        end
        text(ax, 0.01, 0.90, sprintf('\\pm3\\sigma guides: DSN %.3g, Ant %.3g', 3*sigDSN, 3*sigANT), ...
            'Units','normalized', 'VerticalAlignment','top', 'FontSize',9, 'Color',[0.2 0.2 0.2]);
        if showLegend
            hScat = findobj(ax, 'Type', 'Scatter');
            hScat = hScat(~cellfun(@isempty, get(hScat,'DisplayName')));
            if ~isempty(hScat)
                if ~iscell(hScat)
                    hScat = num2cell(hScat);
                end
                names = cellfun(@(h) get(h,'DisplayName'), hScat, 'UniformOutput', false);
                [namesUnique, ia] = unique(names, 'stable');
                hUnique = [hScat{ia}];
                lg = legend(ax, hUnique, namesUnique, 'Location','eastoutside');
                lg.Box = 'off';
            end
        end
    end
    function s = whichMeasurementLabel(yLabelLocal)
        if contains(lower(yLabelLocal),'rate')
            s = 'RR residuals';
        else
            s = 'Range residuals';
        end
    end
    function yl = robustSymLim(y, mustInclude)
        y = y(isfinite(y));
        if isempty(y)
            m = mustInclude;
        else
            a = sort(abs(y));
            k = max(1, round(0.995*numel(a)));
            m = a(k);
            m = max(m, mustInclude);
            if m == 0, m = mustInclude; end
        end
        yl = [-m m];
    end
    function plotIterHistoryBetter(iterH, titleStr, figName)
        it    = iterH(:,1);
        dX    = iterH(:,2);
        cost  = iterH(:,3);
        rrRMS = iterH(:,4);
        rhoRMS= iterH(:,5);
        f = figure('Name', figName, 'Color','w');
        tl = tiledlayout(f, 3, 1, 'TileSpacing','compact', 'Padding','compact');
        ax1 = nexttile(tl,1);
        semilogy(ax1, it, max(dX, eps), '-o', 'LineWidth',1.4, 'MarkerSize',6);
        grid(ax1,'on'); grid(ax1,'minor'); ax1.Box='on';
        ylabel(ax1,'||\Delta X||');
        title(ax1,'State update norm');
        ax2 = nexttile(tl,2);
        semilogy(ax2, it, max(cost, eps), '-o', 'LineWidth',1.4, 'MarkerSize',6);
        grid(ax2,'on'); grid(ax2,'minor'); ax2.Box='on';
        ylabel(ax2,'Cost');
        title(ax2,'Weighted cost');
        ax3 = nexttile(tl,3);
        hold(ax3,'on');
        plot(ax3, it, rrRMS, '-o', 'LineWidth',1.4, 'MarkerSize',6, 'DisplayName','RR RMS [mm/s]');
        if any(isfinite(rhoRMS))
            plot(ax3, it, rhoRMS, '-s', 'LineWidth',1.4, 'MarkerSize',6, 'DisplayName','Range RMS [m]');
        end
        grid(ax3,'on'); grid(ax3,'minor'); ax3.Box='on';
        xlabel(ax3,'Iteration'); ylabel(ax3,'RMS');
        title(ax3,'Residual RMS');
        legend(ax3,'Location','best'); legend(ax3,'boxoff');
        sgtitle(tl, titleStr, 'FontWeight','bold');
    end
    function plotUncertaintyAndRequirement(figName, superTitle, X0, P0, dco, ltm)
        maxDay = max([ltm, max(tDays(isfinite(tDays)))]);
        nPlot  = 250;
        tPlotDays = unique([linspace(0, maxDay, nPlot), dco, ltm]).';
        tPlotET   = ET0 + tPlotDays*86400;
        odeOptsPlot = odeset(odeOpts, 'RelTol',1e-10, 'AbsTol',1e-12);
        X_aug0 = [X0(:); reshape(eye(10),100,1)];
        [~, Xaug] = ode113(odeFun, tPlotET, X_aug0, odeOptsPlot);
        Xs = Xaug(:,1:10);
        R  = Xs(:,1:3);
        V  = Xs(:,4:6);
        [~, iD] = min(abs(tPlotDays - dco));
        rD = R(iD,:).';
        vD = V(iD,:).';
        Ad = rtnA(rD, vD);
        dR = (R - rD.');          % n x 3
        prtn = (Ad * dR.').';
        pR = prtn(:,1);
        pT = prtn(:,2);
        sigmaR = nan(size(tPlotDays));
        sigmaT = nan(size(tPlotDays));
        sigmaN = nan(size(tPlotDays));
        Ppos_all = zeros(3,3,numel(tPlotDays));
        for k = 1:numel(tPlotDays)
            Phi = reshape(Xaug(k,11:end), 10, 10);
            Pt  = Phi * P0 * Phi.';
            Pt  = 0.5*(Pt + Pt.');
            Ppos = Pt(1:3,1:3);
            Ppos_all(:,:,k) = Ppos;
            Ak = rtnA(R(k,:).', V(k,:).');
            P_rtn = Ak * Ppos * Ak.';
            s = sqrt(max(diag(P_rtn), 0));
            sigmaR(k) = s(1);
            sigmaT(k) = s(2);
            sigmaN(k) = s(3);
        end
        [~, iL] = min(abs(tPlotDays - ltm));
        sigT_L = sigmaT(iL);
        f = figure('Name', figName, 'Color','w');
        tl = tiledlayout(f, 1, 2, 'TileSpacing','compact', 'Padding','compact');
        ax1 = nexttile(tl,1);
        hold(ax1,'on'); grid(ax1,'on'); grid(ax1,'minor'); ax1.Box='on'; ax1.Layer='top';
        plot(ax1, pR, pT, '-', 'Color',[0 0 0], 'LineWidth',1.2, 'DisplayName','Mean (relative)');
        scatter(ax1, 0, 0, 35, 'filled', 'MarkerFaceColor',[0 0 0], 'DisplayName','DCO');
        idxE = unique(round(linspace(1, numel(tPlotDays), 6)));
        idxE = unique([idxE(:); iD; iL]);
        for ii = 1:numel(idxE)
            k = idxE(ii);
            mu = [pR(k); pT(k)];
            P_rtn_dco = Ad * Ppos_all(:,:,k) * Ad.';
            P2 = P_rtn_dco(1:2,1:2);
            a = 0.08 + 0.18*(ii/numel(idxE));
            covEllipsePatch(ax1, mu, P2, 3, [0.1 0.4 0.8], a);
            if k == iL
                scatter(ax1, mu(1), mu(2), 28, 'filled', 'MarkerFaceColor',[0.1 0.4 0.8], ...
                    'DisplayName','LTM (relative)');
            end
        end
        axis(ax1,'equal');
        xlabel(ax1,'\DeltaR [km]');
        ylabel(ax1,'\DeltaT (along-track) [km]');
        title(ax1,'3\sigma covariance in DCO RT plane', 'FontWeight','bold');
        legend(ax1,'Location','best'); legend(ax1,'boxoff');
        ax2 = nexttile(tl,2);
        hold(ax2,'on'); grid(ax2,'on'); grid(ax2,'minor'); ax2.Box='on'; ax2.Layer='top';
        plot(ax2, tPlotDays, sigmaR, '-', 'LineWidth',1.3, 'DisplayName','\sigma_R');
        plot(ax2, tPlotDays, sigmaT, '-', 'LineWidth',1.3, 'DisplayName','\sigma_T (along-track)');
        plot(ax2, tPlotDays, sigmaN, '-', 'LineWidth',1.3, 'DisplayName','\sigma_N');
        xline(ax2, dco, ':', 'DCO', 'Color',[0 0 0], 'LineWidth',1.0);
        xline(ax2, ltm, '--', 'LTM', 'Color',[0 0 0], 'LineWidth',1.0);
        yline(ax2, 100, '--', '100 km req', 'Color',[0.2 0.2 0.2], 'LineWidth',1.0);
        text(ax2, 0.02, 0.98, sprintf('\\sigma_T @ LTM (t=%.0f d): %.2f km', ltm, sigT_L), ...
            'Units','normalized', 'VerticalAlignment','top', 'FontSize',10);
        xlabel(ax2,'Days since detection');
        ylabel(ax2,'1\sigma position uncertainty [km]');
        title(ax2,'RTN uncertainties + along-track requirement', 'FontWeight','bold');
        legend(ax2,'Location','best'); legend(ax2,'boxoff');
        sgtitle(tl, superTitle, 'FontWeight','bold');
        if sigT_L < 100
            fprintf('Batch LSQ: ✓ Meets LTM along-track requirement: sigmaT = %.2f km < 100 km\n', sigT_L);
        else
            fprintf('Batch LSQ: ✗ DOES NOT meet LTM along-track requirement: sigmaT = %.2f km >= 100 km\n', sigT_L);
        end
    end
    function plotParamSummary(figName, superTitle, Xhat, Phat, XpriorPhys)
        k_est   = Xhat(7);
        lat_deg = rad2deg(Xhat(8));
        lon_deg = rad2deg(Xhat(9));
        b_mmps  = Xhat(10) * 1e6;
        sig3 = 3*sqrt(max(diag(Phat),0));
        k_3   = sig3(7);
        lat_3 = rad2deg(sig3(8));
        lon_3 = rad2deg(sig3(9));
        b_3   = sig3(10) * 1e6;
        k0   = XpriorPhys(7);
        lat0 = rad2deg(XpriorPhys(8));
        lon0 = rad2deg(XpriorPhys(9));
        b0   = XpriorPhys(10) * 1e6;
        f = figure('Name', figName, 'Color','w');
        ax = axes(f);
        hold(ax,'on'); grid(ax,'on'); grid(ax,'minor'); ax.Box='on'; ax.Layer='top';
        x = 1:4;
        y = [k_est, lat_deg, lon_deg, b_mmps];
        e = [k_3,   lat_3,   lon_3,   b_3];
        errorbar(ax, x, y, e, 'o', 'LineWidth',1.4, 'MarkerSize',7, ...
            'CapSize',10, 'DisplayName','Estimate \pm 3\sigma');
        scatter(ax, x, [k0, lat0, lon0, b0], 45, 'x', 'LineWidth',1.6, ...
            'DisplayName','Prior (mean)');
        set(ax,'XTick',x, 'XTickLabel', {'k_{SRP}','lat_4 [deg]','lon_4 [deg]','RR bias [mm/s]'}, ...
            'XTickLabelRotation',20);
        title(ax, superTitle, 'FontWeight','bold');
        ylabel(ax,'Value');
        legend(ax,'Location','best'); legend(ax,'boxoff');
    end
    function A = rtnA(r, v)
        r = r(:); v = v(:);
        Rhat = r / max(norm(r), eps);
        h = cross(r, v);
        if norm(h) < 1e-12
            tmp = [1;0;0];
            if abs(dot(tmp,Rhat)) > 0.9
                tmp = [0;1;0];
            end
            Nhat = cross(Rhat, tmp);
            Nhat = Nhat / max(norm(Nhat), eps);
        else
            Nhat = h / norm(h);
        end
        That = cross(Nhat, Rhat);
        That = That / max(norm(That), eps);
        A = [Rhat.'; That.'; Nhat.'];
    end
    function covEllipsePatch(ax, mu, P2, nsig, rgb, faceAlpha)
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
    end
end
function [W0, Pspd] = safeWhiten(P)
    Pspd  = 0.5*(P+P.');
    n     = size(Pspd,1);
    jitter= 0;
    for k=1:8
        [L,p] = chol(Pspd,'lower');
        if p==0
            W0 = L \ eye(n);
            return;
        end
        jitter = max(1e-18, 10^k * 1e-18);
        Pspd   = Pspd + jitter*eye(n);
    end
    error('safeWhiten:SPD','Could not Cholesky-factor prior covariance.');
end
function s = measSigma(stID, type, sig, scale, isEarly, opts)
    if stID == 4
        if type=="rho", s = sig.rho_ant*scale.rho_ant;
        else,           s = sig.rr_ant *scale.rr_ant;
        end
    else
        if type=="rho", s = sig.rho_dsn*scale.rho_dsn;
        else,           s = sig.rr_dsn *scale.rr_dsn;
        end
    end
    if opts.earlyMode == "downweight" && isEarly
        s = s * opts.earlyScale;
    end
end
function w = huberSqrtWeight(e, c)
    a = abs(e);
    w = ones(size(e));
    ii = a > c;
    w(ii) = sqrt(c ./ a(ii));
end
function P = invSymPD(N)
    N = 0.5*(N+N.');
    R = chol(N,'lower');
    Ri= R \ eye(size(N,1));
    P = Ri.' * Ri;
    P = 0.5*(P+P.');
end
function labels = stationLegend(uSt)
    labels = cell(1,numel(uSt));
    for i=1:numel(uSt)
        labels{i} = sprintf('#%d %s', uSt(i), stationName(uSt(i)));
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
function plotCovEllipse(mu, P2, nsig)
    P2 = 0.5*(P2+P2.');
    [V,D] = eig(P2);
    D     = max(D, 0);
    A     = V*sqrt(D);
    th   = linspace(0, 2*pi, 200);
    circ = [cos(th); sin(th)];
    ell  = mu(:) + nsig * (A * circ);
    plot(ell(1,:), ell(2,:), 'k--', 'LineWidth',1.2, ...
         'DisplayName', sprintf('%g\\sigma ellipse', nsig));
    scatter(mu(1), mu(2), 25, 'k', 'filled');
end
function saveAllFigures(outDir, prefix)
    figs = findall(0, 'Type', 'figure');
    if isempty(figs)
        fprintf('No figures to save.\n');
        return;
    end
    nums = arrayfun(@(f) f.Number, figs);
    [~,ord] = sort(nums);
    figs = figs(ord);
    for i = 1:numel(figs)
        f = figs(i);
        if isempty(f.Name)
            f.Name = sprintf('Figure%d', f.Number);
        end
        safeName = regexprep(f.Name, '[^a-zA-Z0-9_-]', '_');
        base = sprintf('%s_%02d_%s', prefix, f.Number, safeName);
        pngPath = fullfile(outDir, [base '.png']);
        pdfPath = fullfile(outDir, [base '.pdf']);
        figPath = fullfile(outDir, [base '.fig']);
        set(f, 'Color', 'w');
        if exist('exportgraphics','file') == 2
            exportgraphics(f, pngPath, 'Resolution', 300);
            exportgraphics(f, pdfPath, 'ContentType', 'vector');
        else
            print(f, fullfile(outDir, base), '-dpng', '-r300');
            print(f, fullfile(outDir, base), '-dpdf', '-r300');
        end
        if exist('savefig','file') == 2
            savefig(f, figPath);
        end
    end
    fprintf('Saved %d figures to: %s\n', numel(figs), outDir);
end