function out = runPrelimPassthruEKF(prelimEkf, outDir)
    if nargin < 2 || isempty(outDir)
        outDir = fullfile(pwd,'figs_Passthru_EKF');
    end
    if ~exist(outDir,'dir'); mkdir(outDir); end
    [params, sc, st, ~, ~] = projectInit();
    [ET0, ~, ~] = initTime();
    params.t0_et       = ET0;
    params.bias_end_et = ET0 + 6*params.tday;
    x0     = prelimEkf.x_final;
    t0     = prelimEkf.tk(end);
    file_6to14 = 'ASTE583_Project_LTB_Measurements_6D-14D_Truth.csv';
    meas = parseMeasurementData(file_6to14, ET0);
    [tk, idx] = sort(meas.t(:));
    stationID = meas.stationID(idx);
    z_rho     = meas.range(idx);
    z_rr      = meas.rr(idx);
    N = numel(tk);
    res_rho = nan(N,1);
    res_rr  = nan(N,1);
    odeOpts = odeset('RelTol',1e-11,'AbsTol',1e-13);
    t_prev  = t0;
    x_prev  = x0;
    for k = 1:N
        t_meas = tk(k);
        if t_meas <= params.bias_end_et
            continue;
        end
        dt = t_meas - t_prev;
        if dt < 0
            continue;
        end
        if dt > 0
            odeFun = @(t,X) projectDynamics(t, X, params, sc);
            [~, Xtraj] = ode113(odeFun, [t_prev t_meas], x_prev, odeOpts);
            x_pred = Xtraj(end,:).';
        else
            x_pred = x_prev;
        end
        x_meas = x_pred;
        x_meas(10) = 0;
        [z_pred, ~] = measurementModel(t_meas, x_meas, stationID(k), params, st);
        if isfinite(z_rho(k))
            res_rho(k) = z_rho(k) - z_pred(1);
        end
        res_rr(k)  = z_rr(k)  - z_pred(2);
        x_prev = x_pred;
        t_prev = t_meas;
    end
    t_days = (tk - ET0)/86400;
    uSt    = unique(stationID);
    colors = lines(numel(uSt));
    figure('Name','Passthru_EKF_RR','Color','w');
    hold on; grid on;
    for i = 1:numel(uSt)
        idx = stationID==uSt(i) & isfinite(res_rr);
        if any(idx)
            scatter(t_days(idx), res_rr(idx)*1e6, 5, colors(i,:), 'filled', ...
                'DisplayName', sprintf('#%d %s', uSt(i), stationName(uSt(i))));
        end
    end
    yline(0,'k-');
    xlabel('Days since detection');
    ylabel('Range-rate residual [mm/s]');
    title('Passthru of prelim 0–6d EKF solution on 6–14d data');
    legend('Location','bestoutside');
    figure('Name','Passthru_EKF_Range','Color','w');
    if any(isfinite(res_rho))
        hold on; grid on;
        for i = 1:numel(uSt)
            idx = stationID==uSt(i) & isfinite(res_rho);
            if any(idx)
                scatter(t_days(idx), res_rho(idx)*1e3, 5, colors(i,:), 'filled', ...
                    'DisplayName', sprintf('#%d %s', uSt(i), stationName(uSt(i))));
            end
        end
        yline(0,'k-');
        xlabel('Days since detection');
        ylabel('Range residual [m]');
        title('Passthru of prelim 0–6d EKF solution on 6–14d data');
        legend('Location','bestoutside');
    else
        axis off;
        text(0.5,0.5,'No range measurements in post-DCO dataset.', ...
             'HorizontalAlignment','center','FontSize',12);
    end
    if exist('saveAllFigures','file') == 2
        saveAllFigures(outDir,'Passthru_EKF');
    end
    out.tk       = tk;
    out.t_days   = t_days;
    out.station  = stationID;
    out.res_rho  = res_rho;
    out.res_rr   = res_rr;
end