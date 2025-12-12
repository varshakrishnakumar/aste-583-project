function out = runPrelimPassthruBatch(batch06, span, opts)
    if nargin < 2 || isempty(span), span = '6to14'; end
    if nargin < 3, opts = struct(); end
    if ~isfield(opts,'makePlots'), opts.makePlots = true; end
    span = lower(string(span));
    [params, sc, st, ~, ~] = projectInit();
    [ET0, ~, ~]            = initTime();
    params.t0_et           = ET0;
    params.bias_end_et     = ET0 + 6*params.tday;
    switch span
        case {"6to14","6-14"}
            files = "ASTE583_Project_LTB_Measurements_6D-14D_Truth.csv";
        otherwise
            files = string(span);
    end
    meas       = parseMeasurementData(files, ET0);
    tMeas      = meas.t(:);
    stationIDs = meas.stationID(:);
    obs_rho    = meas.range(:);
    obs_rr     = meas.rr(:);
    N          = numel(tMeas);
    [tU, ~, idxU] = unique(tMeas, 'stable');
    odeFun  = @(t,X) projectDynamics(t, X, params, sc);
    odeOpts = odeset('RelTol',1e-11,'AbsTol',1e-13);
    X0prop  = batch06.X_est(:);
    X_aug0  = [X0prop; reshape(eye(10),100,1)];
    if abs(tU(1) - ET0) < 1e-6
        tspan     = tU;
        dropFirst = false;
    else
        tspan     = [ET0; tU];
        dropFirst = true;
    end
    [~, Xaug] = ode113(odeFun, tspan, X_aug0, odeOpts);
    if dropFirst
        Xaug = Xaug(2:end,:);
    end
    XU = Xaug(:,1:10);
    res_rho = NaN(N,1);
    res_rr  = NaN(N,1);
    for k = 1:N
        ui = idxU(k);
        xk = XU(ui,:).';                      % 10x1 at this time
        biasActive = (tMeas(k) <= params.bias_end_et);
        if ~biasActive
            xk(10) = 0;
        end
        y_pred = measurementModel(tMeas(k), xk, stationIDs(k), params, st);
        if isfinite(obs_rho(k))
            res_rho(k) = obs_rho(k) - y_pred(1);
        end
        if isfinite(obs_rr(k))
            res_rr(k)  = obs_rr(k)  - y_pred(2);
        end
    end
    out.tMeas      = tMeas;
    out.stationIDs = stationIDs;
    out.res_rho    = res_rho;
    out.res_rr     = res_rr;
    out.ET0        = ET0;
    out.files      = files;
    if opts.makePlots
        makePassthruPlots(out);
    end
end
function makePassthruPlots(out)
    ET0        = out.ET0;
    tDays      = (out.tMeas - ET0) / 86400;
    stationIDs = out.stationIDs;
    res_rho    = out.res_rho;
    res_rr     = out.res_rr;
    uSt    = unique(stationIDs);
    colors = lines(numel(uSt));
    figure('Name','Batch_passthru_RR','Color','w'); hold on;
    for i = 1:numel(uSt)
        idx = (stationIDs == uSt(i)) & isfinite(res_rr);
        if any(idx)
            plot(tDays(idx), res_rr(idx)*1e6, '.', ...
                'Color', colors(i,:), ...
                'DisplayName', stationName(uSt(i)));
        end
    end
    yline(0,'k-','LineWidth',1);
    grid on; box on;
    xlabel('Time since detection [days]');
    ylabel('Pass-through RR residual [mm/s]');
    title('Pass-through RR residuals: 0–6 batch fit vs 6–14 data');
    legend('Location','bestoutside');
    if any(isfinite(res_rho))
        figure('Name','Batch_passthru_Range','Color','w'); hold on;
        for i = 1:numel(uSt)
            idx = (stationIDs == uSt(i)) & isfinite(res_rho);
            if any(idx)
                plot(tDays(idx), res_rho(idx)*1e3, '.', ...
                    'Color', colors(i,:), ...
                    'DisplayName', stationName(uSt(i)));
            end
        end
        yline(0,'k-','LineWidth',1);
        grid on; box on;
        xlabel('Time since detection [days]');
        ylabel('Pass-through range residual [m]');
        title('Pass-through range residuals: 0–6 batch fit vs 6–14 data');
        legend('Location','bestoutside');
    end
end