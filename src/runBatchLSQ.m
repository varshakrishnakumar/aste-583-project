function out = runBatchLSQ(dataSpec)
% runBatchLSQ  Batch Gauss-Newton least squares OD (supports 0–6, 6–14, 0–14).
%
% Estimates 10-state at detection epoch ET0:
%   [r; v; k_SRP; lat4; lon4; bias_rr]
%%  - Propagate from ET0 (not first measurement time)
%  - Proper RR bias gating (active only first 6 days)
%  - Station-dependent measurement sigmas (range + rr)
%  - A-priori regularization using P0 (prevents singular normal matrix)
%  - No giant diag(W) allocation
%  - Handles file(s) or span keywords ('0to6','6to14','0to14')
%
% Output:
%   out struct with final state, covariance, iteration history, residual stats.

    % ---------- Inputs ----------
    if nargin < 1 || isempty(dataSpec)
        dataSpec = "0to14";
    end
    dataSpec = string(dataSpec);

    % ---------- Init ----------
    [params, sc, st, X0, P0] = projectInit();
    [ET0, ~, ~] = initTime();

    params.t0_et = ET0;                 % needed for station rotation
    params.bias_end_et = ET0 + 6*params.tday;

    % ---------- Load measurements ----------
    files = resolveFiles(dataSpec);

    measAll.t = []; measAll.stationID = []; measAll.range = []; measAll.rr = [];
    for f = files
        m = parseMeasurementData(f, ET0);
        measAll.t         = [measAll.t;         m.t(:)];
        measAll.stationID = [measAll.stationID; m.stationID(:)];
        measAll.range     = [measAll.range;     m.range(:)];
        measAll.rr        = [measAll.rr;        m.rr(:)];
    end

    % Sort by time (critical)
    [tMeas, ord] = sort(measAll.t(:));
    stationIDs = measAll.stationID(ord);
    obs_rho    = measAll.range(ord);
    obs_rr     = measAll.rr(ord);

    N = numel(tMeas);

    % Unique times for efficient propagation (handles duplicate times cleanly)
    [tU, ~, idxU] = unique(tMeas, 'stable');

    % ---------- Measurement sigmas (project spec) ----------
    sig.rho_dsn = 1e-3;   % km  (1 m)
    sig.rho_ant = 1e-2;   % km  (10 m)
    sig.rr_dsn  = 1e-7;   % km/s (0.1 mm/s)
    sig.rr_ant  = 1e-6;   % km/s (1.0 mm/s)

    % Optional inflation knobs (use if residuals/NIS are huge)
    scale.rho_dsn = 1;
    scale.rho_ant = 1;
    scale.rr_dsn  = 1;
    scale.rr_ant  = 1;

    % ---------- Iteration controls ----------
    maxIter = 8;
    tol_dx  = 1e-6;   % absolute norm; you can tighten later
    odeOpts = odeset('RelTol',1e-11,'AbsTol',1e-13);

    % A-priori regularization (HIGHLY recommended)
    useApriori = true;
    X_apri = X0;
    P_apri = P0;
    W0 = chol(P_apri, 'lower') \ eye(10);   % P0^{-1/2}

    % ---------- Solve loop ----------
    X_est = X0;

    iterHist = zeros(maxIter, 5); % [iter normDx cost rrRMS rhoRMS]
    fprintf('Batch LSQ using %s (%d meas)\n', strjoin(files, ", "), N);
    fprintf('%4s  %12s  %12s  %12s  %12s\n', 'it', '||dX||', 'cost', 'RR_RMS(mm/s)', 'RHO_RMS(m)');

    for it = 1:maxIter
        % ---- Propagate state+STM once to all unique times ----
        X_aug0 = [X_est(:); reshape(eye(10),100,1)];
        tspan  = [ET0; tU];  % ensures Phi is from ET0
        odeFun = @(t,X) projectDynamics(t, X, params, sc);

        [~, Xaug] = ode113(odeFun, tspan, X_aug0, odeOpts);
        Xaug = Xaug(2:end,:); % rows aligned with tU

        XU   = Xaug(:,1:10);              % |tU| x 10
        PhiU = Xaug(:,11:end);            % |tU| x 100

        % ---- Preallocate rows ----
        nRange = sum(isfinite(obs_rho));
        M = N + nRange; % rr always + range when present

        H_all   = zeros(M,10);
        res_all = zeros(M,1);
        w       = zeros(M,1); % 1/sigma
        rowType = zeros(M,1); % 1=range, 2=rr (optional debug)
        rowSta  = zeros(M,1);

        r = 0;

        % ---- Build residuals + design matrix ----
        for k = 1:N
            tu_idx = idxU(k);

            xk   = XU(tu_idx,:).';
            Phi  = reshape(PhiU(tu_idx,:).', 10, 10);

            stID = stationIDs(k);

            biasActive = (tMeas(k) <= params.bias_end_et);

            x_for_meas = xk;
            if ~biasActive
                x_for_meas(10) = 0; % force bias off
            end

            [y_pred, Hloc] = measurementModel(tMeas(k), x_for_meas, stID, params, st);

            % if bias off, force sensitivity off too (robust even if measurementModel always has H(2,10)=1)
            if ~biasActive
                Hloc(2,10) = 0;
            end

            % RANGE (if present)
            if isfinite(obs_rho(k))
                r = r + 1;
                res_all(r) = obs_rho(k) - y_pred(1);
                H_all(r,:) = Hloc(1,:) * Phi;

                w(r) = 1 / measSigma(stID, "rho", sig, scale);
                rowType(r) = 1; rowSta(r) = stID;
            end

            % RANGE-RATE
            r = r + 1;
            res_all(r) = obs_rr(k) - y_pred(2);
            H_all(r,:) = Hloc(2,:) * Phi;

            w(r) = 1 / measSigma(stID, "rr", sig, scale);
            rowType(r) = 2; rowSta(r) = stID;
        end

        % ---- Apply weights (no diag(W)!) ----
        H_w   = H_all .* w;
        res_w = res_all .* w;

        % ---- Add a-priori ----
        if useApriori
            H_w   = [H_w;   W0];
            res_w = [res_w; W0 * (X_apri - X_est)];
        end

        % ---- Solve weighted LSQ via QR ----
        [Q,R] = qr(H_w, 0);
        dX = R \ (Q' * res_w);

        % ---- Update ----
        X_new = X_est + dX;

        % ---- Diagnostics ----
        cost = res_w.' * res_w;

        rr_res  = res_all(rowType==2);
        rho_res = res_all(rowType==1);

        rr_rms_mmps  = sqrt(mean(rr_res.^2)) * 1e6;  % km/s -> mm/s
        if isempty(rho_res)
            rho_rms_m = NaN;
        else
            rho_rms_m = sqrt(mean(rho_res.^2)) * 1e3; % km -> m
        end

        iterHist(it,:) = [it, norm(dX), cost, rr_rms_mmps, rho_rms_m];
        fprintf('%4d  %12.3e  %12.3e  %12.3f  %12.3f\n', it, norm(dX), cost, rr_rms_mmps, rho_rms_m);

        X_est = X_new;

        if norm(dX) < tol_dx
            iterHist = iterHist(1:it,:);
            break;
        end
    end

    % ---------- Final covariance estimate ----------
    % At final iteration, normal matrix approx N = R'R (from QR of H_w)
    % Cov ≈ inv(N) = inv(R) * inv(R)'
    Ri = R \ eye(10);
    P_est = Ri * Ri.';
    P_est = 0.5*(P_est + P_est.');

    % ---------- Prefit vs Postfit residual plots ----------
    res_prefit = computeResiduals(X0);
    res_postfit = computeResiduals(X_est);

    t_days = (tMeas - ET0) / 86400;

    figure('Color','w','Name','Batch LSQ residuals');
    if any(isfinite(obs_rho))
        subplot(2,1,1); hold on; grid on;
        plot(t_days, res_prefit.rho*1e3, '.', 'DisplayName','prefit');  % m
        plot(t_days, res_postfit.rho*1e3, '.', 'DisplayName','postfit');
        ylabel('Range resid [m]'); legend('Location','best');
        title('Range residuals');

        subplot(2,1,2); hold on; grid on;
        plot(t_days, res_prefit.rr*1e6, '.', 'DisplayName','prefit');  % mm/s
        plot(t_days, res_postfit.rr*1e6, '.', 'DisplayName','postfit');
        xlabel('Days since detection'); ylabel('RR resid [mm/s]');
        legend('Location','best'); title('Range-rate residuals');
    else
        plot(t_days, res_prefit.rr*1e6, '.', t_days, res_postfit.rr*1e6, '.'); grid on;
        xlabel('Days since detection'); ylabel('RR resid [mm/s]');
        legend('prefit','postfit'); title('Range-rate residuals');
    end

    % ---------- Print final estimate with 3σ ----------
    sig3 = 3*sqrt(max(diag(P_est),0));
    labels = {'x [km]','y [km]','z [km]','vx [km/s]','vy [km/s]','vz [km/s]',...
              'k_SRP','lat_4 [rad]','lon_4 [rad]','bias_rr [km/s]'};
    fprintf('\nFinal Batch LSQ estimate at ET0:\n');
    for i=1:10
        fprintf('%-12s  %.6e   (± %.6e 3σ)\n', labels{i}, X_est(i), sig3(i));
    end

    % ---------- Outputs ----------
    out.X_est = X_est;
    out.P_est = P_est;
    out.iterHist = iterHist;
    out.files = files;
    out.ET0 = ET0;
    out.res_prefit = res_prefit;
    out.res_postfit = res_postfit;

    % ===== nested helpers =====

    function files = resolveFiles(spec)
        spec = lower(strtrim(spec));
        if spec == "0to6" || spec == "0-6"
            files = "ASTE583_Project_LTB_Measurements_0-6D_Truth.csv";
        elseif spec == "6to14" || spec == "6-14"
            files = "ASTE583_Project_LTB_Measurements_6D-14D_Truth.csv";
        elseif spec == "0to14" || spec == "0-14" || spec == "all"
            files = ["ASTE583_Project_LTB_Measurements_0-6D_Truth.csv", ...
                     "ASTE583_Project_LTB_Measurements_6D-14D_Truth.csv"];
        else
            files = spec; % assume it's a filename
        end
    end

    function s = measSigma(stationID, type, sig, scale)
        if stationID == 4
            if type == "rho"
                s = sig.rho_ant * scale.rho_ant;
            else
                s = sig.rr_ant * scale.rr_ant;
            end
        else
            if type == "rho"
                s = sig.rho_dsn * scale.rho_dsn;
            else
                s = sig.rr_dsn * scale.rr_dsn;
            end
        end
    end

    function res = computeResiduals(X_epoch)
        % Propagate from ET0 to all unique times
        X_aug0 = [X_epoch(:); reshape(eye(10),100,1)];
        tspan  = [ET0; tU];
        [~, XaugR] = ode113(odeFun, tspan, X_aug0, odeOpts);
        XaugR = XaugR(2:end,:);

        XU_r   = XaugR(:,1:10);
        PhiU_r = XaugR(:,11:end); %#ok<NASGU> % not needed here

        res.rho = nan(N,1);
        res.rr  = nan(N,1);

        for kk=1:N
            tu_idx = idxU(kk);
            xk = XU_r(tu_idx,:).';

            biasActive = (tMeas(kk) <= params.bias_end_et);
            if ~biasActive
                xk(10) = 0;
            end

            [y_pred, ~] = measurementModel(tMeas(kk), xk, stationIDs(kk), params, st);

            if isfinite(obs_rho(kk))
                res.rho(kk) = obs_rho(kk) - y_pred(1);
            end
            res.rr(kk) = obs_rr(kk) - y_pred(2);
        end
    end
end
