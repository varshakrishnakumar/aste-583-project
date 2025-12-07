function out = runPrelimEKF(outDir, span)
% runPrelimEKF  EKF orbit determination for ASTE583 project.
%
% Usage:
%   out = runPrelimEKF();                       % default span = '0to6'
%   out = runPrelimEKF([], '0to14');            % run both files together
%   out = runPrelimEKF('figs_EKF_6to14','6to14');
%
% span options:
%   '0to6'  : only the 0–6 day file (RR only, bias active)
%   '6to14' : only the 6–14 day file (range+RR, bias inactive)
%   '0to14' : both files concatenated (recommended final run)

    if nargin < 2 || isempty(span)
        span = '0to6'; % keep legacy behavior
    end
    span = lower(string(span));

    if nargin < 1 || isempty(outDir)
        outDir = fullfile(pwd, "figs_EKF_" + span);
    end
    if ~exist(outDir,'dir'); mkdir(outDir); end

    % ---------- Init ----------
    [params, sc, st, X0, P0] = projectInit();
    [ET0, ~, ~] = initTime();

    % IMPORTANT for station rotation and GST0 meaning:
    params.t0_et = ET0;

    % Project rule: RR bias exists only during first 6 days
    params.bias_end_et = ET0 + 6*params.tday;

    % ---------- Load measurements ----------
    file_0to6  = 'ASTE583_Project_LTB_Measurements_0-6D_Truth.csv';
    file_6to14 = 'ASTE583_Project_LTB_Measurements_6D-14D_Truth.csv';

    switch span
        case "0to6"
            meas = parseMeasurementData(file_0to6, ET0);

        case "6to14"
            meas = parseMeasurementData(file_6to14, ET0);

        case "0to14"
            m1 = parseMeasurementData(file_0to6,  ET0);
            m2 = parseMeasurementData(file_6to14, ET0);
            meas.t         = [m1.t;         m2.t];
            meas.stationID = [m1.stationID; m2.stationID];
            meas.range     = [m1.range;     m2.range];
            meas.rr        = [m1.rr;        m2.rr];

        otherwise
            error('runPrelimEKF:BadSpan', 'span must be ''0to6'', ''6to14'', or ''0to14''.');
    end

    % Sort by time
    [tk, idx] = sort(meas.t(:));
    stk   = meas.stationID(idx);
    z_rho = meas.range(idx); % [km] (NaN in 0–6)
    z_rr  = meas.rr(idx);    % [km/s]

    measDecim = 1;
    if measDecim > 1
        keep = 1:measDecim:numel(tk);
        tk = tk(keep); stk = stk(keep); z_rho = z_rho(keep); z_rr = z_rr(keep);
    end

    N = numel(tk);

    % ---------- Filter initial condition ----------
    xhat_k = X0;
    P_k    = P0;

    % Optional covariance inflation
    P0_scale = 1;
    P_k = (P0_scale^2) * P_k;

    % ---------- Measurement noise model ----------
    % Range:
    %   DSN: 1 m = 0.001 km
    %   ANT: 10 m = 0.01 km
    sigma_rho_dsn = 1e-3;   % km
    sigma_rho_ant = 1e-2;   % km

    % Range-rate:
    %   DSN: 0.1 mm/s = 1e-7 km/s
    %   ANT: 1.0 mm/s = 1e-6 km/s
    sigma_rr_dsn  = 1e-7;   % km/s
    sigma_rr_ant  = 1e-6;   % km/s

    % Optional scaling knobs (if residuals are “too big”)
    Rscale_rr_dsn  = 5;
    Rscale_rr_ant  = 10;
    Rscale_rho_dsn = 5;
    Rscale_rho_ant = 5;


    % ---------- Process noise ----------
    % White acceleration noise (km^2/s^3). Tune as needed.
    q_acc_nom   = 1e-16;
    % q_acc_outg  = 1e-15;    % first 2 days
    % outgas_window = 2*86400;

    % Random walk on parameters (often 0 for these)
    q_k_rw    = 0;
    q_bias_rw = 0;

    % ---------- Storage ----------
    time_rel = zeros(N,1);
    x_hist   = zeros(10,N);
    sig_hist = zeros(10,N);
    P_hist   = zeros(10,10,N);

    % Residual storage (some are NaN when range not present)
    res_pre_rho  = nan(N,1);
    res_post_rho = nan(N,1);
    res_pre_rr   = nan(N,1);
    res_post_rr  = nan(N,1);

    nis = nan(N,1);
    dof = nan(N,1);

    % ---------- Prop options ----------
    odeOpts = odeset('RelTol',1e-11,'AbsTol',1e-13);

    t_prev = ET0;

    % ---------- EKF loop ----------
    for k = 1:N
        t_meas = tk(k);
        dt = t_meas - t_prev;
        if dt < 0, dt = 0; end

        % ---- Propagate state + STM ----
        if dt < 1e-9
            x_pred = xhat_k;
            Phi_k  = eye(10);
            P_pred = P_k;
        else
            X_aug0 = [xhat_k; reshape(eye(10),100,1)];
            odeFun = @(t,X) projectDynamics(t, X, params, sc);
            [~, X_traj] = ode113(odeFun, [t_prev t_meas], X_aug0, odeOpts);
            Xk_aug = X_traj(end,:).';

            x_pred = Xk_aug(1:10);
            Phi_k  = reshape(Xk_aug(11:end), 10, 10);

            % Process noise for this step
            % t_since_det = t_prev - ET0;
            % if t_since_det < outgas_window
            %     q_acc = q_acc_outg;
            % else
            q_acc = q_acc_nom;
            % end
            Qd = makeQd(dt, q_acc, q_k_rw, q_bias_rw);

            P_pred = Phi_k * P_k * Phi_k.' + Qd;
            P_pred = 0.5*(P_pred + P_pred.');
        end

        % ---- Measurement selection: 1D (rr only) or 2D (range+rr) ----
        stationID = stk(k);
        useRange  = isfinite(z_rho(k));   % only true in 6–14 file

        % Station-dependent sigmas
        if stationID == 4
            sig_rho = Rscale_rho_ant * sigma_rho_ant;
            sig_rr  = Rscale_rr_ant  * sigma_rr_ant;
        else
            sig_rho = Rscale_rho_dsn * sigma_rho_dsn;
            sig_rr  = Rscale_rr_dsn  * sigma_rr_dsn;
        end

        % Bias only active in first 6 days. To be robust even if
        % measurementModel always adds bias, force bias=0 after cutoff and
        % zero its sensitivity.
        biasActive = (t_meas <= params.bias_end_et);
        x_for_meas = x_pred;
        if ~biasActive
            x_for_meas(10) = 0;
        end

        [z_pred2, H2] = measurementModel(t_meas, x_for_meas, stationID, params, st);
        if ~biasActive
            H2(2,10) = 0;
        end

        if useRange
            z = [z_rho(k); z_rr(k)];
            h = [z_pred2(1); z_pred2(2)];
            H = [H2(1,:);    H2(2,:)];
            R = diag([sig_rho^2, sig_rr^2]);
        else
            z = z_rr(k);
            h = z_pred2(2);
            H = H2(2,:);
            R = sig_rr^2;
        end

        y = z - h;  % innovation / prefit residual(s)

        if useRange
            res_pre_rho(k) = y(1);
            res_pre_rr(k)  = y(2);
        else
            res_pre_rr(k)  = y;
        end

        % ---- EKF update ----
        S = H * P_pred * H.' + R;
        K = (P_pred * H.') / S;

        % NIS
        if useRange
            nis(k) = y.' * (S \ y);
            dof(k) = 2;
        else
            nis(k) = (y^2) / S;
            dof(k) = 1;
        end

        % Soft gating: inflate R if NIS is huge
        p_gate = 0.999999; % very conservative
        chi2_gate = chi2inv_local(p_gate, dof(k));
        R_eff = R;

        if nis(k) > chi2_gate
            scale = nis(k) / chi2_gate;
            R_eff = R * scale;
            S = H * P_pred * H.' + R_eff;
            K = (P_pred * H.') / S;
        end

        % Freeze bias state from being “updated” by unbiased data
        if ~biasActive
            K(10,:) = 0;
        end

        x_upd = x_pred + K * y;

        % Joseph-form covariance update
        I10 = eye(10);
        P_upd = (I10 - K*H) * P_pred * (I10 - K*H).' + K * R_eff * K.';
        P_upd = 0.5*(P_upd + P_upd.');

        % ---- Postfit residuals (re-evaluate model) ----
        x_post_meas = x_upd;
        if ~biasActive
            x_post_meas(10) = 0;
        end
        z_post2 = measurementModel(t_meas, x_post_meas, stationID, params, st);

        if useRange
            res_post_rho(k) = z_rho(k) - z_post2(1);
            res_post_rr(k)  = z_rr(k)  - z_post2(2);
        else
            res_post_rr(k)  = z_rr(k)  - z_post2(2);
        end

        % ---- Log ----
        xhat_k = x_upd;
        P_k    = P_upd;
        t_prev = t_meas;

        x_hist(:,k)   = xhat_k;
        sig_hist(:,k) = sqrt(max(diag(P_k),0));
        P_hist(:,:,k) = P_k;
        time_rel(k)   = t_meas - ET0;

        if mod(k, 5000) == 0
            fprintf('EKF %s: %d / %d (t=%.2f days)\n', span, k, N, time_rel(k)/86400);
        end
    end

    % ---------- Stats ----------
    t_days = time_rel / 86400;

    rr_pre  = res_pre_rr(isfinite(res_pre_rr));
    rr_post = res_post_rr(isfinite(res_post_rr));

    fprintf('RR RMS prefit  = %.6g km/s (%.3f mm/s)\n', sqrt(mean(rr_pre.^2)),  sqrt(mean(rr_pre.^2))*1e6);
    fprintf('RR RMS postfit = %.6g km/s (%.3f mm/s)\n', sqrt(mean(rr_post.^2)), sqrt(mean(rr_post.^2))*1e6);

    if any(isfinite(res_pre_rho))
        rho_pre  = res_pre_rho(isfinite(res_pre_rho));
        rho_post = res_post_rho(isfinite(res_post_rho));
        fprintf('Rho RMS prefit  = %.6g km (%.3f m)\n', sqrt(mean(rho_pre.^2)),  sqrt(mean(rho_pre.^2))*1e3);
        fprintf('Rho RMS postfit = %.6g km (%.3f m)\n', sqrt(mean(rho_post.^2)), sqrt(mean(rho_post.^2))*1e3);
    end

    % ---------- Plots ----------
    uSt = unique(stk(:)).';
    C = lines(numel(uSt));

    % (A) Range-rate residuals by station
    f1 = figure('Color','w','Name','RR residuals by station');
    subplot(2,1,1); hold on; grid on;
    for i = 1:numel(uSt)
        ii = (stk == uSt(i)) & isfinite(res_pre_rr);
        scatter(t_days(ii), res_pre_rr(ii)*1e3, 8, C(i,:), 'filled');
    end
    yline(0,'k-');
    xlabel('Days since detection'); ylabel('Prefit \rhȯ resid [m/s]');
    legend(stationLegend(uSt), 'Location','bestoutside');
    title('Prefit range-rate residuals');

    subplot(2,1,2); hold on; grid on;
    for i = 1:numel(uSt)
        ii = (stk == uSt(i)) & isfinite(res_post_rr);
        scatter(t_days(ii), res_post_rr(ii)*1e3, 8, C(i,:), 'filled');
    end
    yline(0,'k-');
    xlabel('Days since detection'); ylabel('Postfit \rhȯ resid [m/s]');
    legend(stationLegend(uSt), 'Location','bestoutside');
    title('Postfit range-rate residuals');

    % (B) Range residuals (only where present)
    if any(isfinite(res_pre_rho))
        fR = figure('Color','w','Name','Range residuals by station');
        subplot(2,1,1); hold on; grid on;
        for i = 1:numel(uSt)
            ii = (stk == uSt(i)) & isfinite(res_pre_rho);
            scatter(t_days(ii), res_pre_rho(ii)*1e3, 8, C(i,:), 'filled'); % km->m
        end
        yline(0,'k-');
        xlabel('Days since detection'); ylabel('Prefit \rho resid [m]');
        legend(stationLegend(uSt), 'Location','bestoutside');
        title('Prefit range residuals');

        subplot(2,1,2); hold on; grid on;
        for i = 1:numel(uSt)
            ii = (stk == uSt(i)) & isfinite(res_post_rho);
            scatter(t_days(ii), res_post_rho(ii)*1e3, 8, C(i,:), 'filled');
        end
        yline(0,'k-');
        xlabel('Days since detection'); ylabel('Postfit \rho resid [m]');
        legend(stationLegend(uSt), 'Location','bestoutside');
        title('Postfit range residuals');
    end

    % (C) NIS
    f2 = figure('Color','w','Name','NIS');
    semilogy(t_days, nis, '.'); grid on; hold on;
    xlabel('Days since detection'); ylabel('NIS');
    title('Normalized Innovation Squared');

    % (D) k_SRP +/- 3σ
    f3 = figure('Color','w','Name','kSRP');
    k_est = x_hist(7,:); k_sig = sig_hist(7,:);
    plot(t_days, k_est, 'LineWidth',1.5); hold on; grid on;
    plot(t_days, k_est+3*k_sig,'--');
    plot(t_days, k_est-3*k_sig,'--');
    xlabel('Days since detection'); ylabel('k_{SRP}');
    legend('estimate','\pm 3\sigma','Location','best');
    title('SRP scale factor estimate');

    % (E) Station #4 lat/lon +/- 3σ
    f4 = figure('Color','w','Name','Antarctica latlon');
    lat_est = x_hist(8,:); lat_sig = sig_hist(8,:);
    lon_est = x_hist(9,:); lon_sig = sig_hist(9,:);
    subplot(2,1,1); hold on; grid on;
    plot(t_days, lat_est*180/pi,'LineWidth',1.5);
    plot(t_days,(lat_est+3*lat_sig)*180/pi,'--');
    plot(t_days,(lat_est-3*lat_sig)*180/pi,'--');
    xlabel('Days'); ylabel('Lat_4 [deg]');
    legend('estimate','\pm 3\sigma','Location','best');
    title('Station #4 latitude');

    subplot(2,1,2); hold on; grid on;
    plot(t_days, lon_est*180/pi,'LineWidth',1.5);
    plot(t_days,(lon_est+3*lon_sig)*180/pi,'--');
    plot(t_days,(lon_est-3*lon_sig)*180/pi,'--');
    xlabel('Days'); ylabel('Lon_4 [deg]');
    legend('estimate','\pm 3\sigma','Location','best');
    title('Station #4 longitude');

    % (F) RR bias +/- 3σ (will be flat after day 6 by design)
    f5 = figure('Color','w','Name','RR bias');
    b_est = x_hist(10,:); b_sig = sig_hist(10,:);
    plot(t_days, b_est*1e3,'LineWidth',1.5); hold on; grid on;
    plot(t_days, (b_est+3*b_sig)*1e3,'--');
    plot(t_days, (b_est-3*b_sig)*1e3,'--');
    xline(6,'k--','Bias off');
    xlabel('Days since detection'); ylabel('\rhȯ_{bias} [m/s]');
    legend('estimate','\pm 3\sigma','Location','best');
    title('Range-rate bias estimate');

    % (G) 2D delta-r XY + 3σ ellipse (final)
    f7 = figure('Color','w','Name','2D delta-r XY + ellipse');
    dr = x_hist(1:3,:) - X0(1:3);
    plot(dr(1,:), dr(2,:), 'LineWidth',1.2); hold on; grid on;
    xlabel('\Delta r_x [km]'); ylabel('\Delta r_y [km]');
    title('2D view: \Delta r_x vs \Delta r_y with final 3\sigma ellipse');

    Pxy = P_hist(1:2,1:2,end);
    mu  = [dr(1,end); dr(2,end)];
    plotCovEllipse(mu, Pxy, 3);

    % ---------- Save figures ----------
    % saveAllFigures(outDir, "EKF_" + span);
    % 
    % try
    %     print(f, fullfile(outDir, base), '-dpng', '-r300');
    %     print(f, fullfile(outDir, base), '-dpdf', '-r300');
    % catch ME
    %     warning('Failed to save figure %d: %s', f.Number, ME.message);
    % end


    % ---------- Return useful outputs ----------
    out.span     = span;
    out.outDir   = outDir;
    out.ET0      = ET0;
    out.tk       = tk;
    out.station  = stk;
    out.x_final  = xhat_k;
    out.P_final  = P_k;
    out.x_hist   = x_hist;
    out.P_hist   = P_hist;
    out.res_pre_rho  = res_pre_rho;
    out.res_post_rho = res_post_rho;
    out.res_pre_rr   = res_pre_rr;
    out.res_post_rr  = res_post_rr;
    out.nis      = nis;
end

% -------------------- helpers -------------------- %

function Qd = makeQd(dt, q_acc, q_k_rw, q_bias_rw)
% Discrete process noise for [r;v] with white accel + random-walk params
    Qd = zeros(10,10);
    if dt <= 0; return; end

    I3 = eye(3);
    Qrr = (dt^3/3) * q_acc * I3;
    Qrv = (dt^2/2) * q_acc * I3;
    Qvv = (dt)     * q_acc * I3;

    Qd(1:3,1:3) = Qrr;
    Qd(1:3,4:6) = Qrv;
    Qd(4:6,1:3) = Qrv;
    Qd(4:6,4:6) = Qvv;

    Qd(7,7)   = q_k_rw    * dt;
    Qd(10,10) = q_bias_rw * dt;
end

function x = chi2inv_local(p, k)
%CHI2INV_LOCAL Chi-square inverse CDF without requiring Statistics Toolbox.
% If X ~ chi2(k), then X/2 ~ Gamma(k/2, 1), so:
%   chi2inv(p,k) = 2 * gammaincinv(p, k/2)
    x = 2 * gammaincinv(p, k/2);
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
% Plot nsig covariance ellipse for 2x2 covariance
    P2 = 0.5*(P2+P2.');
    [V,D] = eig(P2);
    D = max(D, 0);
    A = V*sqrt(D);

    th = linspace(0, 2*pi, 200);
    circ = [cos(th); sin(th)];
    ell = mu + nsig * (A * circ);

    plot(ell(1,:), ell(2,:), 'k--', 'LineWidth',1.2);
    scatter(mu(1), mu(2), 25, 'k', 'filled');
end

function saveAllFigures(outDir, prefix)
% Saves all open figures as .png + .pdf + .fig into outDir
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
