function out = runPrelimEKF(outDir, span, varargin)
    if nargin < 2 || isempty(span)
        span = '0to6';
    end
    span = lower(string(span));
    if nargin < 1 || isempty(outDir)
        outDir = fullfile(pwd, "figs_EKF_" + span);
    end
    if ~exist(outDir,'dir'); mkdir(outDir); end
    ip = inputParser;
    ip.CaseSensitive = false;
    ip.addParameter('fitThroughDay', []);
    ip.addParameter('ltmDay', 15);
    ip.addParameter('measDecim', 1);
    ip.addParameter('makePlots', true);
    ip.addParameter('savePlots', true);
    ip.parse(varargin{:});
    opt = ip.Results;
    [params, sc, st, X0, P0] = projectInit();
    [ET0, ~, ~] = initTime();
    params.t0_et = ET0;
    params.bias_end_et = ET0 + 6*params.tday;
    if isempty(opt.fitThroughDay)
        if span == "0to6"
            opt.fitThroughDay = 6;
        else
            opt.fitThroughDay = 14;
        end
    end
    fit_end_et = ET0 + opt.fitThroughDay * params.tday;
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
    [tk, idx] = sort(meas.t(:));
    stk   = meas.stationID(idx);
    z_rho = meas.range(idx);
    z_rr  = meas.rr(idx);
    measDecim = max(1, round(opt.measDecim));
    if measDecim > 1
        keep = 1:measDecim:numel(tk);
        tk    = tk(keep);
        stk   = stk(keep);
        z_rho = z_rho(keep);
        z_rr  = z_rr(keep);
    end
    N = numel(tk);
    xhat_k = X0;
    P_k    = P0;
    P0_scale = 1;
    P_k = (P0_scale^2) * P_k;
    sigma_rho_dsn = 1e-3;
    sigma_rho_ant = 1e-2;
    sigma_rr_dsn  = 1e-7;
    sigma_rr_ant  = 1e-6;
    Rscale_rr_dsn  = 5;
    Rscale_rr_ant  = 10;
    Rscale_rho_dsn = 5;
    Rscale_rho_ant = 5;
    q_acc_nom   = 1e-16;
    q_k_rw    = 0;
    q_bias_rw = 0;
    time_rel = zeros(N,1);
    x_hist   = zeros(10,N);
    P_hist   = zeros(10,10,N);
    res_pre_rho  = nan(N,1);
    res_post_rho = nan(N,1);
    res_pre_rr   = nan(N,1);
    res_post_rr  = nan(N,1);
    nis = nan(N,1);
    dof = nan(N,1);
    inFitMask = false(N,1);
    odeOpts = odeset('RelTol',1e-11,'AbsTol',1e-13);
    t_prev = ET0;
    for k = 1:N
        t_meas = tk(k);
        dt = t_meas - t_prev;
        if dt < 0, dt = 0; end
        inFit = (t_meas <= (fit_end_et + 1e-6));
        inFitMask(k) = inFit;
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
            Qd = makeQd(dt, q_acc_nom, q_k_rw, q_bias_rw);
            P_pred = Phi_k * P_k * Phi_k.' + Qd;
            P_pred = 0.5*(P_pred + P_pred.');
        end
        stationID = stk(k);
        useRange  = isfinite(z_rho(k));
        if stationID == 4
            sig_rho = Rscale_rho_ant * sigma_rho_ant;
            sig_rr  = Rscale_rr_ant  * sigma_rr_ant;
        else
            sig_rho = Rscale_rho_dsn * sigma_rho_dsn;
            sig_rr  = Rscale_rr_dsn  * sigma_rr_dsn;
        end
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
        y = z - h;
        if useRange
            res_pre_rho(k) = y(1);
            res_pre_rr(k)  = y(2);
        else
            res_pre_rr(k)  = y;
        end
        if inFit
            S = H * P_pred * H.' + R;
            K = (P_pred * H.') / S;
            if useRange
                nis(k) = y.' * (S \ y);
                dof(k) = 2;
            else
                nis(k) = (y^2) / S;
                dof(k) = 1;
            end
            p_gate = 0.999999;
            chi2_gate = chi2inv_local(p_gate, dof(k));
            R_eff = R;
            if nis(k) > chi2_gate
                scale = nis(k) / chi2_gate;
                R_eff = R * scale;
                S = H * P_pred * H.' + R_eff;
                K = (P_pred * H.') / S;
            end
            if ~biasActive
                K(10,:) = 0;
            end
            x_upd = x_pred + K * y;
            I10 = eye(10);
            P_upd = (I10 - K*H) * P_pred * (I10 - K*H).' + K * R_eff * K.';
            P_upd = 0.5*(P_upd + P_upd.');
        else
            x_upd = x_pred;
            P_upd = P_pred;
            nis(k) = NaN;
            dof(k) = NaN;
        end
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
        xhat_k = x_upd;
        P_k    = P_upd;
        t_prev = t_meas;
        x_hist(:,k)   = xhat_k;
        P_hist(:,:,k) = P_k;
        time_rel(k)   = t_meas - ET0;
        if mod(k, 5000) == 0
            fprintf('EKF %s: %d / %d (t=%.2f days)\n', span, k, N, time_rel(k)/86400);
        end
    end
    t_days = time_rel / 86400;
    fitMask = inFitMask & isfinite(t_days);
    passMask = ~fitMask;
    rr_pre_fit  = res_pre_rr(isfinite(res_pre_rr) & fitMask);
    rr_post_fit = res_post_rr(isfinite(res_post_rr) & fitMask);
    rr_rms_pre_fit  = sqrt(mean(rr_pre_fit.^2));
    rr_rms_post_fit = sqrt(mean(rr_post_fit.^2));
    fprintf('RR RMS (FIT) prefit  = %.6g km/s (%.3f mm/s)\n', rr_rms_pre_fit,  rr_rms_pre_fit*1e6);
    fprintf('RR RMS (FIT) postfit = %.6g km/s (%.3f mm/s)\n', rr_rms_post_fit, rr_rms_post_fit*1e6);
    rr_rms_pass = NaN;
    rr_pass = res_pre_rr(isfinite(res_pre_rr) & passMask);
    if ~isempty(rr_pass)
        rr_rms_pass = sqrt(mean(rr_pass.^2));
        fprintf('RR RMS (PASSTHRU)    = %.6g km/s (%.3f mm/s)\n', rr_rms_pass, rr_rms_pass*1e6);
    end
    rho_rms_pre_fit  = NaN;
    rho_rms_post_fit = NaN;
    rho_rms_pass     = NaN;
    if any(isfinite(res_pre_rho))
        rho_pre_fit  = res_pre_rho(isfinite(res_pre_rho) & fitMask);
        rho_post_fit = res_post_rho(isfinite(res_post_rho) & fitMask);
        if ~isempty(rho_pre_fit),  rho_rms_pre_fit  = sqrt(mean(rho_pre_fit.^2));  end
        if ~isempty(rho_post_fit), rho_rms_post_fit = sqrt(mean(rho_post_fit.^2)); end
        fprintf('Rho RMS (FIT) prefit  = %.6g km (%.3f m)\n', rho_rms_pre_fit,  rho_rms_pre_fit*1e3);
        fprintf('Rho RMS (FIT) postfit = %.6g km (%.3f m)\n', rho_rms_post_fit, rho_rms_post_fit*1e3);
        rho_pass = res_pre_rho(isfinite(res_pre_rho) & passMask);
        if ~isempty(rho_pass)
            rho_rms_pass = sqrt(mean(rho_pass.^2));
            fprintf('Rho RMS (PASSTHRU)    = %.6g km (%.3f m)\n', rho_rms_pass, rho_rms_pass*1e3);
        end
    end
    t_ltm_et = ET0 + opt.ltmDay * params.tday;
    x_ltm = xhat_k;
    P_ltm = P_k;
    if t_prev < t_ltm_et
        dt = t_ltm_et - t_prev;
        X_aug0 = [xhat_k; reshape(eye(10),100,1)];
        odeFun = @(t,X) projectDynamics(t, X, params, sc);
        [~, Xtraj] = ode113(odeFun, [t_prev t_ltm_et], X_aug0, odeOpts);
        Xf_aug = Xtraj(end,:).';
        x_ltm = Xf_aug(1:10);
        Phi   = reshape(Xf_aug(11:end),10,10);
        Qd    = makeQd(dt, q_acc_nom, q_k_rw, q_bias_rw);
        P_ltm = Phi*P_k*Phi.' + Qd;
        P_ltm = 0.5*(P_ltm + P_ltm.');
    end
    [sigmaT_ltm_km, sigmaR_ltm_km, sigmaN_ltm_km] = sigmaRTN(x_ltm, P_ltm);
    fprintf('LTM uncertainty (1-sigma, RTN): sigmaR=%.2f km, sigmaT=%.2f km, sigmaN=%.2f km\n', ...
        sigmaR_ltm_km, sigmaT_ltm_km, sigmaN_ltm_km);
    if sigmaT_ltm_km < 100
        fprintf('✓ Meets along-track requirement at LTM: sigmaT < 100 km\n');
    else
        fprintf('✗ DOES NOT meet along-track requirement at LTM: sigmaT >= 100 km\n');
    end
    if opt.makePlots
        try
            t_days_plot = [t_days; opt.ltmDay];
            stk_plot    = [stk; NaN];
            res_pre_rr_plot   = [res_pre_rr; NaN];
            res_post_rr_plot  = [res_post_rr; NaN];
            res_pre_rho_plot  = [res_pre_rho; NaN];
            res_post_rho_plot = [res_post_rho; NaN];
            nis_plot = [nis; NaN];
            dof_plot = [dof; NaN];
            x_hist_plot = [x_hist, x_ltm];
            P_hist_plot = cat(3, P_hist, P_ltm);
            fitMask_plot = (t_days_plot <= (opt.fitThroughDay + 1e-12));
            stateLabels = {'x','y','z','v_x','v_y','v_z', ...
               'k_{SRP}','lat_4','lon_4','rho_bias_dot'};
            makeEkfPlots(span, t_days_plot, stk_plot, ...
                         res_pre_rr_plot, res_post_rr_plot, ...
                         res_pre_rho_plot, res_post_rho_plot, ...
                         x_hist_plot, P_hist_plot, ...
                         nis_plot, dof_plot, ...
                         'fitMask', fitMask_plot, ...
                         'dcoDay', opt.fitThroughDay, ...
                         'ltmDay', opt.ltmDay, ...
                         'stateLabels', stateLabels);
        catch ME
            warning('runPrelimEKF:PlottingFailed', 'EKF plotting failed: %s', ME.message);
        end
    end
    out.span           = span;
    out.outDir         = outDir;
    out.fitThroughDay  = opt.fitThroughDay;
    out.ltmDay         = opt.ltmDay;
    out.X0             = X0;
    out.P0             = P0;
    out.ET0            = ET0;
    out.tk             = tk;
    out.station        = stk;
    out.x_final        = xhat_k;
    out.P_final        = P_k;
    out.x_hist         = x_hist;
    out.P_hist         = P_hist;
    out.t_days         = t_days;
    out.inFitMask      = inFitMask;
    out.res_pre_rho    = res_pre_rho;
    out.res_post_rho   = res_post_rho;
    out.res_pre_rr     = res_pre_rr;
    out.res_post_rr    = res_post_rr;
    out.nis            = nis;
    out.dof            = dof;
    out.rr_rms_pre_fit   = rr_rms_pre_fit;
    out.rr_rms_post_fit  = rr_rms_post_fit;
    out.rr_rms_pass      = rr_rms_pass;
    out.rho_rms_pre_fit  = rho_rms_pre_fit;
    out.rho_rms_post_fit = rho_rms_post_fit;
    out.rho_rms_pass     = rho_rms_pass;
    out.x_ltm          = x_ltm;
    out.P_ltm          = P_ltm;
    out.sigmaT_ltm_km  = sigmaT_ltm_km;
    out.sigmaR_ltm_km  = sigmaR_ltm_km;
    out.sigmaN_ltm_km  = sigmaN_ltm_km;
end
function Qd = makeQd(dt, q_acc, q_k_rw, q_bias_rw)
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
    x = 2 * gammaincinv(p, k/2);
end
function [sigmaT, sigmaR, sigmaN] = sigmaRTN(x, P)
    r = x(1:3);
    v = x(4:6);
    Ppos = P(1:3,1:3);
    Rhat = r / max(norm(r), eps);
    h = cross(r, v);
    if norm(h) < 1e-12
        tmp = [1;0;0];
        if abs(dot(tmp,Rhat)) > 0.9; tmp = [0;1;0]; end
        Nhat = cross(Rhat, tmp);
        Nhat = Nhat / max(norm(Nhat), eps);
    else
        Nhat = h / norm(h);
    end
    That = cross(Nhat, Rhat);
    That = That / max(norm(That), eps);
    A = [Rhat.'; That.'; Nhat.']; % inertial -> RTN
    P_rtn = A * Ppos * A.';
    s = sqrt(max(diag(P_rtn), 0));
    sigmaR = s(1);
    sigmaT = s(2);
    sigmaN = s(3);
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