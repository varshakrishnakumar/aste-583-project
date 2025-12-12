function [sigma_RTN, P_r_RTN] = computeAlongTrackSigmaLTM()
    [params, sc, st, X0, P0] = projectInit();
    [ET0, LTM, LCM]          = initTime();
    opts            = struct();
    opts.maxIter    = 11;
    opts.earlyMode  = "drop";
    opts.earlyDays  = 2;
    opts.earlyScale = 10;
    opts.makePlots  = false;
    out = runBatchLSQ("0to14","ekf_k",opts);
    deltaX0 = out.X_est - X0;
    P_final = out.P_est;
    X0fit = X0 + deltaX0;
    X0fit(7) = 0.6996824;
    Pfit = P_final;
    odeOpts = odeset('RelTol',1e-11,'AbsTol',1e-13);
    odeFun  = @(t,X) projectDynamics(t, X, params, sc);
    X_aug0 = [X0fit; reshape(eye(10),100,1)];
    [~, X_aug] = ode45(odeFun, [ET0 LTM], X_aug0, odeOpts);
    X_LTM_aug  = X_aug(end,:).';
    Phi_LTM = reshape(X_LTM_aug(11:end),10,10);
    P_LTM = Phi_LTM * Pfit * Phi_LTM.';
    P_LTM = 0.5 * (P_LTM + P_LTM.');   % enforce symmetry
    P_r_inertial = P_LTM(1:3,1:3);
    xEarth = cspice_spkezr('EARTH', LTM, 'ECLIPJ2000', 'NONE', 'SUN');
    r_e    = xEarth(1:3).';
    v_e    = xEarth(4:6).';
    r_sc = X_LTM_aug(1:3);
    v_sc = X_LTM_aug(4:6);
    r_rel = r_sc - r_e;
    v_rel = v_sc - v_e;
    rhat = r_rel / norm(r_rel);
    h    = cross(r_rel, v_rel);
    zhat = h / norm(h);
    that = cross(zhat, rhat);
    R_RTN_INERTIAL = [rhat.'; that.'; zhat.'];
    P_r_RTN = R_RTN_INERTIAL * P_r_inertial * R_RTN_INERTIAL.';
    P_r_RTN = 0.5 * (P_r_RTN + P_r_RTN.');   % symmetry cleanup
    sigma_RTN = sqrt(max(diag(P_r_RTN),0));
    sigma_R   = sigma_RTN(1);
    sigma_T   = sigma_RTN(2);
    sigma_N   = sigma_RTN(3);
    fprintf('\n---- LTM Position Uncertainty (from 0–14d batch cov) ----\n');
    fprintf('1-sigma radial      = %8.3f km\n', sigma_R);
    fprintf('1-sigma along-track = %8.3f km\n', sigma_T);
    fprintf('1-sigma cross-track = %8.3f km\n', sigma_N);
    fprintf('Transverse requirement (1-sigma) = %6.1f km\n', params.sigmaTransverse);
    if sigma_T < params.sigmaTransverse
        fprintf('=> Requirement MET in linear covariance sense.\n');
    else
        fprintf('=> Requirement NOT met in linear covariance sense.\n');
    end
end