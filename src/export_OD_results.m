function export_OD_results()
    opts = struct();
    opts.maxIter   = 20;
    opts.tol_dx    = 1e-6;
    opts.useLM     = false;
    opts.robust    = false;
    opts.huberC    = 3.0;
    opts.earlyMode = "drop";
    opts.earlyDays = 2;
    opts.earlyScale= 10;
    ekf06 = runPrelimEKF([], '0to6');
    ekf014 = runPrelimEKF([], '0to14');
    batch06 = runBatchLSQ('0to6','ekf_k',opts);
    batch014 = runBatchLSQ('0to14','ekf_k',opts);
    save('OD_all_cases.mat','ekf06','ekf014','batch06','batch014');
    fprintf('Saved all four cases to OD_all_cases.mat\n');
    summary = buildSummaryTable(ekf06,  "EKF",   "0to6",  true);
    summary = [summary;
               buildSummaryTable(ekf014, "EKF",   "0to14", true);
               buildSummaryTable(batch06,"BATCH", "0to6",  false);
               buildSummaryTable(batch014,"BATCH","0to14", false)];
    writetable(summary, 'OD_case_summary.csv');
    fprintf('Wrote OD_case_summary.csv\n');
    exportResidualsEKF(ekf06,   'EKF_0to6_residuals.csv');
    exportResidualsEKF(ekf014,  'EKF_0to14_residuals.csv');
    exportResidualsBatch(batch06, 'BATCH_0to6_residuals.csv');
    exportResidualsBatch(batch014,'BATCH_0to14_residuals.csv');
    fprintf('Wrote residual CSVs for all cases.\n');
end
function T = buildSummaryTable(out, filterName, spanLabel, isEKF)
    if isEKF
        x_prior = out.X0(:);
        x_final = out.x_final(:);
        P_final = out.P_final;
    else
        x_prior = out.X0(:);
        x_final = out.X_est(:);
        P_final = out.P_est;
    end
    labels = {'x','y','z','vx','vy','vz','k_SRP','lat4','lon4','bias_rr'};
    T = table();
    T.filter = string(filterName);
    T.span   = string(spanLabel);
    for i = 1:10
        priorName = sprintf('%s_prior', labels{i});
        estName   = sprintf('%s_est',   labels{i});
        sigName   = sprintf('%s_sig',   labels{i});
        T.(priorName) = x_prior(i);
        T.(estName)   = x_final(i);
        T.(sigName)   = sqrt(max(P_final(i,i),0));
    end
end
function exportResidualsEKF(out, filename)
    t_days = (out.tk(:) - out.ET0) / 86400;
    st     = out.station(:);
    pre_rr  = out.res_pre_rr(:);
    post_rr = out.res_post_rr(:);
    pre_rho = out.res_pre_rho(:);
    post_rho= out.res_post_rho(:);
    T = table();
    T.t_days   = t_days;
    T.station  = st;
    T.prefit_rr_kmps   = pre_rr;
    T.postfit_rr_kmps  = post_rr;
    T.prefit_rho_km    = pre_rho;
    T.postfit_rho_km   = post_rho;
    writetable(T, filename);
end
function exportResidualsBatch(out, filename)
    t_days = (out.tMeas(:) - out.ET0) / 86400;
    st     = out.stationIDs(:);
    pre_rr  = out.res_prefit.rr(:);
    post_rr = out.res_postfit.rr(:);
    pre_rho = out.res_prefit.rho(:);
    post_rho= out.res_postfit.rho(:);
    T = table();
    T.t_days   = t_days;
    T.station  = st;
    T.prefit_rr_kmps   = pre_rr;
    T.postfit_rr_kmps  = post_rr;
    T.prefit_rho_km    = pre_rho;
    T.postfit_rho_km   = post_rho;
    writetable(T, filename);
end