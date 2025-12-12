function runBatchLSQ_0to6_driver()
    max_it = 11;
    tol    = 1e-6;
    [Xhat, info] = prelimBatchLSQ(max_it, tol);
    t_meas_hr = (info.tk - info.ET0)/3600;
    Niter     = size(info.res_pre_all,1);
    prefit_ms  = info.res_pre_all  * 1e3;
    postfit_ms = info.res_post_all * 1e3;
    maxPerFig = 4;
    baseTitle = 'Doppler (range-rate)';
    yUnit     = 'm/s';
    prefix    = 'batchLSQ_rr';
    makeResidualPlots(t_meas_hr, prefit_ms, postfit_ms, ...
                      Niter, maxPerFig, baseTitle, yUnit, prefix);
end