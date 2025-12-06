function debug_prof_case()
    [ET0,~,~] = initTime();                

    [params, sc, st, X0, ~] = projectInit();
    params.t0_et = ET0;                      

    dt_days = 6.060763888888889;
    t_test  = ET0 + dt_days*params.tday;

    opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
    [~, Xout] = ode113(@(t,x) projectDynamics(t,x,params,sc), [ET0 t_test], X0, opts);
    Xtest = Xout(end,1:10).';

    [y, ~] = measurementModel(t_test, Xtest, 1, params, st);

    y_ref = [5.894311564735891e7; 1.885427030618933];

    fprintf('our range      : %.15e km\n', y(1));
    fprintf('Ref  range     : %.15e km\n', y_ref(1));
    fprintf('Diff           : %.6f km\n\n', y(1) - y_ref(1));

    fprintf('our range-rate : %.15e km/s\n', y(2));
    fprintf('Ref  range-rate: %.15e km/s\n', y_ref(2));
    fprintf('Diff           : %.9e km/s (%.6f mm/s)\n', y(2)-y_ref(2), (y(2)-y_ref(2))*1e6);
end
