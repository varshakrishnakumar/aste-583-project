function runPrelimEKF()
% runPrelimEKF  Preliminary OD using EKF on 0-6 day data.
%
% Produces:
%   - Prefit & postfit range-rate residuals vs time
%   - k_SRP estimate with 3-sigma bounds
%   - Antarctica station lat/lon estimates with 3-sigma
%   - Range-rate bias estimate with 3-sigma

    % ---------- Init t_project & time ----------
    [params, sc, st, X0, P0] = projectInit();

    [ET0, ~, ~] = initTime();   % detection epoch ET in sec past J2000

    % ---------- Load measurements (0-6 days) ----------
    meas = parseMeasurementData('ASTE583_Project_LTB_Measurements_0-6D_Truth.csv', ET0);

    N   = numel(meas.t);
    tk  = meas.t;           % ET times
    stk = meas.stationID;   % station IDs
    z_rr = meas.rr;         % truth range-rate [km/s]

    % ---------- Filter settings ----------
    xhat = X0;              % 10x1
    P    = P0;              % 10x10

    % Process noise Q (start with small or zero; tune later)
    Q = zeros(10,10);
    % Example: allow a small random walk on k_SRP and bias:
    q_k    = 1e-12;  % (unitless^2 / s)
    q_bias = 1e-12;  % ((km/s)^2 / s)
    Q(7,7)   = q_k;
    Q(10,10) = q_bias;

    % Measurement noise variance for range-rate 
    sigma_rr = 1e-3;          % 1 m/s
    R_rr     = sigma_rr^2;    % scalar

    % ---------- Storage for analysis ----------
    time_rel = zeros(N,1);    % seconds since detection
    x_hist   = zeros(10,N);
    sig_hist = zeros(10,N);   % 1-sigma from P
    res_pre  = zeros(N,1);    % prefit residual (range-rate)
    res_post = zeros(N,1);    % postfit residual (range-rate)

    t_prev = tk(1);
    xhat_k = xhat;
    P_k    = P;

    % ---------- Main EKF loop ----------
    for k = 1:N

        t_meas   = tk(k);
        dt       = t_meas - t_prev;
    
        % ----- Propagate state + STM from t_prev to t_meas -----
        if abs(dt) < 1e-9
            x_pred = xhat_k;
            Phi_k  = eye(10);
            P_pred = P_k;   % no process noise yet
        else
            X_aug0 = [xhat_k; reshape(eye(10),100,1)];
    
            odeFun = @(t,X) projectDynamics(t,X,params,sc);
            [~, X_traj] = ode113(odeFun, [t_prev t_meas], X_aug0);
            Xk_aug = X_traj(end,:).';
    
            x_pred = Xk_aug(1:10);
            Phi_k  = reshape(Xk_aug(11:end), 10, 10);
    
            % Covariance propagation
            P_pred = Phi_k * P_k * Phi_k.' + Q * max(dt,0);
        end

        % ----- Measurement prediction -----
        stationID = stk(k);

        [z_pred, Hk] = measurementModel(t_meas, x_pred, stationID, params, st);

        % We only have range-rate for 0-6 days => use 2nd component
        z_rr_pred = z_pred(2);
        H_rr      = Hk(2,:);      % 1x10

        % Prefit residual
        y_rr = z_rr(k) - z_rr_pred;
        res_pre(k) = y_rr;

        % Kalman gain: scalar measurement
        S   = H_rr * P_pred * H_rr.' + R_rr;
        K   = (P_pred * H_rr.') / S;   % 10x1

        % State update
        x_upd = x_pred + K * y_rr;

        % Covariance update (Joseph form optional; simple form here)
        I10 = eye(10);
        P_upd = (I10 - K*H_rr) * P_pred;

        % Postfit residual (optional, using updated state)
        [z_pred_post, ~] = measurementModel(t_meas, x_upd, stationID, params, st);
        res_post(k) = z_rr(k) - z_pred_post(2);

        % Log
        xhat_k           = x_upd;
        P_k              = P_upd;
        x_hist(:,k)      = xhat_k;
        sig_hist(:,k)    = sqrt(diag(P_k));
        time_rel(k)      = t_meas - ET0;

        t_prev = t_meas;
    end

    % ---------- Plots / analysis ----------

    t_days = time_rel / 86400;    % days since detection

    % 1) Prefit vs postfit range-rate residuals
    figure;
    subplot(2,1,1);
    plot(t_days, res_pre*1e3, '.-'); hold on;
    plot(t_days, res_post*1e3, '.-');
    xlabel('Days since detection');
    ylabel('Range-rate residual [m/s]');
    legend('Prefit','Postfit');
    grid on;
    title('Range-rate residuals (0–6 days)');

    % 2) k_SRP estimate + 3-sigma
    k_est   = x_hist(7,:);
    k_sig   = sig_hist(7,:);
    subplot(2,1,2);
    plot(t_days, k_est, 'LineWidth', 1.5); hold on;
    plot(t_days, k_est + 3*k_sig, '--');
    plot(t_days, k_est - 3*k_sig, '--');
    xlabel('Days since detection');
    ylabel('k_{SRP}');
    legend('estimate','\pm 3\sigma');
    grid on;
    title('SRP scale factor estimate');

    % 3) Antarctica lat/lon + 3-sigma
    figure;
    subplot(2,1,1);
    lat_est = x_hist(8,:);
    lat_sig = sig_hist(8,:);
    plot(t_days, lat_est*180/pi, 'LineWidth', 1.5); hold on;
    plot(t_days, (lat_est + 3*lat_sig)*180/pi, '--');
    plot(t_days, (lat_est - 3*lat_sig)*180/pi, '--');
    xlabel('Days since detection');
    ylabel('Latitude_4 [deg]');
    legend('estimate','\pm 3\sigma');
    grid on; title('Antarctica latitude estimate');

    subplot(2,1,2);
    lon_est = x_hist(9,:);
    lon_sig = sig_hist(9,:);
    plot(t_days, lon_est*180/pi, 'LineWidth', 1.5); hold on;
    plot(t_days, (lon_est + 3*lon_sig)*180/pi, '--');
    plot(t_days, (lon_est - 3*lon_sig)*180/pi, '--');
    xlabel('Days since detection');
    ylabel('Longitude_4 [deg]');
    legend('estimate','\pm 3\sigma');
    grid on; title('Antarctica longitude estimate');

    % 4) Range-rate bias + 3-sigma
    figure;
    bias_est = x_hist(10,:);
    bias_sig = sig_hist(10,:);
    plot(t_days, bias_est*1e3, 'LineWidth', 1.5); hold on;
    plot(t_days, (bias_est + 3*bias_sig)*1e3, '--');
    plot(t_days, (bias_est - 3*bias_sig)*1e3, '--');
    xlabel('Days since detection');
    ylabel('Bias_{rr} [m/s]');
    legend('estimate','\pm 3\sigma');
    grid on; title('Range-rate bias estimate');
end
