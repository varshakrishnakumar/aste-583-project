function out = runBatchLSQ(span, mode, opts)
% runBatchLSQ  Batch Gauss-Newton LSQ OD for ASTE583.
%
% State at ET0:
%   X = [r(3); v(3); k_SRP; lat4; lon4; bias_rr]
%
% span: '0to6' | '6to14' | '0to14' | filename | string array of filenames
%
% mode:
%   'handout'          (default) : use projectInit() prior, estimate k freely (can go negative)
%   'handout_logk'                : enforce k>0 using eta=log(k) update (no dynamics changes)
%   'ekf_k'                       : take k from EKF and hold k fixed (column 7 = 0)
%   'ekf_prior'                   : use EKF k as prior mean, estimate k (physical coords)
%   'ekf_prior_logk'              : same but estimate eta=log(k) (enforce k>0)
%
% opts (optional struct):
%   opts.maxIter    (default 8)
%   opts.tol_dx     (default 1e-6)
%   opts.useLM      (default false)  % Levenberg-Marquardt damping on normal matrix
%   opts.lmLambda   (default 1e-6)  % starting damping
%   opts.robust     (default false) % Huber IRLS
%   opts.huberC     (default 3.0)   % threshold in sigmas
%   opts.earlyMode  = "use"|"downweight"|"drop" (default "use")
%   opts.earlyDays  (default 2)
%   opts.earlyScale (default 10)   % sigma multiplier during early window if downweight
%   opts.makePlots  (default true)
%
% Outputs:
%   out.X_est, out.P_est at ET0 (physical units)
%   out.iterHist
%   out.res_prefit, out.res_postfit (structs with .rho, .rr)
%   out.tMeas, out.stationIDs, out.files, out.ET0

    if nargin < 1 || isempty(span), span = '0to14'; end
    if nargin < 2 || isempty(mode), mode = 'handout'; end
    if nargin < 3, opts = struct(); end

    span = string(span);
    mode = lower(string(mode));

    % ---------- Defaults ----------
    if ~isfield(opts,'maxIter'),    opts.maxIter = 8; end
    if ~isfield(opts,'tol_dx'),     opts.tol_dx  = 1e-6; end
    if ~isfield(opts,'useLM'),      opts.useLM   = false; end
    if ~isfield(opts,'lmLambda'),   opts.lmLambda= 1e-6; end
    if ~isfield(opts,'robust'),     opts.robust  = false; end
    if ~isfield(opts,'huberC'),     opts.huberC  = 3.0; end
    if ~isfield(opts,'earlyMode'),  opts.earlyMode = "use"; end
    if ~isfield(opts,'earlyDays'),  opts.earlyDays = 2; end
    if ~isfield(opts,'earlyScale'), opts.earlyScale= 10; end
    % if ~isfield(opts,'makePlots'),  opts.makePlots = true; end
    % if ~isfield(opts,'saveFigures'), opts.saveFigures = false; end
    % if ~isfield(opts,'outDir')
    %     opts.outDir = fullfile(pwd, sprintf("figs_Batch_%s_%s", span, mode));
    % end
    % if ~isfield(opts,'figPrefix')
    %     opts.figPrefix = sprintf("Batch_%s_%s", span, mode);
    % end


    useLogK     = (mode == "handout_logk" || mode == "ekf_prior_logk");
    freezeK     = (mode == "ekf_k");   % only this mode truly holds k fixed
    useEkfPrior = (mode == "ekf_k" || mode == "ekf_prior" || mode == "ekf_prior_logk");

    % ---------- Init ----------
    [params, sc, st, X0_phys, P0_phys] = projectInit();
    [ET0, ~, ~] = initTime();

    params.t0_et       = ET0;
    params.bias_end_et = ET0 + 6*params.tday;

    % ---------- Optionally get k_SRP from EKF ----------
    if useEkfPrior
        try
            fprintf('runBatchLSQ: getting k_SRP from EKF (recommend 6to14)...\n');
            ekf   = runPrelimEKF([], '0to6');   % less bias/outgas contamination than 0to6
            k_ekf = ekf.x_final(7);
            fprintf('  EKF k_SRP = %.6f\n', k_ekf);

            % EKF value used as prior mean for k
            X0_phys(7) = k_ekf;

            % Prior variance on k (e.g., 1σ ≈ 0.3 on k)
            P0_phys(7,7) = max(P0_phys(7,7), (0.3)^2);
        catch ME
            warning('runBatchLSQ:EKFfail','EKF k fetch failed: %s. Using handout prior k instead. (%s)', ME.message);
            if freezeK
                % If EKF fails, "ekf_k" mode falls back to a handout-style estimate.
                freezeK = false;
            end
        end
    end

    % ---------- Prior (possibly in eta=log(k) coordinates) ----------
    if useLogK
        if X0_phys(7) <= 0
            error('runBatchLSQ:BadPriorK','Need positive prior k to use log-k mode.');
        end
        X_apri = X0_phys;
        X_apri(7) = log(X0_phys(7));                    % eta prior
        P_apri = P0_phys;
        P_apri(7,7) = P0_phys(7,7) / (X0_phys(7)^2);    % var(eta) ≈ var(k)/k^2
    else
        X_apri = X0_phys;
        P_apri = P0_phys;
    end

    % Prior whitening P^{-1/2} (with jitter safety)
    [W0, P_apri] = safeWhiten(P_apri); %#ok<NASGU>

    % ---------- Load measurements ----------
    files = resolveFiles(span);

    measAll.t = []; measAll.stationID = []; measAll.range = []; measAll.rr = [];
    for f = files
        m = parseMeasurementData(f, ET0);
        measAll.t         = [measAll.t;         m.t(:)];
        measAll.stationID = [measAll.stationID; m.stationID(:)];
        measAll.range     = [measAll.range;     m.range(:)];
        measAll.rr        = [measAll.rr;        m.rr(:)];
    end

    [tMeas, ord] = sort(measAll.t(:));
    stationIDs   = measAll.stationID(ord);
    obs_rho      = measAll.range(ord);
    obs_rr       = measAll.rr(ord);
    N            = numel(tMeas);

    % Early window handling
    earlyCut = ET0 + opts.earlyDays*86400;
    if opts.earlyMode == "drop"
        keep      = (tMeas > earlyCut);
        tMeas     = tMeas(keep);
        stationIDs= stationIDs(keep);
        obs_rho   = obs_rho(keep);
        obs_rr    = obs_rr(keep);
        N         = numel(tMeas);
    end

    [tU, ~, idxU] = unique(tMeas, 'stable');

    % ---------- Noise model ----------
    sig.rho_dsn = 1e-3;   % km
    sig.rho_ant = 1e-2;   % km
    sig.rr_dsn  = 1e-7;   % km/s
    sig.rr_ant  = 1e-6;   % km/s

    scale.rho_dsn = 5;
    scale.rho_ant = 5;
    scale.rr_dsn  = 5;
    scale.rr_ant  = 5;

    maxIter = opts.maxIter;
    tol_dx  = opts.tol_dx;
    odeOpts = odeset('RelTol',1e-11,'AbsTol',1e-13);

    odeFun = @(t,X) projectDynamics(t, X, params, sc);

    % X_est internal coordinates: either physical k or eta=log(k)
    X_est = X_apri;

    iterHist = zeros(maxIter,5);  % [it, ||dX||, cost, RR_RMS(mm/s), Rho_RMS(m)]
    fprintf('Batch LSQ using %s (%d measurements), mode=%s\n', strjoin(files,", "), N, mode);
    fprintf('%4s  %12s  %12s  %14s  %14s\n','it','||dX||','cost','RR_RMS [mm/s]','Rho_RMS [m]');

    % We'll also keep the final H/r/sigma from the last iteration
    H_all = []; r_all = []; sigma_all = [];

    for it = 1:maxIter

        % ---- propagate from ET0 to all unique times ----
        X_prop0 = X_est(:);
        if useLogK
            % convert eta -> k for propagation only
            X_prop0(7) = exp(X_est(7));
        end

        X_aug0 = [X_prop0; reshape(eye(10),100,1)];

        if abs(tU(1) - ET0) < 1e-6
            tspan     = tU;
            dropFirst = false;
        else
            tspan     = [ET0; tU];
            dropFirst = true;
        end

        [~, Xaug] = ode113(odeFun, tspan, X_aug0, odeOpts);
        if dropFirst, Xaug = Xaug(2:end,:); end

        XU   = Xaug(:,1:10);     % physical states (k in slot 7)
        PhiU = Xaug(:,11:end);   % STM w.r.t. physical initial state at ET0

        k0_phys = X_prop0(7);    % physical k at ET0 used this iteration

        % ---- build H and residual vector r ----
        nRange = sum(isfinite(obs_rho));
        M      = N + nRange;

        H_all     = zeros(M,10);
        r_all     = zeros(M,1);
        sigma_all = zeros(M,1);
        rowType   = zeros(M,1); % 1=rho, 2=rr

        row = 0;

        for k = 1:N
            ui  = idxU(k);

            xk_phys = XU(ui,:).';
            Phi     = reshape(PhiU(ui,:).',10,10);

            stID = stationIDs(k);

            % Bias only active first 6 days
            biasActive = (tMeas(k) <= params.bias_end_et);
            x_meas = xk_phys;
            if ~biasActive
                x_meas(10) = 0;
            end

            [y_pred, Hloc] = measurementModel(tMeas(k), x_meas, stID, params, st);
            if ~biasActive
                Hloc(2,10) = 0;
            end

            % Early downweighting option
            isEarly = (tMeas(k) <= earlyCut);

            % RANGE
            if isfinite(obs_rho(k))
                row         = row + 1;
                r_all(row)  = obs_rho(k) - y_pred(1);

                H_row = Hloc(1,:) * Phi;

                % convert column 7 from dk to deta if using logk
                if useLogK
                    H_row(7) = H_row(7) * k0_phys;
                end

                % freeze k if requested
                if freezeK
                    H_row(7) = 0;
                end

                H_all(row,:)   = H_row;
                sigma_all(row) = measSigma(stID,"rho",sig,scale,isEarly,opts);
                rowType(row)   = 1;
            end

            % RANGE-RATE
            row        = row + 1;
            r_all(row) = obs_rr(k) - y_pred(2);

            H_row = Hloc(2,:) * Phi;

            if useLogK
                H_row(7) = H_row(7) * k0_phys;
            end
            if freezeK
                H_row(7) = 0;
            end

            H_all(row,:)   = H_row;
            sigma_all(row) = measSigma(stID,"rr",sig,scale,isEarly,opts);
            rowType(row)   = 2;   % row type flag for range-rate residual
        end

        % trim
        H_all     = H_all(1:row,:);
        r_all     = r_all(1:row);
        sigma_all = sigma_all(1:row);
        rowType   = rowType(1:row);

        % ---- weights ----
        w = 1 ./ sigma_all;

        % Optional robust IRLS (Huber)
        if opts.robust
            e = r_all .* w;  % normalized residuals
            w = w .* huberSqrtWeight(e, opts.huberC);
        end

        Hw = H_all .* w;
        rw = r_all .* w;

        % add prior rows
        Hw_tot = [Hw; W0];
        rw_tot = [rw; W0*(X_apri - X_est)];

        % ---- solve ----
        if opts.useLM
            % LM on normal matrix (10x10)
            Nmat = Hw_tot.'*Hw_tot;
            g    = Hw_tot.'*rw_tot;
            D    = diag(max(diag(Nmat),1e-12));
            dX   = (Nmat + opts.lmLambda*D) \ g;
        else
            [Q,R] = qr(Hw_tot,0);
            dX    = R \ (Q' * rw_tot);
        end

        % ---- update ----
        X_est = X_est + dX;

        % diagnostics
        cost = rw_tot.' * rw_tot;

        rr_rms_mmps  = sqrt(mean(r_all(rowType==2).^2)) * 1e6;
        if any(rowType==1)
            rho_rms_m = sqrt(mean(r_all(rowType==1).^2)) * 1e3;
        else
            rho_rms_m = NaN;
        end

        iterHist(it,:) = [it, norm(dX), cost, rr_rms_mmps, rho_rms_m];
        fprintf('%4d  %12.3e  %12.3e  %14.3f  %14.3f\n', it, norm(dX), cost, rr_rms_mmps, rho_rms_m);

        if norm(dX) < tol_dx
            iterHist = iterHist(1:it,:);
            break;
        end
    end

    % ---------- Final covariance (posterior, internal coords) ----------
    if opts.useLM
        % Rebuild weights at final iterate using last H_all/r_all/sigma_all
        w = 1 ./ sigma_all;
        if opts.robust
            e = r_all .* w;
            w = w .* huberSqrtWeight(e, opts.huberC);
        end
        Hw   = H_all .* w;
        Nmat = Hw.'*Hw + (W0.'*W0);
        P_int = invSymPD(Nmat);
    else
        % Use R from last QR (10x10 upper-triangular)
        Ri    = R \ eye(10);
        P_int = Ri * Ri.';
        P_int = 0.5*(P_int + P_int.');
    end

    % Convert internal eta-cov to physical k-cov if logk mode
    X_final = X_est;
    P_final = P_int;
    if useLogK
        eta    = X_est(7);
        k_phys = exp(eta);
        X_final(7) = k_phys;

        J       = eye(10);
        J(7,7)  = k_phys;          % dk = k * deta
        P_final = J * P_int * J.';
        P_final = 0.5*(P_final + P_final.');
    end

    % ---------- Residuals pre/post ----------
    res_prefit  = computeResiduals(X_apri);
    res_postfit = computeResiduals(X_est);

    % ---------- Print final ----------
    sig3 = 3*sqrt(max(diag(P_final),0));
    labels = {'x [km]','y [km]','z [km]','vx [km/s]','vy [km/s]','vz [km/s]',...
              'k_SRP','lat_4 [rad]','lon_4 [rad]','bias_rr [km/s]'};
    fprintf('\nFinal Batch LSQ estimate at ET0 (mode=%s):\n', mode);
    for i=1:10
        fprintf('%-12s  %.6e   (± %.6e 3σ)\n', labels{i}, X_final(i), sig3(i));
    end
    % 
    % % ---------- Plots (optional) ----------
    % if opts.makePlots
    %     makeBatchPlots();
    % end
    % 
    %     % ---------- Save figures (optional) ----------
    % if opts.saveFigures
    %     if ~exist(opts.outDir,'dir')
    %         mkdir(opts.outDir);
    %     end
    %     saveAllFigures(opts.outDir, char(opts.figPrefix));
    % end


    % ---------- outputs ----------
    out.X_est      = X_final;
    out.P_est      = P_final;
    out.iterHist   = iterHist;
    out.files      = files;
    out.ET0        = ET0;
    out.res_prefit = res_prefit;
    out.res_postfit= res_postfit;
    out.tMeas      = tMeas;
    out.stationIDs = stationIDs;
    out.X0 = X0_phys;   % prior mean at ET0
    out.P0 = P0_phys;   % prior covariance (optional but nice)
    out.span = span;   % remember the data span: '0to6','6to14','0to14',...



% ================= nested helpers =================

    function files = resolveFiles(spec)
        spec = string(spec);
        if numel(spec)>1, files = spec; return; end
        s = lower(strtrim(spec));
        if s=="0to6" || s=="0-6"
            files = "ASTE583_Project_LTB_Measurements_0-6D_Truth.csv";
        elseif s=="6to14" || s=="6-14"
            files = "ASTE583_Project_LTB_Measurements_6D-14D_Truth.csv";
        elseif s=="0to14" || s=="0-14" || s=="all"
            files = ["ASTE583_Project_LTB_Measurements_0-6D_Truth.csv", ...
                     "ASTE583_Project_LTB_Measurements_6D-14D_Truth.csv"];
        else
            files = spec;
        end
    end


    function res = computeResiduals(X_epoch_int)
        % X_epoch_int is in internal coordinates (k or eta).
        X0prop = X_epoch_int(:);
        if useLogK
            X0prop(7) = exp(X_epoch_int(7));
        end

        % Propagate state only (STM not used here; augmented identity is reused
        % for compatibility with projectDynamics)
        X_aug0 = [X0prop; reshape(eye(10),100,1)];

        if abs(tU(1)-ET0)<1e-6
            tspanR     = tU; 
            dropFirstR = false;
        else
            tspanR     = [ET0; tU]; 
            dropFirstR = true;
        end

        [~, XaugR] = ode113(odeFun, tspanR, X_aug0, odeOpts);
        if dropFirstR, XaugR = XaugR(2:end,:); end

        XU_r = XaugR(:,1:10);

        res.rho = nan(N,1);
        res.rr  = nan(N,1);

        for kk=1:N
            xk = XU_r(idxU(kk),:).';

            biasActive = (tMeas(kk) <= params.bias_end_et);
            if ~biasActive
                xk(10) = 0;
            end

            y_pred = measurementModel(tMeas(kk), xk, stationIDs(kk), params, st);

            if isfinite(obs_rho(kk))
                res.rho(kk) = obs_rho(kk) - y_pred(1);
            end
            res.rr(kk) = obs_rr(kk) - y_pred(2);
        end
    end

end

% ================= helper functions outside =================

function [W0, Pspd] = safeWhiten(P)
    Pspd  = 0.5*(P+P.');
    n     = size(Pspd,1);
    jitter= 0;
    for k=1:8
        [L,p] = chol(Pspd,'lower');
        if p==0
            W0 = L \ eye(n);
            return;
        end
        jitter = max(1e-18, 10^k * 1e-18);
        Pspd   = Pspd + jitter*eye(n);
    end
    error('safeWhiten:SPD','Could not Cholesky-factor prior covariance.');
end

function s = measSigma(stID, type, sig, scale, isEarly, opts)
    if stID == 4
        if type=="rho", s = sig.rho_ant*scale.rho_ant;
        else,           s = sig.rr_ant *scale.rr_ant;
        end
    else
        if type=="rho", s = sig.rho_dsn*scale.rho_dsn;
        else,           s = sig.rr_dsn *scale.rr_dsn;
        end
    end

    if opts.earlyMode == "downweight" && isEarly
        s = s * opts.earlyScale;
    end
end

function w = huberSqrtWeight(e, c)
    a = abs(e);
    w = ones(size(e));
    ii = a > c;
    w(ii) = sqrt(c ./ a(ii));
end

function P = invSymPD(N)
    % Small 10x10 symmetric PD inversion via Cholesky solve
    N = 0.5*(N+N.');
    R = chol(N,'lower');
    Ri= R \ eye(size(N,1));
    P = Ri.' * Ri;
    P = 0.5*(P+P.');
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
    % Plot nsig covariance ellipse for 2x2 covariance P2 around mean mu
    P2 = 0.5*(P2+P2.');
    [V,D] = eig(P2);
    D     = max(D, 0);
    A     = V*sqrt(D);

    th   = linspace(0, 2*pi, 200);
    circ = [cos(th); sin(th)];
    ell  = mu(:) + nsig * (A * circ);

    plot(ell(1,:), ell(2,:), 'k--', 'LineWidth',1.2, ...
         'DisplayName', sprintf('%g\\sigma ellipse', nsig));
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
