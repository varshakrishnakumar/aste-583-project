function [y, H] = measurementModel(t, X, stationID, params, st, measType)
%MEASUREMENTMODEL Predict range / range-rate and measurement Jacobian.
%
% State (uses first 10 elements of X; X may optionally include an STM):
%   1:3   r_sc      [km]    Spacecraft position wrt Sun (ECLIPJ2000 / EMO)
%   4:6   v_sc      [km/s]  Spacecraft velocity wrt Sun
%   7     k_SRP     [-]     SRP scale factor (no direct measurement sensitivity)
%   8     lat_4     [rad]   Antarctica station latitude
%   9     lon_4     [rad]   Antarctica station longitude (positive-east)
%   10    bias_rr   [km/s]  range-rate bias (only active first 6 days)
%
% Inputs:
%   t         - ET seconds past J2000 (SPICE time)
%   X         - state (10x1 or 110x1; first 10 used)
%   stationID - 1..4
%   params    - constants struct (needs mue, tday (optional), t0_et (optional))
%   st        - station struct with .lat, .long for stations 1..3
%   measType  - (optional) 'both' (default) | 'range' | 'rr'
%
% Outputs:
%   y - predicted measurement(s)
%       'both'  -> [range; rr]  (2x1)
%       'range' -> range        (1x1)
%       'rr'    -> rr           (1x1)
%   H - Jacobian wrt 10-state (rows match y)

    if nargin < 6 || isempty(measType)
        measType = 'both';
    end
    measType = lower(string(measType));

    % Basic checks
    if stationID < 1 || stationID > 4
        error('measurementModel:BadStationID', 'stationID must be 1..4.');
    end

    X = X(:);
    if numel(X) < 10
        error('measurementModel:BadState', 'X must have at least 10 elements.');
    end
    x = X(1:10);

    % Unpack state
    r_sc    = x(1:3);
    v_sc    = x(4:6);
    lat_4   = x(8);
    lon_4   = x(9);
    bias_rr = x(10);

    % Earth state wrt Sun in ECLIPJ2000 (EMO2000)
    xEarth  = cspice_spkezr('EARTH', t, 'ECLIPJ2000', 'NONE', 'SUN');
    r_earth = xEarth(1:3);
    v_earth = xEarth(4:6);

    % Select station lat/lon
    if stationID == 4
        lat_sta = lat_4;
        lon_sta = lon_4;
    else
        lat_sta = st.lat(stationID);
        lon_sta = st.long(stationID);
    end

    % Get station state wrt Earth (and optionally partials wrt lat/lon)
    if stationID == 4
        [r_loc, v_loc, dr_dlat, dr_dlon, dv_dlat, dv_dlon] = ...
            stationInertialState(t, lat_sta, lon_sta, params);
    else
        [r_loc, v_loc] = stationInertialState(t, lat_sta, lon_sta, params);
        dr_dlat = []; dr_dlon = []; dv_dlat = []; dv_dlon = [];
    end

    % Station state wrt Sun
    r_sta = r_earth + r_loc;
    v_sta = v_earth + v_loc;

    % Line-of-sight geometry
    rho_vec = r_sc - r_sta;
    rho     = norm(rho_vec);
    u       = rho_vec / rho;     % LOS unit vector
    v_rel   = v_sc - v_sta;

    % Bias is only active during first 6 days after detection (if params supports it)
    biasActive = true;
    if isfield(params, 'bias_end_et')
        biasActive = (t <= params.bias_end_et);
    elseif isfield(params, 't0_et')
        tday = 86400;
        if isfield(params, 'tday'); tday = params.tday; end
        biasActive = (t <= params.t0_et + 6*tday);
    end

    if biasActive
        rr_pred = dot(v_rel, u) + bias_rr;
        drr_db  = 1;
    else
        rr_pred = dot(v_rel, u);
        drr_db  = 0;
    end

    range_pred = rho;

    % -------- Full measurement vector + Jacobian (2x10) -------- %
    Hfull = zeros(2,10);
    I3 = eye(3);

    % Range partials
    Hfull(1,1:3) = u.';     % d(rho)/d(r_sc)

    % Range-rate partials
    Hfull(2,4:6) = u.';     % d(rr)/d(v_sc)
    Hfull(2,10)  = drr_db;  % d(rr)/d(bias_rr) (0 after day 6)

    % d(rr)/d(r_sc) = v_rel' * d(u)/d(r_sc)
    % d(u)/d(r_sc) = (I - u u') / rho
    Hfull(2,1:3) = (v_rel.' / rho) * (I3 - u*u.');

    % Station 4 lat/lon partials
    if stationID == 4
        % rho_vec = r_sc - (r_earth + r_loc)  => d(rho_vec)/dp = -d(r_loc)/dp
        drho_dlat = -dr_dlat;
        drho_dlon = -dr_dlon;

        % Range: d(rho)/dp = u' * d(rho_vec)/dp
        Hfull(1,8) = u.' * drho_dlat;
        Hfull(1,9) = u.' * drho_dlon;

        % Range-rate:
        % v_rel = v_sc - (v_earth + v_loc) => d(v_rel)/dp = -d(v_loc)/dp
        dvrel_dlat = -dv_dlat;
        dvrel_dlon = -dv_dlon;

        % d(u)/dp = (I - u u') * d(rho_vec)/dp / rho
        du_dlat = (I3 - u*u.') * drho_dlat / rho;
        du_dlon = (I3 - u*u.') * drho_dlon / rho;

        % d(rr)/dp = d(v_rel)/dp · u + v_rel · d(u)/dp
        Hfull(2,8) = dvrel_dlat.' * u + v_rel.' * du_dlat;
        Hfull(2,9) = dvrel_dlon.' * u + v_rel.' * du_dlon;
    end

    % -------- Return requested measurement type -------- %
    switch measType
        case {"both","all"}
            y = [range_pred; rr_pred];
            H = Hfull;

        case {"range","rho"}
            y = range_pred;
            H = Hfull(1,:);

        case {"rr","rangerate","range-rate","rhodot"}
            y = rr_pred;
            H = Hfull(2,:);

        otherwise
            error('measurementModel:BadMeasType', ...
                'measType must be ''both'', ''range'', or ''rr''.');
    end
end
