function [y, H] = measurementModel(t, X, stationID, params, st, measType)
    if nargin < 6 || isempty(measType)
        measType = 'both';
    end
    measType = lower(string(measType));
    if stationID < 1 || stationID > 4
        error('measurementModel:BadStationID', 'stationID must be 1..4.');
    end
    X = X(:);
    if numel(X) < 10
        error('measurementModel:BadState', 'X must have at least 10 elements.');
    end
    x = X(1:10);
    r_sc    = x(1:3);
    v_sc    = x(4:6);
    lat_4   = x(8);
    lon_4   = x(9);
    bias_rr = x(10);
    xEarth  = cspice_spkezr('EARTH', t, 'ECLIPJ2000', 'NONE', 'SUN');
    r_earth = xEarth(1:3);
    v_earth = xEarth(4:6);
    if stationID == 4
        lat_sta = lat_4;
        lon_sta = lon_4;
    else
        lat_sta = st.lat(stationID);
        lon_sta = st.long(stationID);
    end
    if stationID == 4
        [r_loc, v_loc, dr_dlat, dr_dlon, dv_dlat, dv_dlon] = ...
            stationInertialState(t, lat_sta, lon_sta, params);
    else
        [r_loc, v_loc] = stationInertialState(t, lat_sta, lon_sta, params);
        dr_dlat = []; dr_dlon = []; dv_dlat = []; dv_dlon = [];
    end
    r_sta = r_earth + r_loc;
    v_sta = v_earth + v_loc;
    rho_vec = r_sc - r_sta;
    rho     = norm(rho_vec);
    u       = rho_vec / rho;
    v_rel   = v_sc - v_sta;
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
    Hfull = zeros(2,10);
    I3 = eye(3);
    Hfull(1,1:3) = u.';     % d(rho)/d(r_sc)
    Hfull(2,4:6) = u.';     % d(rr)/d(v_sc)
    Hfull(2,10)  = drr_db;
    Hfull(2,1:3) = (v_rel.' / rho) * (I3 - u*u.');
    if stationID == 4
        drho_dlat = -dr_dlat;
        drho_dlon = -dr_dlon;
        Hfull(1,8) = u.' * drho_dlat;
        Hfull(1,9) = u.' * drho_dlon;
        dvrel_dlat = -dv_dlat;
        dvrel_dlon = -dv_dlon;
        du_dlat = (I3 - u*u.') * drho_dlat / rho;
        du_dlon = (I3 - u*u.') * drho_dlon / rho;
        Hfull(2,8) = dvrel_dlat.' * u + v_rel.' * du_dlat;
        Hfull(2,9) = dvrel_dlon.' * u + v_rel.' * du_dlon;
    end
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