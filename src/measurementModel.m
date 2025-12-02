function [y, H] = measurementModel(T, X, stationID, params, st)
% measurementModel
%   Computes predicted range and range-rate measurements and
%   the corresponding Jacobian H wrt the 10x1 state X.
%
% State X layout:
%   1:3   r_sc      spacecraft position (Sun-centered EMO, km)
%   4:6   v_sc      spacecraft velocity (Sun-centered EMO, km/s)
%   7     k_SRP     SRP scale factor (unitless)
%   8     lat_4     Antarctica station latitude [rad]
%   9     lon_4     Antarctica station longitude [rad]
%   10    bias_rr   range-rate bias [km/s]
%
% Inputs:
%   T         - ephemeris time (sec past J2000)
%   X         - 10x1 state vector
%   stationID - integer 1..4 (1-3 DSN, 4 Antarctica)
%   params    - params struct from projectInit
%   st        - station struct from projectInit (lat, long arrays)
%
% Outputs:
%   y - 2x1 [range; rangeRate]
%   H - 2x10 measurement Jacobian wrt X

    % Unpack spacecraft state
    r_sc    = X(1:3);    % km
    v_sc    = X(4:6);    % km/s
    k_SRP   = X(7);      %#ok<NASGU>  % no direct dependence in measurement
    lat_4   = X(8);      % rad (Antarctica)
    lon_4   = X(9);      % rad
    bias_rr = X(10);     % km/s

    % Get Earth state in Sun-centered EMO (using SPICE; must have kernels loaded)
    xEarth  = cspice_spkezr('Earth', T, 'ECLIPJ2000', 'NONE', 'Sun');
    r_earth = xEarth(1:3);   % km
    v_earth = xEarth(4:6);   % km/s

    % Select station latitude/longitude
    if stationID == 4
        % Antarctica station uses estimated lat/lon from state
        lat_sta = lat_4;
        lon_sta = lon_4;
    else
        % DSN stations: fixed, truth lat/lon from st
        lat_sta = st.lat(stationID);
        lon_sta = st.long(stationID);
    end

    % Station position/velocity relative to Earth center in EMO frame
    [r_loc, v_loc, dr_loc_dlat, dr_loc_dlon, dv_loc_dlat, dv_loc_dlon] = ...
        stationInertialState(T, lat_sta, lon_sta, params);

    % Station inertial state in Sun-centered EMO
    r_sta = r_earth + r_loc;  % km
    v_sta = v_earth + v_loc;  % km/s

    % Relative geometry
    rho_vec = r_sc - r_sta;   % LOS vector (km)
    rho     = norm(rho_vec);  % range (km)
    u       = rho_vec / rho;  % LOS unit vector

    v_rel   = v_sc - v_sta;   % relative velocity (km/s)

    % Predicted measurements
    range     = rho;
    rangeRate = dot(v_rel, u) + bias_rr;  % km/s

    y = [range; rangeRate];

    % ------------------ Build H (2x10) ------------------ %
    H = zeros(2,10);
    I3 = eye(3);

    % Range partials:
    %   h1 = rho = ||rho_vec||
    %   ∂h1/∂r_sc = u'
    H(1,1:3) = u';       % wrt r_sc
    % ∂h1/∂v_sc = 0
    % no direct sensitivity to k_SRP, bias, etc. in measurement function

    % Range-rate partials:
    %   h2 = v_rel · u + bias_rr
    %   ∂h2/∂r_sc = (v_rel'/rho) * (I - u*u')
    dhr2_dr = (v_rel'/rho) * (I3 - u*u');    % 1x3
    H(2,1:3) = dhr2_dr;                      % wrt r_sc

    %   ∂h1/∂v_sc = 0,  ∂h2/∂v_sc = u'
    H(2,4:6) = u';                           % wrt v_sc

    %   ∂h1/∂k_SRP = 0, ∂h2/∂k_SRP = 0 (only via dynamics)
    %   so H(:,7) = 0 is fine.

    %   ∂h2/∂bias_rr = 1
    H(2,10) = 1;

    % ------------------ Partials wrt Antarctica lat/lon ------------------ %
    % Only apply for stationID == 4; for other stations those columns are zero.
    if stationID == 4
        % rho_vec = r_sc - (r_earth + r_loc)
        % ∂rho_vec/∂lat = -∂r_loc/∂lat,   ∂rho_vec/∂lon = -∂r_loc/∂lon
        drho_dlat_vec = -dr_loc_dlat;
        drho_dlon_vec = -dr_loc_dlon;

        % Range:
        %   ∂rho/∂p = u' * ∂rho_vec/∂p
        drange_dlat = u' * drho_dlat_vec;
        drange_dlon = u' * drho_dlon_vec;

        H(1,8) = drange_dlat;
        H(1,9) = drange_dlon;

        % Range-rate:
        %   v_rel = v_sc - (v_earth + v_loc)
        %   ∂v_rel/∂p = -∂v_loc/∂p
        dvrel_dlat = -dv_loc_dlat;
        dvrel_dlon = -dv_loc_dlon;

        %   ∂u/∂p = (I - u*u') * (∂rho_vec/∂p) / rho
        du_dlat = (I3 - u*u') * drho_dlat_vec / rho;
        du_dlon = (I3 - u*u') * drho_dlon_vec / rho;

        %   ∂h2/∂p = (∂v_rel/∂p)' * u + v_rel' * ∂u/∂p
        drate_dlat = dvrel_dlat' * u + v_rel' * du_dlat;
        drate_dlon = dvrel_dlon' * u + v_rel' * du_dlon;

        H(2,8) = drate_dlat;
        H(2,9) = drate_dlon;
    end
end
