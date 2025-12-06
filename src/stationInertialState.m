function [r_loc, v_loc, dr_loc_dlat, dr_loc_dlon, dv_loc_dlat, dv_loc_dlon] = ...
    stationInertialState(t_et, lat, lon, params)
%STATIONINERTIALSTATE Station state wrt Earth center in EMO (ecliptic J2000).
%
% Inputs:
%   t_et  - ephemeris time (ET seconds past J2000), SPICE-compatible
%   lat   - geodetic latitude [rad]  (treated as geocentric here; altitude=0)
%   lon   - east longitude [rad]     (positive-east)
%   params.Re        - Earth radius [km]
%   params.we        - Earth rotation rate [rad/s]
%   params.GST0      - phi_G0 at detection epoch [rad] (00:10:43)
%   params.t0_et     - detection epoch ET seconds past J2000
%   params.R_EMEtoEMO- rotation EME2000 -> EMO2000
%
% Outputs (3x1):
%   r_loc, v_loc         - station position/velocity in EMO, Earth-centered
%   dr_loc_dlat/dlon     - partials of r_loc wrt lat/lon
%   dv_loc_dlat/dlon     - partials of v_loc wrt lat/lon

    % ---- Required params ----
    if ~isfield(params,'t0_et')
        error('stationInertialState:MissingT0', ...
              'params.t0_et is required (ET at detection epoch) to use GST0 correctly.');
    end

    Re   = params.Re;
    we   = params.we;
    GST0 = params.GST0;
    R_EMEtoEMO = params.R_EMEtoEMO;

    % ---- Earth rotation angle about z (anchored at detection epoch) ----
    dt = t_et - params.t0_et;                 % seconds since detection
    theta = GST0 + we * dt;                   % rad
    theta = mod(theta, 2*pi);                 % keep trig numerically happy

    cth = cos(theta);
    sth = sin(theta);

    % ECEF -> ECI (EME2000) rotation about z
    R3 = [ cth -sth 0;
           sth  cth 0;
            0    0  1 ];

    % ---- Station position in ECEF (spherical Earth, altitude=0) ----
    clat = cos(lat); slat = sin(lat);
    clon = cos(lon); slon = sin(lon);

    r_ecef = Re * [ clat*clon;
                    clat*slon;
                    slat ];

    % Partials in ECEF
    if nargout > 2
        dr_ecef_dlat = Re * [ -slat*clon;
                              -slat*slon;
                               clat ];

        dr_ecef_dlon = Re * [ -clat*slon;
                               clat*clon;
                               0 ];
    end

    % ---- Rotate to ECI (EME2000) ----
    r_eci = R3 * r_ecef;

    if nargout > 2
        dr_eci_dlat = R3 * dr_ecef_dlat;
        dr_eci_dlon = R3 * dr_ecef_dlon;
    end

    % ---- Velocity from Earth rotation: v = Ω × r ----
    Omega = [0; 0; we];
    v_eci = cross(Omega, r_eci);

    if nargout > 2
        dv_eci_dlat = cross(Omega, dr_eci_dlat);
        dv_eci_dlon = cross(Omega, dr_eci_dlon);
    end

    % ---- Convert to EMO (ecliptic) frame ----
    r_loc = R_EMEtoEMO * r_eci;
    v_loc = R_EMEtoEMO * v_eci;

    if nargout > 2
        dr_loc_dlat = R_EMEtoEMO * dr_eci_dlat;
        dr_loc_dlon = R_EMEtoEMO * dr_eci_dlon;
        dv_loc_dlat = R_EMEtoEMO * dv_eci_dlat;
        dv_loc_dlon = R_EMEtoEMO * dv_eci_dlon;
    end
end
