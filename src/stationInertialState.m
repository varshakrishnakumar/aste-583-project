function [r_loc, v_loc, dr_loc_dlat, dr_loc_dlon, dv_loc_dlat, dv_loc_dlon] = ...
    stationInertialState(T, lat, lon, params)
% stationInertialState
%   Computes station position/velocity in Sun-centered EMO frame,
%   relative to Earth's center, and partials wrt lat and lon.
%
% Inputs:
%   T     - ephemeris time (sec past J2000)
%   lat   - station latitude [rad]
%   lon   - station longitude [rad]
%   params - struct from projectInit (Re, we, GST0, R_EMEtoEMO)
%
% Outputs (all 3x1 vectors):
%   r_loc       - station position in EMO2000, Earth-centered (km)
%   v_loc       - station velocity in EMO2000, Earth-centered (km/s)
%   dr_loc_dlat - ∂r_loc/∂lat
%   dr_loc_dlon - ∂r_loc/∂lon
%   dv_loc_dlat - ∂v_loc/∂lat
%   dv_loc_dlon - ∂v_loc/∂lon

    Re   = params.Re;
    we   = params.we;
    GST0 = params.GST0;
    R_EMEtoEMO = params.R_EMEtoEMO;

    % Earth rotation angle about z at time T (inertial)
    theta = GST0 + we*T;
    cth = cos(theta);
    sth = sin(theta);
    R3 = [ cth -sth 0;
           sth  cth 0;
            0    0  1 ];

    % Station position in Earth-fixed (ECEF)
    clat = cos(lat); slat = sin(lat);
    clon = cos(lon); slon = sin(lon);

    r_ECEF = Re * [ clat*clon;
                    clat*slon;
                    slat      ];

    % Partials wrt lat, lon in ECEF
    dr_ECEF_dlat = Re * [ -slat*clon;
                          -slat*slon;
                           clat     ];

    dr_ECEF_dlon = Re * [ -clat*slon;
                           clat*clon;
                           0         ];

    % Rotate to ECI (EME2000)
    r_ECI       = R3 * r_ECEF;
    dr_ECI_dlat = R3 * dr_ECEF_dlat;
    dr_ECI_dlon = R3 * dr_ECEF_dlon;

    % Transform to EMO (ecliptic) frame
    r_loc       = R_EMEtoEMO * r_ECI;
    dr_loc_dlat = R_EMEtoEMO * dr_ECI_dlat;
    dr_loc_dlon = R_EMEtoEMO * dr_ECI_dlon;

    % Velocity from Earth's rotation: v_ECI = Ω × r_ECI
    Omega_eci   = [0; 0; we];
    v_ECI       = cross(Omega_eci, r_ECI);
    dv_ECI_dlat = cross(Omega_eci, dr_ECI_dlat);
    dv_ECI_dlon = cross(Omega_eci, dr_ECI_dlon);

    v_loc       = R_EMEtoEMO * v_ECI;
    dv_loc_dlat = R_EMEtoEMO * dv_ECI_dlat;
    dv_loc_dlon = R_EMEtoEMO * dv_ECI_dlon;
end
