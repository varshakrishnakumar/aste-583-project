function [r_loc, v_loc, dr_loc_dlat, dr_loc_dlon, dv_loc_dlat, dv_loc_dlon] = ...
    stationInertialState(t_et, lat, lon, params)
    if ~isfield(params,'t0_et')
        error('stationInertialState:MissingT0', ...
              'params.t0_et is required (ET at detection epoch) to use GST0 correctly.');
    end
    Re   = params.Re;
    we   = params.we;
    GST0 = params.GST0;
    R_EMEtoEMO = params.R_EMEtoEMO;
    dt = t_et - params.t0_et;
    theta = GST0 + we * dt;
    theta = mod(theta, 2*pi);
    cth = cos(theta);
    sth = sin(theta);
    R3 = [ cth -sth 0;
           sth  cth 0;
            0    0  1 ];
    clat = cos(lat); slat = sin(lat);
    clon = cos(lon); slon = sin(lon);
    r_ecef = Re * [ clat*clon;
                    clat*slon;
                    slat ];
    if nargout > 2
        dr_ecef_dlat = Re * [ -slat*clon;
                              -slat*slon;
                               clat ];
        dr_ecef_dlon = Re * [ -clat*slon;
                               clat*clon;
                               0 ];
    end
    r_eci = R3 * r_ecef;
    if nargout > 2
        dr_eci_dlat = R3 * dr_ecef_dlat;
        dr_eci_dlon = R3 * dr_ecef_dlon;
    end
    Omega = [0; 0; we];
    v_eci = cross(Omega, r_eci);
    if nargout > 2
        dv_eci_dlat = cross(Omega, dr_eci_dlat);
        dv_eci_dlon = cross(Omega, dr_eci_dlon);
    end
    r_loc = R_EMEtoEMO * r_eci;
    v_loc = R_EMEtoEMO * v_eci;
    if nargout > 2
        dr_loc_dlat = R_EMEtoEMO * dr_eci_dlat;
        dr_loc_dlon = R_EMEtoEMO * dr_eci_dlon;
        dv_loc_dlat = R_EMEtoEMO * dv_eci_dlat;
        dv_loc_dlon = R_EMEtoEMO * dv_eci_dlon;
    end
end