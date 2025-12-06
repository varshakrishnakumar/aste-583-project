function [params,sc,st,X0,P0] = projectInit()
% ----------------------------------------------------------------------%
% This function declares the constants, initial conditions, and all other
% starting information which defines the ASTE583 final project.
%
% Inputs: None
% Outputs:  params = Struct of General Constants
%           sc = Struct of Spacecraft Constants
%           st = Struct of Observation Stations
%           X0 = Initial State Vector (10x1)
%           P0 = Initial Covariance Matrix (10x10)
% ----------------------------------------------------------------------%

% -------------------- General Constants -------------------- %
    % Earth
    params.mue  = 398600;          % km^3/s^2
    params.Re   = 6378;            % km
    params.J2   = 0.0010826;
    params.we   = 7.292116e-5;     % rad/s
    params.AU   = 149597870.7;     % km

    % Sun
    params.mus  = 132712000000;    % km^3/s^2

    % SRP
    params.psrp   = (1e-3)*(4.54e-6)/(1e-6); % kN/km^2
    params.ksrp   = 1;          % SRP scale factor (state 7 initial guess)
    params.sigsrp = 1/3;        % SRP a priori 3-sigma value (unitless)

    % Time
    params.tday  = 86400;                      % s
    params.tyear = 365.25*params.tday;         % s
    params.GST0  = (10*60 + 43) * params.we;   % rad

    % Rotation Matrix EME2000 to EMO2000 (equatorial -> ecliptic)
    params.eclipticTilt = (23 + 26/60 + 21.448/3600) * (pi/180); % rad
    cT = cos(params.eclipticTilt);
    sT = sin(params.eclipticTilt);
    params.R_EMEtoEMO = [ 1  0   0;
                          0  cT  sT;
                          0 -sT  cT ];

    % Requirements
    params.maxDeltaV      = 1.139;  % km/s
    params.sigmaTransverse = 100;   % km

% -------------------- Spacecraft Constants -------------------- %
    sc.m   = 200;        % kg
    sc.A   = 1.5e-6;     % km^2
    sc.rho = 0.1;        % reflectivity coefficient

% -------------------- Observation Station Constants -------------------- %
    % DSN stations are indices 1-3. Antarctica station is index 4.
    % DSN station locations are truth. Antarctica station is secret, so it
    % is approximate.
    st.lat  = [ 35.244352,  -35.220919,  40.241355,  -80 ] * (pi/180); % rad
    st.long = [-116.889538, 148.981267,  -4.2480085,   0 ] * (pi/180); % rad (+ East)

% -------------------- Initial State Vector (Sun-centered EMO2000) ------- %
    % State: [R0; V0; k_SRP; lat_stn4; lon_stn4; rangeRateBias]
    R0 = [ 1.067623147085261;
           1.148757045773147;
          -0.000321627221208 ] * 1e8;    % km

    V0 = [ -22.148505873534173;
            18.814312217049999;
            -0.098774507382220 ];        % km/s

    X0 = [ R0;
           V0;
           params.ksrp;   % k_SRP
           st.lat(4);     % Antarctica lat
           st.long(4);    % Antarctica lon
           0 ];           % range-rate bias (km/s), initial 0

% -------------------- Initial Covariance Matrix (10x10) ----------------- %
    % Position/velocity 3-sigma:
    %   pos: 100 km (spherical)
    %   vel: 0.001 km/s (1 m/s)
    %
    % Augmented parameters (example 3-sigma values; tune as needed):
    %   k_SRP: params.sigsrp (3-sigma)
    %   lat_4, lon_4: ~5 deg (3-sigma)
    %   bias_rr: 1e-4 km/s (3-sigma) ~ 0.1 m/s

    P0 = zeros(10,10);

    % pos
    for i = 1:3
        P0(i,i) = (100)^2;       % (km)^2
    end

    % vel
    for i = 4:6
        P0(i,i) = (0.001)^2;     % (km/s)^2
    end

    % k_SRP (unitless)
    sig_k = params.sigsrp;      % already 1-sigma
    P0(7,7) = sig_k^2;

    % Antarctica lat/lon (rad)
    sig_lat4 = (5*pi/180) / 3;   % 5 deg 3-sigma -> 1-sigma (example)
    sig_lon4 = (5*pi/180) / 3;
    P0(8,8) = sig_lat4^2;
    P0(9,9) = sig_lon4^2;

    % range-rate bias (km/s)
    sig_bias = (1e-4)/3;         % 3-sigma = 1e-4 km/s (0.1 m/s)
    P0(10,10) = sig_bias^2;
end
