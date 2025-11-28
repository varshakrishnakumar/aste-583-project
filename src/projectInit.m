function [params,sc,st,X0,P0] = projectInit()
% ----------------------------------------------------------------------%
% This function declares the constants, initial conditions, and all other
% starting information which defines the ASTE583 final project.
%
% Inputs: None
% Outputs:  params = Struct of General Constants
%           sc = Struct of Spacecraft Constants
%           st = Struct of Observation Stations
%           X0 = Initial State Vector
%           P0 = Initial Covariance Matrix
% ----------------------------------------------------------------------%
% General Constants
    %Earth
    params.mue = 398600; %km^3/s^2
    params.Re = 6378; %km
    params.J2 = 0.0010826;
    params.we = 7.292116*10^-5; %rad/s
    params.AU = 149597870.7; %AU in km
    
    %Sun
    params.mus = 132712000000; %km^3/s^2

    %SRP
    params.psrp = (10^-3)*(4.54*10^-6)/(10^-6); %kN/km^2
    params.ksrp = 1; %SRP scale factor - UPDATE BASED ON MEASUREMENTS
    params.sigsrp = 1/3; %SRP a priori 3-sigma value

    %Time
    params.tday = 86400; %s
    params.tyear = 365.25*params.tday; %s
    params.GST0 = (10*60 + 43) * params.we; % GST0 in rad

    %Rotation Matrix EME2000 to EMO2000
    params.eclipticTilt = (23 + 26/60 + 21.448/3600) * (pi/180); %rad
    params.R_EMEtoEMO = [1 0 0; 0 cos(params.eclipticTilt) sin(params.eclipticTilt); 0 -sin(params.eclipticTilt) cos(params.eclipticTilt)];

    %Requirements
    params.maxDeltaV = 1.139; %km/s
    params.sigmaTransverse = 100; %km

% Spacecraft Constants
    sc.m = 200; %kg
    sc.A = 1.5*10^-6; %km^2
    sc.rho = 0.1; %reflectivity coefficient

% Observation Station Constants
    % DSN stations are indices 1-3. Antarctica station is index 4.
    % DSN station locations are truth. Antarctica station is secret, so it
    % is approximate.
    st.lat = [35.244352 -35.220919 40.241355 -80] * (pi/180); %Latitudes in rad
    st.long = [-116.889538 148.981267 -4.2480085 0] * (pi/180); %Longitudes in rad (+ is East)

% Initial State Vector (EMO2000, Sun Centered)
    % State Vector to begin will include R0, V0, ksrp. station 4 latitude,
    % station 4 longitude, and range rate bias (initially set to zero).
    R0 = [1.067623147085261; 1.148757045773147; -0.000321627221208]*10^8; %km
    V0 = [-22.148505873534173; 18.814312217049999; -0.098774507382220]; %km/s
    X0 = [R0;V0; params.ksrp; st.lat(4); st.long(4); 0];

% Initial Covariance Matrix
    P0 = zeros(6,6);
    for(i = 1:1:3)
        P0(i,i) = 100^2; %100km spherical position covariance
    end
    for(i = 4:1:6)
        P0(i,i) = 0.001^2; %1m/s spherical velocity covariance
    end 


end