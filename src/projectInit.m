function [params,sc,st,X0,P0] = projectInit()
    params.mue  = 398600;
    params.Re   = 6378;
    params.J2   = 0.0010826;
    params.we   = 7.292116e-5;
    params.AU   = 149597870.7;
    params.mus  = 132712000000;
    params.psrp   = (1e-3)*(4.54e-6)/(1e-6);
    params.ksrp   = 1;
    params.sigsrp = 1/3;
    params.tday  = 86400;
    params.tyear = 365.25*params.tday;
    params.GST0 = (10/60 + 43/3600) * (pi/180);
    params.eclipticTilt = (23 + 26/60 + 21.448/3600) * (pi/180);
    cT = cos(params.eclipticTilt);
    sT = sin(params.eclipticTilt);
    params.R_EMEtoEMO = [ 1  0   0;
                          0  cT  sT;
                          0 -sT  cT ];
    params.maxDeltaV      = 1.139;
    params.sigmaTransverse = 100;
    sc.m   = 200;
    sc.A   = 1.5e-6;
    sc.rho = 0.1;
    st.lat  = [ 35.244352,  -35.220919,  40.241355,  -80 ] * (pi/180);
    st.long = [-116.889538, 148.981267,  -4.2480085,   0 ] * (pi/180);
    R0 = [ 1.067623147085261;
           1.148757045773147;
          -0.000321627221208 ] * 1e8;
    V0 = [ -22.148505873534173;
            18.814312217049999;
            -0.098774507382220 ];
    X0 = [ R0;
           V0;
           params.ksrp;
           st.lat(4);
           st.long(4);
           0 ];
    P0 = zeros(10,10);
    for i = 1:3
        P0(i,i) = (100)^2;
    end
    for i = 4:6
        P0(i,i) = (0.001)^2;
    end
    sig_k = params.sigsrp;
    P0(7,7) = sig_k^2;
    sig_lat4 = (5*pi/180) / 3;
    sig_lon4 = (5*pi/180) / 3;
    P0(8,8) = sig_lat4^2;
    P0(9,9) = sig_lon4^2;
    sig_bias = (1e-4)/3;
    P0(10,10) = sig_bias^2;
end