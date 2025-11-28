function xdot = projectDynamics(T,X,params,sc)
    %Position and Velocity Vecotrs/Magnitudes
    r = X(1:3);
    rmag = norm(r);
    v = X(4:6);

    %STM
    STM = reshape(X(11:end),10,10);

    %SRP Values
    P_SRP = params.psrp; %kN/km^2
    k_SRP = X(7); %unitless

    % mu parameters
    mus = params.mus; %km^3/s^2
    mue = params.mue; %km^3/s^2
    AU = params.AU; % km^2

    % sc parameters
    m = sc.m;
    rho = sc.rho;
    A = sc.A;

    % Earth Vector
    xEarth = cspice_spkezr('Earth',T, 'ECLIPJ2000', 'NONE','Sun');
    r_earth = xEarth(1:3);
    r_earthmag = norm(r_earth);

    % Relative Spacecraft/Earth Vector
    r32 = r - r_earth;
    r32mag = norm(r32);
    
    % Jacobian Matrix
    Amat = zeros(10,10); % Set to All Zeros
    Amat(1:3,4:6) = eye(3); % Identity in Velocity Partials
    G = -1*((mus/rmag^3)*eye(3) - ((A*P_SRP*k_SRP*(1+rho)*AU^2)/(m*rmag^3))*eye(3)... 
        + ((3*A*P_SRP*k_SRP*(1+rho)*AU^2)/(m*rmag^5))*r*transpose(r) - ((3*mus)/rmag^5)*r*transpose(r)...
        + (mue/r32mag^3)*eye(3) - ((3*mue)/r32mag^5)*r32*transpose(r32));
    Amat(4:6,1:3) = G;
    S = ((P_SRP*A*(1+rho)*AU^2)/m)*(r/rmag^3);
    Amat(4:6,7) = S;

    %STM
    STMdot = Amat*STM;
    
    % Dynamics
    f = (-mus/rmag^3)*r - mue*(r32/r32mag^3 + r_earth/r_earthmag^3) + ((k_SRP*P_SRP*(1+rho)*A*AU^2)/m)*(r/rmag^3);

    xdot = [v;f;0;0;0;0;reshape(STMdot,100,1)];
end