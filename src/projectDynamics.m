function xdot = projectDynamics(t, X, params, sc)
%PROJECTDYNAMICS Sun-centered dynamics (+ optional STM) for ASTE583 project.
%
% State (10x1):
%   X(1:3)   r_sc_sun   [km]    Spacecraft position wrt Sun (ecliptic J2000)
%   X(4:6)   v_sc_sun   [km/s]  Spacecraft velocity wrt Sun
%   X(7)     k_srp      [-]     SRP scale factor
%   X(8)     lat_stn4   [rad]   (constant in dynamics)
%   X(9)     lon_stn4   [rad]   (constant in dynamics)
%   X(10)    rr_bias    [km/s]  (constant in dynamics)
%
% Optional (appended):
%   X(11:end) = Phi(:) where Phi is 10x10 STM (column-stacked).
%
% Time:
%   t = ephemeris time seconds past J2000 (ET), compatible with SPICE.
%
% Output:
%   xdot has the same length as X.

    X = X(:); % ensure column vector

    n = numel(X);
    hasSTM = (n == 110);
    if ~(n == 10 || hasSTM)
        error('projectDynamics:BadStateLength', ...
            'Expected X length 10 (state) or 110 (state+STM). Got %d.', n);
    end

    % -------------------- Unpack state -------------------- %
    r     = X(1:3);
    v     = X(4:6);
    k_srp = X(7);

    rmag = norm(r);

    % -------------------- Ephemeris (Earth wrt Sun) -------------------- %
    % NOTE: SPICE calls inside an ODE RHS can get slow for long propagations/
    % Monte Carlo. If runtime becomes painful, precompute Earth ephemeris and
    % interpolate instead.
    xEarth = cspice_spkezr('EARTH', t, 'ECLIPJ2000', 'NONE', 'SUN'); % km, km/s
    r_e    = xEarth(1:3);
    r_emag = norm(r_e);

    r_sc_e    = r - r_e;          % spacecraft wrt Earth (in same inertial frame)
    r_sc_emag = norm(r_sc_e);

    % -------------------- Accelerations (km/s^2) -------------------- %
    mus = params.mus;
    mue = params.mue;

    % Sun gravity
    a_sun = -mus * r / rmag^3;

    % Earth third-body (includes indirect term so a_earth = 0 when r = 0)
    a_earth = -mue * ( r_sc_e / r_sc_emag^3 + r_e / r_emag^3 );

    % SRP (no shadowing)
    % params.psrp is in kN/km^2; with A in km^2 and m in kg, this yields km/s^2
    P0 = params.psrp;       % kN/km^2 at 1 AU
    AU = params.AU;         % km
    Cr = 1 + sc.rho;        % unitless (simple reflectivity model)

    C0 = (P0 * sc.A * Cr * AU^2) / sc.m;  % km/s^2 * (AU/r)^2 handled via r/r^3
    a_srp = k_srp * C0 * r / rmag^3;

    a = a_sun + a_earth + a_srp;

    % 10-state derivative (augmented params are constants in the dynamics)
    xdot_state = [v;
                  a;
                  0; 0; 0; 0];

    % -------------------- Optional STM propagation -------------------- %
    if hasSTM
        Phi = reshape(X(11:end), 10, 10);

        I3 = eye(3);

        % Helper: d/dr ( r / ||r||^3 ) = I/r^3 - 3 rr^T / r^5
        rrT = r * r.';
        r3  = rmag^3;
        r5  = rmag^5;
        D_r = (I3 / r3) - (3 * rrT / r5);

        % Gradients wrt position
        G_sun  = -mus * D_r;
        G_srp  = (k_srp * C0) * D_r;

        r32    = r_sc_e;
        r32mag = r_sc_emag;
        r32r32T = r32 * r32.';
        r32_3  = r32mag^3;
        r32_5  = r32mag^5;
        D_r32  = (I3 / r32_3) - (3 * r32r32T / r32_5);
        G_earth = -mue * D_r32;   % indirect term has no r-dependence

        % Assemble A = df/dx (10x10)
        A = zeros(10,10);
        A(1:3,4:6) = I3;
        A(4:6,1:3) = G_sun + G_earth + G_srp;
        A(4:6,7)   = C0 * r / rmag^3; % da/dk_srp

        Phi_dot = A * Phi;

        xdot = [xdot_state;
                Phi_dot(:)];
    else
        xdot = xdot_state;
    end
end
