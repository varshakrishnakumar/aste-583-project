function xdot = projectDynamics(t, X, params, sc)
    X = X(:);
    n = numel(X);
    hasSTM = (n == 110);
    if ~(n == 10 || hasSTM)
        error('projectDynamics:BadStateLength', ...
            'Expected X length 10 (state) or 110 (state+STM). Got %d.', n);
    end
    r     = X(1:3);
    v     = X(4:6);
    k_srp = X(7);
    rmag = norm(r);
    xEarth = cspice_spkezr('EARTH', t, 'ECLIPJ2000', 'NONE', 'SUN');
    r_e    = xEarth(1:3);
    r_emag = norm(r_e);
    r_sc_e    = r - r_e;
    r_sc_emag = norm(r_sc_e);
    mus = params.mus;
    mue = params.mue;
    a_sun = -mus * r / rmag^3;
    a_earth = -mue * ( r_sc_e / r_sc_emag^3 + r_e / r_emag^3 );
    P0 = params.psrp;
    AU = params.AU;
    Cr = 1 + sc.rho;
    C0 = (P0 * sc.A * Cr * AU^2) / sc.m;
    a_srp = k_srp * C0 * r / rmag^3;
    a = a_sun + a_earth + a_srp;
    xdot_state = [v;
                  a;
                  0; 0; 0; 0];
    if hasSTM
        Phi = reshape(X(11:end), 10, 10);
        I3 = eye(3);
        rrT = r * r.';
        r3  = rmag^3;
        r5  = rmag^5;
        D_r = (I3 / r3) - (3 * rrT / r5);
        G_sun  = -mus * D_r;
        G_srp  = (k_srp * C0) * D_r;
        r32    = r_sc_e;
        r32mag = r_sc_emag;
        r32r32T = r32 * r32.';
        r32_3  = r32mag^3;
        r32_5  = r32mag^5;
        D_r32  = (I3 / r32_3) - (3 * r32r32T / r32_5);
        G_earth = -mue * D_r32;
        A = zeros(10,10);
        A(1:3,4:6) = I3;
        A(4:6,1:3) = G_sun + G_earth + G_srp;
        A(4:6,7)   = C0 * r / rmag^3;
        Phi_dot = A * Phi;
        xdot = [xdot_state;
                Phi_dot(:)];
    else
        xdot = xdot_state;
    end
end