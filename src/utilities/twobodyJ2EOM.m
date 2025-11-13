% ~~~~~~~~~~~~~~~~~~~~~~~~
function dXdt = twobodyJ2EOM(t, X, mu, R_E, J2)
% ------------------------
%{
  This function calculates first and second time derivatives of r
  governed by the equation of two-body 3D motion + J2 perturbation

  Inputs:
    t   - time (s)
    X   - 6x1 state vector [r; v] (km, km/s)_
    mu  - gravitational parameter (km^3/s^2)
    R_E - Earth equatorial radius (km)
    J2  - zonal harmonic coefficient

  Outputs:
    dXdt - 6x1 time derivative [dr/dt; dv/dt]
%}
% ~~~~~~~~~~~~~~~~~~~~~~~~

% State vector
rvec = X(1:3);
vvec = X(4:6);

% Distance
x = rvec(1); y = rvec(2); z = rvec(3);
r = norm(rvec);

% Two-body acceleration
a_tb = -mu/r^3 * rvec;

% J2 perturbation acceleration
z2 = z^2; 
r2 = r^2;
factor = 1.5 * J2 * mu * R_E^2 / r^5;

ax = factor * x * (5*z2/r2 - 1);
ay = factor * y * (5*z2/r2 - 1);
az = factor * z * (5*z2/r2 - 3);

a_J2 = [ax; ay; az];

% Total acceleration
vdotvec = a_tb + a_J2;

% State derivative
rdotvec = vvec;
dXdt = [rdotvec; vdotvec];

end
