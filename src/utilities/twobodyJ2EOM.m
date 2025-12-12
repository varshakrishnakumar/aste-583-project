function dXdt = twobodyJ2EOM(t, X, mu, R_E, J2)
rvec = X(1:3);
vvec = X(4:6);
x = rvec(1); y = rvec(2); z = rvec(3);
r = norm(rvec);
a_tb = -mu/r^3 * rvec;
z2 = z^2;
r2 = r^2;
factor = 1.5 * J2 * mu * R_E^2 / r^5;
ax = factor * x * (5*z2/r2 - 1);
ay = factor * y * (5*z2/r2 - 1);
az = factor * z * (5*z2/r2 - 3);
a_J2 = [ax; ay; az];
vdotvec = a_tb + a_J2;
rdotvec = vvec;
dXdt = [rdotvec; vdotvec];
end