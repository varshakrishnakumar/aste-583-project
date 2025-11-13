function dXdt = twoBodySTM(~, X, mu)
% X = [r(3); v(3); vec(Phi(6x6))]

r   = X(1:3);
v   = X(4:6);
Phi = reshape(X(7:end),6,6);

rn  = norm(r);
a   = -mu * r / rn^3;

% System matrix A = [ 0 I ; dadr  0 ]
I3 = eye(3); Z3 = zeros(3);
% gravity-gradient (∂a/∂r)
mu_over_r5 = mu / rn^5;
dadr = mu_over_r5 * (3*(r*r.') - (rn^2)*eye(3));   % equivalent to -mu/r^3*(I - 3 r r^T/r^2)

A = [Z3, I3;
     dadr, Z3];

% STM differential
dPhidt = A * Phi;

% pack
dXdt = [v; a; reshape(dPhidt,36,1)];
end
