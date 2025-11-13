function [rI, vI] = oe2rv(a,e,i,RAAN,w,nu,mu)
% classical elements to inertial r,v (km, km/s)
p = a*(1 - e^2);
rpf = p / (1 + e*cos(nu));
rPQW = [ rpf*cos(nu); rpf*sin(nu); 0 ];
vPQW = [-sqrt(mu/p)*sin(nu);
         sqrt(mu/p)*(e + cos(nu));
         0];

% rotation: PQW -> IJK
cO=cos(RAAN); sO=sin(RAAN);
ci=cos(i);    si=sin(i);
cw=cos(w);    sw=sin(w);

RzO = [ cO -sO 0; sO cO 0; 0 0 1 ];
Rxi = [ 1 0 0; 0 ci -si; 0 si ci ];
Rzw = [ cw -sw 0; sw cw 0; 0 0 1 ];
Q = RzO*Rxi*Rzw;              

rI = Q*rPQW;
vI = Q*vPQW;
end
