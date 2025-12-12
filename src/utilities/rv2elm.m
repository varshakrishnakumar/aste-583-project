function [elm, vec] = rv2elm(r, v, mu)
    if nargin < 3, mu = 398600.4418; end
    r = r(:); v = v(:);
    R = norm(r); V = norm(v);
    h = cross(r,v); H = norm(h);
    n = cross([0;0;1],h); N = norm(n);
    evec = ((V^2 - mu/R)*r - dot(r,v)*v) / mu;
    e = norm(evec);
    energy = V^2/2 - mu/R;
    if abs(e-1) < 1e-8
        a = Inf;
    else
        a = -mu / (2*energy);
    end
    i = acos(h(3)/H);
    Om = atan2(n(2), n(1));
    if Om < 0, Om = Om + 2*pi; end
    if evec(3) >= 0
        w = acos(dot(n, evec)/(N * e));
    else
        w = 2*pi - acos(dot(n, evec)/(N * e));
    end
    if e > 1e-8
        sinf = dot(h, cross(evec, r)) / (e*H*R);
        cosf = dot(evec, r) / (e*R);
        f = atan2(sinf, cosf);
    else
        sinu = dot(h, cross(n, r)) / (N*H*R);
        cosu = dot(n, r) / (N*R);
        f = atan2(sinu, cosu);
    end
    if f < 0, f = f + 2*pi; end
    if e < 1 - 1e-8
        cosE = (e + cos(f)) / (1 + e*cos(f));
        sinE = sqrt(1 - e^2) * sin(f) / (1 + e*cos(f));
        E = atan2(sinE, cosE);
        if E < 0, E = E + 2*pi; end
        M = mod(E - e*sin(E), 2*pi);
    else
        E = NaN;
        M = NaN;
    end
    elm = struct('a',a,'e',e,'i',i,'Om',Om,'w',w,'f',f,'E',E,'M',M);
    vec = [a e i Om w f E M];
end