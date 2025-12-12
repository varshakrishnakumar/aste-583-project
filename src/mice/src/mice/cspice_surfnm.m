function [normal] = cspice_surfnm(a, b, c, point)
   switch nargin
      case 4
         a     = zzmice_dp(a);
         b     = zzmice_dp(b);
         c     = zzmice_dp(c);
         point = zzmice_dp(point);
      otherwise
         error ( ['Usage: [_normal(3)_] = '             ...
                  'cspice_surfnm( a, b, c, _point(3)_ )'] )
   end
   try
      [normal] = mice( 'surfnm_c',  a, b, c, point);
   catch spiceerr
      rethrow(spiceerr)
   end