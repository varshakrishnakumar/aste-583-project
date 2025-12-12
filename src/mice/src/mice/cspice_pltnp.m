function [pnear, dist] = cspice_pltnp(point, v1, v2, v3)
   switch nargin
      case 4
         point = zzmice_dp(point);
         v1    = zzmice_dp(v1);
         v2    = zzmice_dp(v2);
         v3    = zzmice_dp(v3);
      otherwise
         error ( 'Usage: [pnear, dist] = cspice_pltnp(point, v1, v2, v3)' )
   end
   try
      [pnear, dist] = mice( 'pltnp_c', point, v1, v2, v3);
   catch spiceerr
      rethrow(spiceerr)
   end