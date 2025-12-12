function [r] = cspice_vrotv( v, axis, theta)
   switch nargin
      case 3
         v     = zzmice_dp(v);
         axis  = zzmice_dp(axis);
         theta = zzmice_dp(theta);
      otherwise
         error ( ['Usage: [r(3)] = cspice_vrotv( v(3), axis(3), theta)'] )
   end
   try
      [r] = mice('vrotv_c', v, axis, theta);
   catch spiceerr
      rethrow(spiceerr)
   end