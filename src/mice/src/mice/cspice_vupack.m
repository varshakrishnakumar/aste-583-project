function [x, y, z] = cspice_vupack( v )
   switch nargin
      case 1
         v = zzmice_dp(v);
      otherwise
         error ( 'Usage: [x, y, z] = cspice_vupack( v(3) )' )
   end
   try
      [x, y, z] = mice('vupack_c', v);
   catch spiceerr
      rethrow(spiceerr)
   end