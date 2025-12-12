function [range, az, el] = cspice_recazl( rectan, azccw, elplsz )
   switch nargin
      case 3
         rectan = zzmice_dp(rectan);
         azccw  = zzmice_int(azccw);
         elplsz = zzmice_int(elplsz);
      otherwise
         error ( [ 'Usage: [range, az, el] = '                              ...
                   'cspice_recazl( rectan(3), azccw, elplsz )' ] )
   end
   try
      [range, az, el] = mice('recazl_c', rectan, azccw, elplsz);
   catch spiceerr
      rethrow(spiceerr)
   end