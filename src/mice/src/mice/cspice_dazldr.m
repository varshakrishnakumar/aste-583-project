function [jacobi] = cspice_dazldr( x, y, z, azccw, elplsz )
   switch nargin
      case 5
         x      = zzmice_dp(x);
         y      = zzmice_dp(y);
         z      = zzmice_dp(z);
         azccw  = zzmice_int(azccw);
         elplsz = zzmice_int(elplsz);
      otherwise
         error ( [ 'Usage: [jacobi(3,3)] = '                                ...
                   'cspice_dazldr( x, y, z, azccw, elplsz )' ] )
   end
   try
      [jacobi] = mice('dazldr_c', x, y, z, azccw, elplsz);
   catch spiceerr
      rethrow(spiceerr)
   end