function [jacobi] = cspice_drdazl( range, az, el, azccw, elplsz )
   switch nargin
      case 5
         range  = zzmice_dp(range);
         az     = zzmice_dp(az);
         el     = zzmice_dp(el);
         azccw  = zzmice_int(azccw);
         elplsz = zzmice_int(elplsz);
      otherwise
         error ( [ 'Usage: [jacobi(3,3)] = '                                ...
                   'cspice_drdazl( range, az, el, azccw, elplsz )' ] )
   end
   try
      [jacobi] = mice('drdazl_c', range, az, el, azccw, elplsz);
   catch spiceerr
      rethrow(spiceerr)
   end