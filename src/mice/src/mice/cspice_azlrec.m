function [rectan] = cspice_azlrec( range, az, el, azccw, elplsz )
   switch nargin
      case 5
         range  = zzmice_dp(range);
         az     = zzmice_dp(az);
         el     = zzmice_dp(el);
         azccw  = zzmice_int(azccw);
         elplsz = zzmice_int(elplsz);
      otherwise
         error ( [ 'Usage: [rectan(3)] = '                                  ...
                   'cspice_azlrec( range, az, el, azccw, elplsz )' ] )
   end
   try
      [rectan] = mice('azlrec_c', range, az, el, azccw, elplsz);
   catch spiceerr
      rethrow(spiceerr)
   end