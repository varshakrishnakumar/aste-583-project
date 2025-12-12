function [r, clon, z] = cspice_reccyl(rectan)
   switch nargin
      case 1
         rectan = zzmice_dp(rectan);
      otherwise
         error ( 'Usage: [_r_, _clon_, _z_] = cspice_reccyl(_rectan(3)_)' )
   end
   try
      [r, clon, z] = mice('reccyl_c',rectan);
   catch spiceerr
      rethrow(spiceerr)
   end