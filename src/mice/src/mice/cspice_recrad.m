function [range, ra, dec] = cspice_recrad(rectan)
   switch nargin
      case 1
         rectan = zzmice_dp(rectan);
      otherwise
         error ( 'Usage: [_range_, _ra_, _dec_] = cspice_recrad(_rectan(3)_)' )
   end
   try
      [range, ra, dec] = mice('recrad_c',rectan);
   catch spiceerr
      rethrow(spiceerr)
   end