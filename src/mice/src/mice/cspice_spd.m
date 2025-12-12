function [spd] = cspice_spd
   switch nargin
      case 0
         ;
      otherwise
         error ( 'Usage: [spd] = cspice_spd' )
   end
   try
      [spd] =  mice('spd_c');
   catch spiceerr
      rethrow(spiceerr)
   end