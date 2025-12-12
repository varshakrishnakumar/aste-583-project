function [clight] = cspice_clight
   switch nargin
      case 0
         ;
      otherwise
         error ( 'Usage: [clight] = cspice_clight' )
   end
   try
      [clight] =  mice('clight_c');
   catch spiceerr
      rethrow(spiceerr)
   end