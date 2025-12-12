function [rpd] = cspice_rpd
   switch nargin
      case 0
         ;
      otherwise
         error ( 'Usage: [rpd] = cspice_rpd' )
   end
   try
      [rpd] =  mice('rpd_c');
   catch spiceerr
      rethrow(spiceerr)
   end