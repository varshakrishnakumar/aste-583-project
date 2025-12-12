function [tyear] = cspice_tyear
   switch nargin
      case 0
         ;
      otherwise
         error ( 'Usage: [tyear] = cspice_tyear' )
   end
   try
      [tyear] =  mice('tyear_c' );
   catch spiceerr
      rethrow(spiceerr)
   end