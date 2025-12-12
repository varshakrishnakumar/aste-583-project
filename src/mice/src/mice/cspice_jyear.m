function [jyear] = cspice_jyear
   switch nargin
      case 0
         ;
      otherwise
         error ( 'Usage: [jyear] = cspice_jyear' )
   end
   try
      [jyear] =  mice('jyear_c' );
   catch spiceerr
      rethrow(spiceerr)
   end