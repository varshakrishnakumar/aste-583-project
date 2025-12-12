function [j2000] = cspice_j2000
   switch nargin
      case 0
         ;
      otherwise
         error ( 'Usage: [j2000] = cspice_j2000' )
   end
   try
      [j2000] =  mice('j2000_c' );
   catch spiceerr
      rethrow(spiceerr)
   end