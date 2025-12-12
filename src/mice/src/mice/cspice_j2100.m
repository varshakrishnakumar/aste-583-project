function [j2100] = cspice_j2100
   switch nargin
      case 0
         ;
      otherwise
         error ( 'Usage: [j2100] = cspice_j21000' )
   end
   try
      [j2100] =  mice('j2100_c');
   catch spiceerr
      rethrow(spiceerr)
   end