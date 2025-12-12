function [b1950] = cspice_b1950
   switch nargin
      case 0
         ;
      otherwise
         error ( 'Usage: [b1950] = cspice_b1950' )
   end
   try
      b1950 =  mice('b1950_c');
   catch spiceerr
      rethrow(spiceerr)
   end