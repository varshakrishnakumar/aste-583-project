function [j1950] = cspice_j1950
   switch nargin
      case 0
         ;
      otherwise
         error ( 'Usage: [j1950] = cspice_j1950' )
   end
   try
      [j1950] =  mice('j1950_c');
   catch spiceerr
      rethrow(spiceerr)
   end