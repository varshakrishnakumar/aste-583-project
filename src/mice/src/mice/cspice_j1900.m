function [j1900] = cspice_j1900
   switch nargin
      case 0
         ;
      otherwise
         error ( 'Usage: [j1900] = cspice_j1900' )
   end
   try
      [j1900] =  mice('j1900_c');
   catch spiceerr
      rethrow(spiceerr)
   end