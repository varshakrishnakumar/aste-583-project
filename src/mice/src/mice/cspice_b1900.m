function [b1900] = cspice_b1900
   switch nargin
      case 0
         ;
      otherwise
         error ( 'Usage: [b1900] = cspice_b1900' )
   end
   try
      b1900 =  mice('b1900_c');
   catch spiceerr
      rethrow(spiceerr)
   end