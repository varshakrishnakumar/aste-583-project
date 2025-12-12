function [handle] = cspice_dasops
   switch nargin
      case 0
         ;
      otherwise
         error ( 'Usage: [handle] = cspice_dasops' )
   end
   try
      [handle] = mice('dasops_c');
   catch spiceerr
      rethrow(spiceerr)
   end