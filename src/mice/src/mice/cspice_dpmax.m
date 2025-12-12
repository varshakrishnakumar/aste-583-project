function [dpmax] = cspice_dpmax
   switch nargin
      case 0
         ;
      otherwise
         error ( 'Usage: [dpmax] = cspice_dpmax' )
   end
   try
      [dpmax] = mice('dpmax_c');
   catch spiceerr
      rethrow(spiceerr)
   end