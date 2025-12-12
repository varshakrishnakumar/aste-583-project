function [intmax] = cspice_intmax
   switch nargin
      case 0
         ;
      otherwise
         error ( 'Usage: [intmax] = cspice_intmax' )
   end
   try
      [intmax] = mice('intmax_c');
   catch spiceerr
      rethrow(spiceerr)
   end