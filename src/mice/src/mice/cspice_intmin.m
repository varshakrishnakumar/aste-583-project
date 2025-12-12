function [intmin] = cspice_intmin
   switch nargin
      case 0
         ;
      otherwise
         error ( 'Usage: [intmin] = cspice_intmin' )
   end
   try
      [intmin] = mice('intmin_c');
   catch spiceerr
      rethrow(spiceerr)
   end