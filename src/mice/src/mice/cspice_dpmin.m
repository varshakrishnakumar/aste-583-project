function [dpmin] = cspice_dpmin
   switch nargin
      case 0
         ;
      otherwise
         error ( 'Usage: [dpmin] = cspice_dpmin' )
   end
   try
      [dpmin] = mice('dpmin_c');
   catch spiceerr
      rethrow(spiceerr)
   end