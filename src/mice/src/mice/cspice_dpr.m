function [dpr] = cspice_dpr
   switch nargin
      case 0
         ;
      otherwise
         error ( 'Usage: [dpr] = cspice_dpr' )
   end
   try
      [dpr] =  mice('dpr_c');
   catch spiceerr
      rethrow(spiceerr)
   end