function [twopi] = cspice_twopi
   switch nargin
      case 0
         ;
      otherwise
         error ( 'Usage: [twopi] = cspice_twopi' )
   end
   try
      [twopi] =  mice('twopi_c');
   catch spiceerr
      rethrow(spiceerr)
   end