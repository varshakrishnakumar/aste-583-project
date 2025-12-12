function [onepi] = cspice_pi
   switch nargin
      case 0
         ;
      otherwise
         error ( 'Usage: [onepi] = cspice_pi' )
   end
   try
      [onepi] =  mice('pi_c');
   catch spiceerr
      rethrow(spiceerr)
   end