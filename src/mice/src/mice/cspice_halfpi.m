function [halfpi] = cspice_halfpi
   switch nargin
      case 0
         ;
      otherwise
         error ( 'Usage: [halfpi] = cspice_halfpi' )
   end
   try
      [halfpi] =  mice('halfpi_c');
   catch spiceerr
      rethrow(spiceerr)
   end