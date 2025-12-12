function [found] = cspice_daffpa
   switch nargin
      case 0
         ;
      otherwise
         error ( 'Usage: [found] = cspice_daffpa' )
   end
   try
      [found] = mice( 'daffpa_c' );
      found   = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end