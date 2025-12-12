function [found] = cspice_daffna
   switch nargin
      case 0
         ;
      otherwise
         error ( 'Usage: [found] = cspice_daffna' )
   end
   try
      [found] = mice( 'daffna_c' );
      found   = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end