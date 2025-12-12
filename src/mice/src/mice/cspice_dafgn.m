function [name] = cspice_dafgn
   switch nargin
      case 0
         ;
      otherwise
         error ( 'Usage: [`name`] = cspice_dafgn' )
   end
   try
      [name] = mice( 'dafgn_c' );
   catch spiceerr
      rethrow(spiceerr)
   end