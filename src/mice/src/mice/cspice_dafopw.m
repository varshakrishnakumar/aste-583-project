function [handle] = cspice_dafopw( fname )
   switch nargin
      case 1
         fname  = zzmice_str(fname);
      otherwise
         error ( 'Usage: [handle] = cspice_dafopw(`fname`)' )
   end
   try
      [handle] = mice( 'dafopw_c', fname );
   catch spiceerr
      rethrow(spiceerr)
   end