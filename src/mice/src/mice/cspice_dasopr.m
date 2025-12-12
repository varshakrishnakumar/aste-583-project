function [handle] = cspice_dasopr( fname )
   switch nargin
      case 1
         fname  = zzmice_str(fname);
      otherwise
         error ( 'Usage: [handle] = cspice_dafopr(`fname`)' )
   end
   try
      [handle] = mice( 'dasopr_c', fname );
   catch spiceerr
      rethrow(spiceerr)
   end