function [handle] = cspice_dafopr( fname )
   switch nargin
      case 1
         fname  = zzmice_str(fname);
      otherwise
         error ( 'Usage: [handle] = cspice_dafopr(`fname`)' )
   end
   try
      [handle] = mice( 'dafopr_c', fname );
   catch spiceerr
      rethrow(spiceerr)
   end