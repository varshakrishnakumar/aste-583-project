function [handle] = cspice_spklef( fname )
   switch nargin
      case 1
         fname = zzmice_str(fname);
      otherwise
         error ( 'Usage: [handle] = cspice_spklef( `fname` )' )
   end
   try
      [handle] = mice('spklef_c', fname);
   catch spiceerr
      rethrow(spiceerr)
   end