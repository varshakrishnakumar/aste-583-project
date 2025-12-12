function [handle] = cspice_dasopw( fname )
   switch nargin
      case 1
         fname = zzmice_str(fname);
      otherwise
         error ( 'Usage: [handle] = cspice_dasopw( `fname` )' )
   end
   try
      [handle] = mice('dasopw_c', fname);
   catch spiceerr
      rethrow(spiceerr)
   end