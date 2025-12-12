function [dskdsc] = cspice_dskgd( handle, dladsc )
   switch nargin
      case 2
         handle = zzmice_int( handle );
         dladsc = zzmice_int( dladsc );
      otherwise
         error ( 'Usage: [dskdsc(24)] = cspice_dskgd( handle, dladsc )' )
   end
   try
      [dskdsc] = mice('dskgd_c', handle, dladsc );
   catch spiceerr
      rethrow(spiceerr)
   end