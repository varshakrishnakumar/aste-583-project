function [data] = cspice_dasrdd( handle, first, last )
   switch nargin
      case 3
         handle = zzmice_int(handle);
         first  = zzmice_int(first);
         last   = zzmice_int(last);
      otherwise
         error ( 'Usage: [data] = cspice_dasrdd( handle, first, last )' )
   end
   try
      [data] = mice('dasrdd_c', handle, first, last);
   catch spiceerr
      rethrow(spiceerr)
   end