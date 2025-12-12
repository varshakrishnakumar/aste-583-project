function cspice_dasllc( handle )
   switch nargin
      case 1
         handle = zzmice_int(handle);
      otherwise
         error ( 'Usage: cspice_dasllc( handle )' )
   end
   try
      mice('dasllc_c', handle);
   catch spiceerr
      rethrow(spiceerr)
   end