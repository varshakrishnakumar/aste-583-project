function cspice_daswbr( handle )
   switch nargin
      case 1
         handle = zzmice_int(handle);
      otherwise
         error ( 'Usage: cspice_daswbr( handle )' )
   end
   try
      mice('daswbr_c', handle);
   catch spiceerr
      rethrow(spiceerr)
   end