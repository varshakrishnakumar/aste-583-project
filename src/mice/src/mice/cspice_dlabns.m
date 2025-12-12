function cspice_dlabns( handle )
   switch nargin
      case 1
         handle = zzmice_int(handle);
      otherwise
         error ( 'Usage: cspice_dlabns( handle )' )
   end
   try
      mice('dlabns_c', handle);
   catch spiceerr
      rethrow(spiceerr)
   end