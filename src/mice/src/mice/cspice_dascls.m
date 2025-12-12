function cspice_dascls( handle )
   switch nargin
      case 1
         handle  = zzmice_int(handle);
      otherwise
         error ( 'Usage: cspice_dascls( handle )' )
   end
   try
      mice( 'dascls_c', handle );
   catch spiceerr
      rethrow(spiceerr)
   end