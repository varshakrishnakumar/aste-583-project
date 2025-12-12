function cspice_dafcls( handle )
   switch nargin
      case 1
         handle  = zzmice_int(handle);
      otherwise
         error ( 'Usage: cspice_dafcls(handle)' )
   end
   try
      mice( 'dafcls_c', handle );
   catch spiceerr
      rethrow(spiceerr)
   end