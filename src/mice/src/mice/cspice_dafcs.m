function cspice_dafcs( handle )
   switch nargin
      case 1
         handle  = zzmice_int(handle);
      otherwise
         error ( 'Usage: cspice_dafcs(handle)' )
   end
   try
      mice( 'dafcs_c', handle );
   catch spiceerr
      rethrow(spiceerr)
   end