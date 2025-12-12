function cspice_dafbfs( handle )
   switch nargin
      case 1
         handle  = zzmice_int(handle);
      otherwise
         error ( 'Usage: cspice_dafbfs(handle)' )
   end
   try
      mice( 'dafbfs_c', handle );
   catch spiceerr
      rethrow(spiceerr)
   end