function cspice_dafbbs( handle )
   switch nargin
      case 1
         handle  = zzmice_int(handle);
      otherwise
         error ( 'Usage: cspice_dafbbs(handle)' )
   end
   try
      mice( 'dafbbs_c', handle );
   catch spiceerr
      rethrow(spiceerr)
   end