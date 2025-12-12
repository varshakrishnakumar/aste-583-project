function cspice_ckcls( handle)
   switch nargin
      case 1
         handle = zzmice_int( handle );
      otherwise
         error ( 'Usage: cspice_ckcls(handle)' )
   end
   try
      mice( 'ckcls_c', handle);
   catch spiceerr
      rethrow(spiceerr)
   end