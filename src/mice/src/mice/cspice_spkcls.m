function cspice_spkcls( handle)
   switch nargin
      case 1
         handle = zzmice_int( handle );
      otherwise
         error ( 'Usage: cspice_spkcls(handle)' )
   end
   try
      mice( 'spkcls_c', handle);
   catch spiceerr
      rethrow(spiceerr)
   end