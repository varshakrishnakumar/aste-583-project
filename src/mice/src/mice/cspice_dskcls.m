function cspice_dskcls( handle, optmiz )
   switch nargin
      case 2
         handle = zzmice_int( handle );
         optmiz = zzmice_int( optmiz );
      otherwise
         error ( 'Usage: cspice_dskcls( handle, optmiz )' )
   end
   try
      mice( 'dskcls_c', handle, optmiz);
   catch spiceerr
      rethrow(spiceerr)
   end