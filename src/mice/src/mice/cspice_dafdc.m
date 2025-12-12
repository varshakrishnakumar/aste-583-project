function cspice_dafdc( handle )
   switch nargin
      case 1
         handle  = zzmice_int(handle);
      otherwise
         error ( 'Usage: cspice_dafdc(handle)' )
   end
   try
      mice( 'dafdc_c', handle );
   catch spiceerr
      rethrow(spiceerr)
   end