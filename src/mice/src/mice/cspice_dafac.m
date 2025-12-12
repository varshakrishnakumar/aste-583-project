function cspice_dafac( handle, buffer )
   switch nargin
      case 2
         handle  = zzmice_int(handle);
         buffer  = zzmice_str(buffer);
      otherwise
         error ( 'Usage: cspice_dafac( handle, buffer )' )
   end
   try
      mice( 'dafac_c', handle, buffer );
   catch spiceerr
      rethrow(spiceerr)
   end