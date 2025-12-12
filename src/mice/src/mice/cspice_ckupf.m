function cspice_ckupf( handle )
   switch nargin
      case 1
         handle = zzmice_int(handle);
      otherwise
         error ( 'Usage: cspice_ckupf( handle )' )
   end
   try
      mice('ckupf_c', handle);
   catch spiceerr
      rethrow(spiceerr)
   end