function cspice_spkuef( handle )
   switch nargin
      case 1
         handle = zzmice_int(handle);
      otherwise
         error ( 'Usage: cspice_spkuef( handle )' )
   end
   try
      mice('spkuef_c', handle);
   catch spiceerr
      rethrow(spiceerr)
   end